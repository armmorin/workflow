from os import environ
from rich.console import Console
import shutil
from pathlib import Path
from typing import Optional, Tuple
from ase import Atoms
from ase.io import read, write
from ase.db import connect
from ase.db.sqlite import SQLite3Database
from forcecurve import fit_images
import numpy as np
from ase.optimize import FIRE, BFGS
from ase.neb import NEB, NEBTools
import matplotlib.pyplot as plt
from herculestools.dft import (
    RunConfiguration as RC,
    create_Vasp_calc,
    set_magnetic_moments,
    create_neb_path
)

here = Path(__file__).resolve().parent.parent
rc_home = RC.home
structures = RC.structures_dir
c = Console()
nnodes = int(environ['SLURM_NNODES'])
ncore = int(environ['NCORE'])

def main(db_id: Optional[int] = None,
         N_images: int = 3, climb: bool = False, fmax: float = 0.05, parallel: bool = False,
         fire: bool = True, vasp: dict = {}, **kwargs) -> Tuple[bool, Optional[dict]]:
    """
    Perform the main NEB (Nudged Elastic Band) calculation.

    Args:
        db_id (Optional[int]): The ID of the NEB in the database. Default is None.
        N_images (int): The number of images in the NEB calculation. Default is 3.
        climb (bool): Whether to use the climbing image NEB method. Default is True.
        fmax (float): The maximum force tolerance for the optimization. Default is 0.05.
        parallel (bool): Whether to run the calculation in parallel. Default is False.
        fire (bool): Whether to use the FIRE optimizer. Default is True.
        vasp (dict): Additional parameters for the VASP calculator. Default is an empty dictionary.

    Returns:
        Tuple[bool, Optional[dict]]: A tuple containing a boolean indicating whether the calculation was successful,
        and an optional dictionary with the trajectory file path and NEB ID.
    """

    # Set parallel execution
    shared_calc = not parallel

    # Database connection
    db_path = structures / "hexag_perovs_wdiscards.db"
    db: SQLite3Database = connect(db_path)

    # Get the name of the structure from the database.
    try:
        iID = kwargs['initial_id']
        c.log(f'Initial ID: {iID}')
        fID = kwargs['final_id']
        c.log(f'Final ID: {fID}')
        vi = kwargs['initial_vac']
        c.log(f'Initial vac: {vi}')
        vf = kwargs['final_vac']
        c.log(f'Final vac: {vf}')
        row = db.get(iID)
        db_name = row.name
        names_list = db_name.split('_')
        name = '_'.join(names_list[:-1])
    except Exception as e:
        c.log(f"Error: {e}")
        # If 'climb' is True, then the db_id is found in the neb_id kwargs
        try:
            db_id = kwargs['neb_id']
        except KeyError:
            db_id = db_id
        if climb:
            c.log(f"Using the climbing image NEB method for db_id: {db_id}")
        else:
            c.log(f"Using the normal NEB method for db_id: {db_id}")
        
        #c.log(f"No initial or final ID provided. Using db_id: {db_id}")
        db_name = db.get(db_id).name
        names_list = db_name.split('_')
        name = '_'.join(names_list[:-1])
        iID = db.get(selection=f"name={name}_vi").id
        fID = db.get(selection=f"name={name}_vf").id
    
    if climb:
        neb_dir = Path(rc_home / 'CINEB' / f"{name}")
        neb_name = f"{name}_cineb"
    else:
        neb_dir = Path(rc_home / 'NEB' / f"{name}")
        neb_name = f"{name}_neb"
        
    neb_dir.mkdir(parents=True, exist_ok=True)
    # The trajectory file must be written only once.
    trajs = []
    path = "traj"
    start_traj = neb_dir / f"{name}_start.{path}"
    
    # If climb is True, check if the directory exists and if the trajectory file doesn't exist.
    try:
        trajs = [i for i in neb_dir.glob(f"*{path}") if i.is_file() and i.stat().st_size > 0]
        traj = max(trajs, key=lambda a: a.stat().st_mtime)
        # Use this trajectory to check for any newer WAVECAR files created after the trajectory's last modification
        traj_mtime = traj.stat().st_mtime
        c.log(f"Found the most recent trajectory file in {traj}.")
        for wavecar in neb_dir.rglob("WAVECAR"):
            if wavecar.stat().st_mtime > traj_mtime:
                c.log(f"Found a newer WAVECAR file in {wavecar.parent}.")
                wavecar.unlink()
            else:
                continue
        images = read(traj.as_posix(),  index=f"-{N_images+2}:")
        c.log(f"trajectory file: {traj.stem} has {len(images)} images.")

        neb = NEB(images, climb=climb, parallel=parallel,
                method="improvedtangent", allow_shared_calculator=shared_calc)
    
    except Exception as e:
            c.log(f"Error: {e}")
            if climb:
                # Copy the trajectory and WAVECAR files, if they exist, to the new directory
                old_neb_dir = Path(rc_home / 'NEB' / name)
                old_trajs = [i for i in old_neb_dir.glob(f"*{path}") if i.is_file() and i.stat().st_size > 0]
                # Sort the trajectories by modification time, and take the most recent one
                old_traj = max(old_trajs, key=lambda a: a.stat().st_mtime)
                # Write the trajectory file to the new directory
                images = read(old_traj.as_posix(), index=f"-{N_images+2}:")
                write(str(start_traj), images)
                c.log(f"Copied the trajectory file from {old_neb_dir} to {neb_dir}")
                # Copy the WAVECAR files and their parent directories to the new directory
                for wavecar in old_neb_dir.rglob("WAVECAR"):
                    new_wavecar = neb_dir / wavecar.relative_to(old_neb_dir)
                    new_wavecar.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy(wavecar, new_wavecar)
                c.log(f"Copied the WAVECAR files from {old_neb_dir} to {neb_dir}")
                neb = create_neb_path(images[0], images[-1], N_images, climb=climb, parallel=parallel)
                
            else:
                c.log(f"A new NEB calculation for {start_traj.stem} will be performed.")
            
                initial_row = db.get(iID)
                initial: Atoms = initial_row.toatoms()
                final_row = db.get(fID)
                final: Atoms = final_row.toatoms()

                # Order the endpoints correctly
                initial.append(initial.pop(vf-1))
                final.append(final.pop(vi))

                # Create images and master directory
                neb = create_neb_path(initial, final, N_images, climb=climb, parallel=parallel)

                write(str(start_traj), neb.images)
                c.log(f"Created the trajectory file {start_traj.stem} for the NEB calculation.")
       
    counter = len(trajs) + 1 
    traj_file = neb_dir / f"{neb_dir.name}_{counter}.{path}"
    
    c.log(f"Trajectory file {traj_file.stem} will be used for the NEB calculation.")
    c.log(f"Starting NEB calculation for {name}")
    c.log(f"Using climbing image: {climb}, FIRE optimizer: {fire}, fmax: {fmax}")

    # Single point calculations for endpoints
    initial = setup_run(neb.images[0], neb_dir / 'init', vasp)
    initial.get_potential_energy()

    final = setup_run(neb.images[-1], neb_dir / 'final', vasp)
    final.get_potential_energy()

    # Create calculator for internal images
    for i, img in enumerate(neb.images[1:-1]):
        img = setup_run(img, neb_dir / f"iimg{i+1:02}", vasp)
        img.get_potential_energy()

    # Run path relaxation with checker
    traj_log = neb_dir / f"{neb_dir.stem}.log"
    run_neb_with_checker(neb, traj_file, traj_log, fire, fmax)

    # Post-calculation processing
    
    converged_traj = read(f"{Path(traj_file)}@-{N_images+2}:")
    energies = np.asarray([image.get_potential_energy() for image in converged_traj])

    nebtool = NEBTools(converged_traj)
    nebtool.plot_bands(label=str(neb_dir / neb_dir.stem))

    # Post-calculation processing
    _, _, _, interpolated_energies, _ = fit_images(converged_traj)
            
    # Get the initial and final energies
    ini_energy = interpolated_energies[:1][0]
    fin_energy = interpolated_energies[-1:][0]
    delta_e = fin_energy - ini_energy
    ts_energy = np.amax(interpolated_energies)
    ts_index = np.argmax(interpolated_energies)
    min_energy = np.amin(interpolated_energies)
    min_index = np.argmin(interpolated_energies)
    c.log(f"Initial energy: {ini_energy:.2f} eV, Final energy: {fin_energy:.2f} eV, Delta E: {delta_e:.2f} eV")
    c.log(f"Min energy: {min_energy:.2f} eV at index {min_index}")
    c.log(f"Transition state energy: {ts_energy:.2f} eV at index {ts_index}")

    #Is an INTERmediate image the local MINimum?
    inter_min = 'yes'
    # Account for the cases where the minima or the transition states are very close to the endpoints
    rtol = 1e-4
    if np.isclose(ini_energy, min_energy, rtol=rtol) or np.isclose(fin_energy, min_energy, rtol=rtol):
        inter_min = 'no'
        c.log(f"The local minima are at an endpoint")
        Ef = ts_energy - ini_energy
        Er = ts_energy - fin_energy
        c.log(f'{Ef:.2f} = {ts_energy:.2f} - {ini_energy:.2f}')
        c.log(f'{Er:.2f} = {ts_energy:.2f} - {fin_energy:.2f}')
    
    elif np.isclose(ini_energy, ts_energy, rtol=rtol) or np.isclose(fin_energy, ts_energy, rtol=rtol):
        c.log(f"The transition state is at an endpoint")
        Ef = fin_energy - min_energy 
        Er = ini_energy - min_energy 
        c.log(f'{Ef:.2f} = {fin_energy:.2f} - {min_energy:.2f}')
        c.log(f'{Er:.2f} = {ini_energy:.2f} - {min_energy:.2f}')
        
    else:
        c.log(f"The minima are not at the endpoints, nor the transition state")
        from_left = np.amin(interpolated_energies[:ts_index])
        from_right = np.amin(interpolated_energies[ts_index+1:])
        Ef = ts_energy - from_left
        Er = ts_energy - from_right
        c.log(f'{Ef:.2f} = {ts_energy:.2f} - {from_left:.2f}')
        c.log(f'{Er:.2f} = {ts_energy:.2f} - {from_right:.2f}')
        
    c.log(f"Forward energy barrier: {Ef:.2f}, Reverse energy barrier: {Er:.2f}")
    barrier = max(Ef, Er)

    # If the calculated barrier is positive, compute KRA correction
    kra_energy = abs(barrier - 0.5 * abs(abs(Ef) - abs(Er)))
    c.log(f"The energy barrier is : {barrier:.3f}")
    c.log(f"KRA barrier: {kra_energy:.2f} eV")
    
    ts = converged_traj[2]     
    db_id = update_or_write(db, ts, name=neb_name, dir=neb_dir.as_posix(),
                            forward_e=Ef, delta_e=delta_e, reverse_e=Er, inter_min=inter_min,
                            barrier=barrier, kra_energy=kra_energy)

    # Copy the database back to the home directory
    home_db = here / "structures/hexag_perovs_wdiscards.db"
    shutil.copy(db_path, home_db)

    # Move the important files from scratch to home
    for file in neb_dir.glob("*"):
        if file.is_file() and file.name != 'WAVECAR':
            # Take the current path and move it to the home directory
            file_parts = list(file.parts)
            new_path = file_parts[:2] + ['energy'] + file_parts[3:]
            new_file = Path(*new_path)
            new_file.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(file, new_file)
    
    # Store the geometries of each image in the neb_geometries folder
    if climb:
        neb_geometries = here / "cineb_geometries" / name
    else:
        neb_geometries = here / "neb_geometries" / name
    
    neb_geometries.mkdir(exist_ok=True, parents=True)
    # Store the images as individual cif files
    write(f"{neb_geometries}/{name}_ini.cif", converged_traj[0])
    write(f"{neb_geometries}/{name}_fin.cif", converged_traj[-1])
    for i, image in enumerate(converged_traj[1:-1]):
        write(f"{neb_geometries}/{name}_img0{i+1}.cif", image)

    fig, ax = plt.subplots()
    ts_energy -=  min_energy
    fit_energies = np.array(interpolated_energies) - min_energy
    energies = np.array(energies - np.amin(energies))

    # Normalize the path
    x0 = np.linspace(0, 1, len(fit_energies))
    x1 = np.linspace(0, 1, len(energies))
    ax.plot(x1, energies, 'k-', linewidth=1, label='Calculated energies')
    ax.plot(x0, fit_energies, 'b--', linewidth=2, label='Fitted energies')
    # Plot the values of the energies for each of the images
    for i, energy in enumerate(energies):
        ax.plot(x1[i], energy, 'go')
    ax.plot(x0[ts_index], ts_energy, 'ro', label=f'Transition state')
    ax.set_xlabel('Reaction coordinate')
    ax.set_ylabel('Energy (eV)')
    ax.legend()
    plt.title(f"{name}, Barrier: {barrier:.2f} eV, Delta E: {delta_e:.2f} eV, KRA: {kra_energy:.2f} eV")
    plt.savefig(f"{neb_geometries}/{name}.png", bbox_inches='tight', dpi=300)
    plt.savefig(f"{neb_dir}/{name}.png", bbox_inches='tight', dpi=300)
    
    return True, {'trajectory': str(traj_file.relative_to(neb_dir)), 'neb_id': int(db_id)}

def monitor_log(logfile: Path, n_last=10, threshold=1e-3) -> bool:
    """
    Monitor the last `n_last` iterations in the log file to check for convergence.
    If the change in energy or forces is below `threshold`, return True for stagnation.
    """
    energies, forces = [], []

    with open(logfile, 'r') as f:
        lines = f.readlines()
        for line in lines[-n_last:]:
            if 'Step' in line:
                continue
            parts = line.split()
            # Some forces add a consistency check at the end
            energy = parts[-2][:-1] if not parts[-2][-1].isnumeric() else parts[-2] 
            energies.append(float(energy))
            forces.append(float(parts[-1]))

    if len(energies) < n_last:
        return False  # Not enough data yet

    energy_change = max(abs(energies[i] - energies[i+1]) for i in range(len(energies) - 1)) 
    force_change = max(abs(forces[i] - forces[i+1]) for i in range(len(forces) - 1))

    return energy_change < threshold and force_change < threshold

def run_neb_with_checker(neb: NEB, traj_file: Path, logfile: Path, fire: bool, fmax: float,
                         n_iter_check: int = 10, n_checks: int = 0, max_restarts: int = 3) -> None:
    """
    Run the NEB calculation with a convergence checker.
    The checker monitors the last `n_iter_check` iterations in the log file.
    It adjusts the parameters or switches the optimizer if stagnation is detected.
    """
    def setup_optimizer(fire, neb):
        dt = 0.05 #/ 2
        maxstep = 0.2 #/ 2
        dtmax = 0.3 #/ 2
        Nmin = 10
        finc = 1.05 #* 0.9
        fdec = 0.4 #* 1.1
        return FIRE(neb, trajectory=str(traj_file), logfile=logfile.as_posix(), dt=dt, maxstep=maxstep, dtmax=dtmax, 
                Nmin=Nmin, finc=finc, fdec=fdec, force_consistent=False) if fire else BFGS(neb, trajectory=str(traj_file), logfile=logfile.as_posix())

    optimizer = setup_optimizer(fire, neb)
    optimizer.run(fmax=fmax)

    check_counter = 0
    while check_counter < max_restarts:
        if monitor_log(logfile, n_last=n_iter_check):#, threshold=fmax / 10):
            c.log("Stagnation detected, adjusting optimizer parameters.")
            fire = not fire
            optimizer = setup_optimizer(fire, neb)
            optimizer.run(fmax=fmax)
            check_counter += 1
        else:
            break

    if check_counter >= max_restarts:
        c.log(f"Calculation stagnated after {max_restarts} restarts.")

def setup_run(atoms: Atoms, directory: Path, vasp_params: dict) -> Atoms:
    """
    Setup and run the DFT calculation for the given atoms object.

    Args:
        atoms (Atoms): The atomic structure to be calculated.
        directory (Path): The directory where the calculation will be performed.
        vasp_params (dict): Additional parameters for the VASP calculator.

    Returns:
        Atoms: The atoms object with the calculator attached.
    """
    directory.mkdir(parents=True, exist_ok=True)
    calc = create_Vasp_calc(atoms, 'PBEsol', directory)
    vasp_settings = {
            "ibrion" : -1,
            "nsw" : 1,
            "kpar" : nnodes,
            "ncore" : ncore,
    }
    vasp_settings.update(vasp_params)
    calc.set(**vasp_settings)
    set_magnetic_moments(atoms)
    return atoms

def update_or_write(db: SQLite3Database, atoms: Atoms, name: str, **kwargs):
    if db.count(name=name) > 0:
        ID = next(db.select(name=name)).id
        db.update(ID, atoms, name=name, **kwargs)
        return ID
    return db.write(atoms, name=name, **kwargs)

# if __name__ == '__main__':
#     from argparse import ArgumentParser
#     parser = ArgumentParser()
#     parser.add_argument('db_id', type=int, default=None)
#     parser.add_argument('vasp', type=dict, default={})
    
#     main(**vars(parser.parse_args()))