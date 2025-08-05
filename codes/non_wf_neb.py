#from argparse import ArgumentParser
from ase.db.sqlite import SQLite3Database
from os import environ
from ase.atoms import Atoms
from pathlib import Path
from rich.console import Console
from ase.io import read, write
from ase.optimize import FIRE, BFGS
from ase.neb import NEB
from ase.neb import NEBTools
from ase.calculators.vasp import Vasp
from ase.db import connect
from typing import Optional
import shutil
from ase.io.trajectory import Trajectory
import numpy as np
import matplotlib.pyplot as plt
from forcecurve import fit_images
#nnodes = int(environ['SLURM_NNODES'])

from herculestools.dft import (
    RunConfiguration,
    create_Vasp_calc,
    set_magnetic_moments,
    create_neb_path
)

here = Path(__file__).parent
rc_home = RunConfiguration.home
db_path = RunConfiguration.structures_dir / 'hexag_perovs_wdiscards.db'
c = Console()
# For the combined NEB, freeze the first and second layers.
def main(db_id: Optional[int], N_images: int=3, climb: bool=False, parallel: bool= False, fire: bool=True, fmax: float=0.05, vasp: dict = None):

    with connect(db_path) as db:
        row = db.get(id=db_id)
        if row is None:
            c.log(f"ID {db_id} not found in the database.")
            return
        neb_name = row.name[:-4]
    
    # Set the shared calculator flag
    shared_calc = not parallel

    
    # Check if the NEB directory exists, and if it does, remove empty trajectory files.
    neb_dir = rc_home / 'NEB' / neb_name
    path = "traj"
    for img_dir in neb_dir.glob(f"*{path}"):
        if img_dir.stat().st_size == 0:
            img_dir.unlink(missing_ok=True)

    trajs = [img for img in neb_dir.glob(f"*{path}") if img.is_file() and img.stat().st_size > 0]
    start_traj = neb_dir / f"{neb_name}_start.{path}"
    counter = len(trajs)
    vf = 31
    vi = 30
    if counter == 0:
        ini_row = db.get(name=f"{neb_name}_vi")
        initial = ini_row.toatoms()
        final_row = db.get(name=f"{neb_name}_vf")
        final = final_row.toatoms()

        # Order the endpoints correctly
        initial.append(initial.pop(vf-1))
        final.append(final.pop(vi))

        # Create images and master directory
        neb = create_neb_path(initial, final, N_images, climb=climb, parallel=parallel)

        write(str(start_traj), neb.images)
        
    else:
        traj = max(trajs, key= lambda a: a.stat().st_mtime)
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
        
    counter = len(trajs) + 1 
    traj_file = neb_dir / f"{neb_name}_{counter}.{path}"
    
    c.log(f"Trajectory file {traj_file.stem} will be used for the NEB calculation.")
    c.log(f"Starting NEB calculation for {neb_name}")
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
    logfile = neb_dir / f"{neb_name}.log"
        
    optimizer = setup_optimizer(traj_file, logfile, fire, neb)
    optimizer.run(fmax=fmax)
  
    # Post-calculation processing
    
    converged_traj = read(f"{Path(traj_file)}@-{N_images+2}:")
    energies = np.asarray([image.get_potential_energy() for image in converged_traj])

    nebtool = NEBTools(converged_traj)
    nebtool.plot_bands(label=str(neb_dir / neb_name))

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
    c.log(f"KRA barrier: {kra_energy:.2f} eV")
    
    ts = converged_traj[2]     
    db_id = update_or_write(db, ts, name=f"{neb_name}_neb", dir=neb_dir.as_posix(),
                            forward_e=Ef, delta_e=delta_e, reverse_e=Er, inter_min=inter_min,
                            barrier=barrier, kra_energy=kra_energy)

    # Copy the database back to the home directory
    home_db = here / "structures/hexag_perovs_wdiscards.db"
    shutil.copy(db_path, home_db)
    
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
    plt.title(f"{neb_name}, Barrier: {barrier:.2f} eV, Delta E: {delta_e:.2f} eV, KRA: {kra_energy:.2f} eV")
    plt.savefig(f"{neb_dir}/{neb_name}.png", bbox_inches='tight', dpi=300)

    # Move the important files from scratch to home
    for file in neb_dir.glob("*"):
        if file.is_file() and file.name != 'WAVECAR':
            # Take the current path and move it to the home directory
            file_parts = list(file.parts)
            new_path = file_parts[:2] + ['energy'] + file_parts[3:]
            new_file = Path(*new_path)
            new_file.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(file, new_file)
    
    c.log(f"The energy barrier is : {barrier:.3f}")
    
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
            #"kpar" : nnodes,
            #"ncore" : ncore,
    }
    vasp_settings.update(vasp_params)
    calc.set(**vasp_settings)
    set_magnetic_moments(atoms)
    return atoms

def setup_optimizer(traj_file, logfile, fire, neb):
    dt = 0.05 #/ 2
    maxstep = 0.2 #/ 2
    dtmax = 0.3 #/ 2
    Nmin = 10
    finc = 1.05 #* 0.9
    fdec = 0.4 #* 1.1
    return FIRE(neb, trajectory=str(traj_file), logfile=logfile.as_posix(), dt=dt, maxstep=maxstep, dtmax=dtmax, 
            Nmin=Nmin, finc=finc, fdec=fdec, force_consistent=False) if fire else BFGS(neb, trajectory=str(traj_file), logfile=logfile.as_posix())

def update_or_write(db: SQLite3Database, atoms: Atoms, name: str, **kwargs):
    if db.count(name=name) > 0:
        ID = next(db.select(name=name)).id
        db.update(ID, atoms, name=name, **kwargs)
        return ID
    return db.write(atoms, name=name, **kwargs)

# Call the main function in accordance with the command line arguments
if __name__ == "__main__":
    import sys
    from ast import literal_eval
    db_id = int(sys.argv[1])
    vasp = literal_eval(sys.argv[2])
    main(db_id, vasp=vasp)
    # parser = ArgumentParser(description="Perform NEB calculations for a given ID.")
    # parser.add_argument("args", nargs="*")
    # parser.add_argument("--db_id", type=int, help="The ID of the structure to perform the NEB calculation on.")
    # parser.add_argument("--N_images", type=int, help="The number of images in the NEB calculation.")
    # parser.add_argument("--climb", action="store_true", help="Use the climbing image method.")
    # parser.add_argument("--parallel", action="store_true", help="Use parallel calculations.")
    # parser.add_argument("--fire", action="store_true", help="Use the FIRE optimizer.")
    # parser.add_argument("--fmax", type=float, help="The maximum force convergence criterion.")
    # parser.add_argument("--vasp", type=dict, help="Additional VASP parameters.")
    # args = parser.parse_args()
    # main(*args.args, db_id=args.db_id, N_images=args.N_images, climb=args.climb, parallel=args.parallel, fire=args.fire, fmax=args.fmax, vasp=args.vasp)