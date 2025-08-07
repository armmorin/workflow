from os import environ
from pathlib import Path
from typing import Optional, Tuple
from xml.etree.ElementTree import ParseError
from rich.console import Console
from ase.calculators.calculator import ReadError
from ase import Atoms
from ase.db import connect
from ase.db.sqlite import SQLite3Database
from ase.io import read
from typing import List
import shutil

from herculestools.dft import (
    RunConfiguration,
    create_Vasp_calc,
    set_magnetic_moments,
)

here = Path(__file__).resolve().parent.parent

c = Console()
nnodes = int(environ['SLURM_NNODES'])
ncore = int(environ['NCORE'])

# UPDATED for newer versions of HERCULEStools and PerQueue
def main(db_id: Optional[int] = None, initial_vac: Optional[int] = None,
         final_vac: Optional[int] = None, vasp: dict = {} ,**kwargs) -> Tuple[bool, Optional[dict]]:
    """
    Perform the preNEB calculation for a given structure.

    Args:
        db_id (Optional[int]): The database id of the structure.
        initial_vac (Optional[int]): The index of the initial vacancy.
        final_vac (Optional[int]): The index of the final vacancy.
        **kwargs: Additional keyword arguments.

    Returns:
        Tuple[bool, Optional[dict]]: A tuple containing a boolean indicating whether the calculation was successful,
        and an optional dictionary with the initial and final ids, initial and final vacancies.

    """
    
    db_path = RunConfiguration.structures_dir / "hexag_perovs_wdiscards.db"
    db: SQLite3Database

    db = connect(db_path)

    c.log(f"Starting preNEB calculation for {db_id}")
    # Get the name of the structure from the database.
    db_name = db.get(db_id).name
    names_list = db_name.split('_')
    name = '_'.join(names_list[:2])

    # Set up directories:
    i_direc = Path(RunConfiguration.home / f"preNEB/{name}_init")
    f_direc = Path(i_direc.parent / f"{name}_final")

    atoms = db.get(f"name={name}_r").toatoms()

    # This has to run when the initial and final vacancies are given.
    if initial_vac is not None and final_vac is not None:
        if initial_vac > final_vac:
            print("For the NEB we require that the index of the initial vacancy "
                "is smaller than the final. This is not the case, so they have "
                "been swapped!!")
            # Pythonic swapping
            initial_vac, final_vac = final_vac, initial_vac

    init = start_run(atoms=atoms, direc=i_direc, vacancy=initial_vac, vasp=vasp, **kwargs)
    i_energy = init.get_potential_energy()

    final = start_run(atoms=atoms, direc=f_direc, vacancy=final_vac, vasp=vasp, **kwargs)
    f_energy = final.get_potential_energy()

    # Move the important files from scratch to home
    direc_list = [i_direc, f_direc]
    for direc in direc_list:
        for file in direc.glob("*"):
            if file.is_file() and file.name != 'WAVECAR':
                # Take the current path and move it to the home directory
                file_parts = list(file.parts)
                new_path = file_parts[:2] + ['energy'] + file_parts[3:]
                new_file = Path(*new_path)
                new_file.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy(file, new_file)

    # Get the energy difference from the two positions.
    dE = f_energy - i_energy
    
    # Get the oxygen vacancy formation energy
    o_energy_vi, o_energy_vf = get_oxy_vac_energy(name, db, i_energy, f_energy)
    
    # Save the result to the database
    iID = update_or_write(db, Atoms(init),  name+"_vi", dir=str(i_direc), delta_e = dE, o_vac_energy = o_energy_vi)
    fID = update_or_write(db, Atoms(final), name+"_vf", dir=str(f_direc), delta_e = dE, o_vac_energy = o_energy_vf)

    # Copy the database back to the home directory
    home_db = here / "structures/hexag_perovs_wdiscards.db"
    shutil.copy(db_path, home_db)
    
    print(f"The energy difference between images is :{dE:.3f} eV")

    return True, {
        'initial_id': iID,
        'final_id': fID,
        'initial_vac': initial_vac,
        'final_vac': final_vac
        }

def read_poscar(poscar: Path) -> Atoms | None: 
    try:
        atoms: Atoms = read(poscar) 
    except Exception as e:
        raise POSCARError(f"Error in {poscar}") from e

    return atoms if atoms is not None else None

class POSCARError(Exception): 
    pass

## DEFINING A FUNCTION TO REGULATE THE PROCESS.
def start_run(atoms: Atoms | List[Atoms] | None, direc: Path, 
              vacancy: Optional[int] = None, vasp: dict={}) -> Atoms: 

    # Making a dictionary of the paths to the files.
    files = {
        'contcar' : direc / "CONTCAR",
        'poscar'  : direc / "POSCAR",
        'vasprun' : direc / "vasprun.xml"
    }

    try:
        atoms = read_poscar(files['poscar'])
        atoms = read(files['contcar']) 
        atoms = read(files['vasprun']) 
        print(f"A structure with {len(atoms)} atoms will start from {direc.name}")
        setup_run(atoms, direc, vasp)

    except (POSCARError, IndexError):
        atoms = start_from_scratch(atoms, direc, vacancy, vasp) # type: ignore

    except (ReadError, FileNotFoundError, ParseError):
        print(f"No CONTCAR or vasprun.xml found in {direc.name}")
        setup_run(atoms, direc, vasp)

    return atoms

# When starting from scratch make sure you are already dealing with a supercell.
def start_from_scratch(atoms: Atoms, direc: Path, vacancy: Optional[int] = None, vasp: dict = {}) -> Atoms:
    # Here. Add a check to see if the system already has vacancies from the atoms object.
    if len(atoms) == 32:
        print(f"Supercell is generated in {direc.name}")
        sc = atoms.repeat((2, 2, 1))
        atoms = sc.copy()
        atoms.pop(vacancy) # type: ignore
    else:
        print(f"Restarting from supercell structure in {direc.name}")

    setup_run(atoms, direc, vasp)

    return atoms

def setup_run(atoms: Atoms | List[Atoms] | None, direc: Path, vasp: dict = {}) -> Atoms:
    set_magnetic_moments(atoms)
    calc = create_Vasp_calc(atoms, 'PBEsol', direc)
    calc.set(
            ediffg=-0.05,
            ibrion=1,
            nsw=500,
            kpar=nnodes,
            ncore=ncore,
            **vasp)

    return atoms

def update_or_write(db: SQLite3Database, atoms: Atoms, name: str, **kwargs):
    if db.count(name=name) > 0:
        ID = next(db.select(name=name)).id
        db.update(ID, atoms, name=name, **kwargs)
        return ID

    return db.write(atoms, name=name, **kwargs)

def get_oxy_vac_energy(entry: str, db: SQLite3Database, ini_energy: float, fin_energy: float) -> Tuple[float, float]:
    # Get the oxygen vacancy formation energy
    relax_row = db.get(name=f"{entry}_r")
    energy_r = relax_row.energy

    # The formation energy of the oxygen vacancy at the different sites
    o_energy_vi = ini_energy - (4 * energy_r) + (-4.52) # This would be with the energy of a single oxygen atom
    o_energy_vf = fin_energy - (4 * energy_r) + (-4.52) 
    
    return o_energy_vi, o_energy_vf

if __name__ == '__main__':
    from sys import argv
    index = int(argv[1])
    args = [int(arg) for arg in argv[2:]]

    print(main(index, *args))
