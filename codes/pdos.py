from sys import argv
#from threading import Timer
from os import environ
from rich.console import Console
from pathlib import Path
from ase import Atoms
import shutil
from ase.db import connect
from ase.db.row import AtomsRow
from ase.db.sqlite import SQLite3Database
from herculestools.dft import (
    RunConfiguration as RC,
    create_Vasp_calc,
    set_magnetic_moments,
)

c = Console()
nnodes = int(environ['SLURM_NNODES'])
ncore = int(environ['NCORE'])

#here = Path.cwd() # Where the code is located
here = Path(__file__).resolve().parent.parent
home = RC.home # Where the calculations are stored
structures = RC.structures_dir # Where the structures are stored

# Connect to the database
db_path = structures / 'hexag_perovs_wdiscards.db'

# --- SCRIPT BODY ---

# Get the specified structure from the database
def main(db_id: int, vasp: dict= {}, **kwargs):
    db: SQLite3Database = connect(db_path)
    with db:
        row = db.get(db_id)
        name = row.name
        direc = Path(row.dir)
        structure = row.toatoms()
        #hsh = row.hash
    
    extent = [-15, 20]
    dos_spacing = 0.01


    # Create the Vasp calculator
    dos_dir = home / direc / "dos"
    dos_dir.mkdir(parents=True, exist_ok=True)
    c.log(f"Calculating DOS for {name} in {dos_dir}")

    [f.unlink() for f in dos_dir.iterdir() if not f.name == 'WAVECAR']

    # Extract and calculate parameters
    emin, emax, *_ = sorted(extent)
    nedos = round((emax-emin) / dos_spacing) + 1

    # NSCF run
    calc = create_Vasp_calc(structure, 'PBEsol', dos_dir, density=40)
    vasp_settings = {
        'istart':1,
        'nsw':0,
        'lorbit':10,
        'nedos':nedos,
        'emin':emin,
        'emax':emax,
        'ismear':-5
    }
    vasp_settings.update(vasp)
    calc.set(**vasp_settings)
    set_magnetic_moments(structure)

    structure.get_potential_energy()

    # Delete the WAVECAR
    (dos_dir / "WAVECAR").unlink(True)

    c.log(f"Calculated DOS for {name} in directory {dos_dir}")
    
    # Move the important files from scratch to home
    for file in dos_dir.glob("*"):
        if file.is_file() and file.name != 'WAVECAR':
            # Take the current path and move it to the home directory
            file_parts = list(file.parts)
            new_path = file_parts[:2] + ['energy'] + file_parts[3:]
            new_file = Path(*new_path)
            new_file.parent.mkdir(parents=True, exist_ok=True)
            try:
                shutil.copy(file, new_file)
            except shutil.SameFileError:
                pass
        
     # Save the result to the database
    with connect(db_path) as db:
        db.update(db_id, dos='done')
        
    # Copy the database back to the home directory
    home_db = here / "structures/hexag_perovs_wdiscards.db"
    try:
        shutil.copy(db_path, home_db)
    except shutil.SameFileError:
        pass
    return True, {'db_id': db_id}

if __name__ == '__main__':
    from sys import argv
    #args = [float(arg) for arg in argv[1:]]
    #print(main(*args))

    print(main(int(argv[1])))
    