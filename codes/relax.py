from os import environ
from pathlib import Path
from typing import Optional, Tuple
from rich.console import Console
from ase import Atoms
from ase.constraints import FixAtoms
from ase.db import connect
from ase.db.sqlite import SQLite3Database
from ase.spacegroup import get_spacegroup
from ase.io import read
import shutil

from herculestools.dft import (
    RunConfiguration,
    create_VASP_calc,
    set_magnetic_moments,
)

here = Path(__file__).resolve().parent.parent

c = Console()
nnodes = int(environ['SLURM_NNODES'])
ncore = int(environ['NCORE'])

def main(db_id: Optional[int] = None, rattle_std: float=0.05,
         sym_precision: float=0.1, distance: float=0.5, vasp: dict={}, **kwargs
         ) -> Tuple[bool, Optional[dict]]:
    
    # Load the structure from the database
    db_path = RunConfiguration.structures_dir / "hexag_perovs_wdiscards.db"
    db: SQLite3Database
    # Load the structure and extract fields
    print(f"db_id: {db_id}")
    with connect(db_path) as db:
        entry = db.get(db_id)
        atoms: Atoms = entry.toatoms()
        name = entry.name

    # Set up calculator and auxilliary things
    lets_run = True
    direc: Path = RunConfiguration.home / f"relaxations/{name}"
    if ((outcar := direc/"vasp.out")).exists():
        outcar_txt = outcar.read_text()
        if "reached required accuracy - stopping structural energy minimisation" in outcar_txt:
            lets_run = False
    
    if ((fl := direc/'vasprun.xml')).exists():
        try:
            atoms = read(fl)
        except:
            pass
    calc = create_VASP_calc(atoms, 'PBEsol', direc.name, direc)
    set_magnetic_moments(atoms)
    calc.set(
        nelm = 250,
        ibrion = 2,
        isif = 3,
        ediffg = -0.05,
        nsw = 250,
        **vasp)

    # Apply calculator and rattle atoms
    atoms.calc = calc

    # Rattle all atoms but the specific oxygen
    cn = FixAtoms([31])
    atoms.set_constraint(cn)
    atoms.rattle(rattle_std)  # <--- Tolerance #1
    atoms.set_constraint()
    # \/ Parameter #3
    atoms[31].position += distance * atoms.cell[0] / atoms.cell.lengths()[0]

    # Does the relaxation as well
    if lets_run:
        atoms.get_potential_energy()
    atoms.wrap()

    # Check that the structure is stable in the spacegroup
    spg = get_spacegroup(atoms, sym_precision)  # <--- Tolerance #2
    
    # Move the important files from scratch to home
    for file in direc.glob("*"):
        if file.is_file() and file.name != 'WAVECAR':
            # Take the current path and move it to the home directory
            file_parts = list(file.parts)
            new_path = file_parts[:2] + ['energy'] + file_parts[3:]
            new_file = Path(*new_path)
            new_file.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(file, new_file)
    
    # Save the result to the database
    with connect(db_path) as db:
        db_id = update_or_write(db, atoms, name+'_r', spacegroup=spg.symbol, dir=str(direc))

    # Copy the database back to the home directory
    home_db = here / "structures/hexag_perovs_wdiscards.db"
    shutil.copy(db_path, home_db)
    
    return True, {'db_id': db_id}

def update_or_write(db: SQLite3Database, atoms: Atoms, name: str, **kwargs):
    if db.count(name=name) > 0:
        ID = next(db.select(name=name)).id
        db.update(ID, atoms, name=name, **kwargs)
        return ID
    return db.write(atoms, name=name, **kwargs)

if __name__ == '__main__':
    from sys import argv
    index = int(argv[1])
    args = [float(arg) for arg in argv[2:]]

    print(main(index, *args))
