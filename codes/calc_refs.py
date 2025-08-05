from ase.io import read, write
from ase.build import bulk
from ase.db import connect
#from ase.calculators.vasp import Vasp
from ase.db.sqlite import SQLite3Database
from ase import Atoms
from herculestools.dft import (
    RunConfiguration,
    create_Vasp_calc,
    set_magnetic_moments
)

# Set up the database connection.
db_path = RunConfiguration.structures_dir / "references.db"
db : SQLite3Database = connect(db_path)

# Set the home directory.
home = RunConfiguration.home
refs_dir = home / "references"
refs_dir.mkdir(exist_ok=True, parents=True)

# Define the main function. It reads from the database and writes the new reference structures.
def main(vasp: dict={}, **kwargs):
    # Get the element from the arguments.
    element = kwargs['element']
    atoms = db.get(id=kwargs['id']).toatoms()
    direc = refs_dir / element
    direc.mkdir(exist_ok=True, parents=True)
    atoms.write(direc / f"{element}.cif")
    set_magnetic_moments(atoms)
    calc = create_Vasp_calc(atoms, 'pbesol', direc, kwargs['element'])
    calc.set(**vasp,
            nelm=250,
            isif=3,
            ediffg=-0.05,
            nsw=250,
            ibrion=2
            )
    
    atoms.calc = calc
    atoms.get_potential_energy()

    db_id = update_or_write(db, atoms, name=kwargs['element'], directory=str(direc))  
    
    return True, {'id': int(db_id)}

def update_or_write(db: SQLite3Database, atoms: Atoms, name: str, **kwargs):
    if db.count(name=name) > 0:
        ID = next(db.select(name=name)).id
        db.update(ID, atoms, name=name, **kwargs)
        return ID

    return db.write(atoms, name=name, **kwargs)

if __name__ == "__main__":
    from sys import argv
    args = argv[1:]
    print(main(element=args[0], symmetry=args[1]))