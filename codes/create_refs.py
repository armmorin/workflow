from mp_api.client import MPRester
from pymatgen.io.ase import AseAtomsAdaptor
from ase.db import connect
from ase.atoms import Atoms
from ase.db.sqlite import SQLite3Database
from perqueue.queue import PersistentQueue
from perqueue.task_classes.task import Task
from herculestools.dft import RunConfiguration
from time import sleep

# Set up the MPRester object.
API_KEY = "X7D2MtOD2Z7lrrcPlKJv6XIVOXRtCiYH"

# Set up the database connection.
from_db = RunConfiguration.structures_dir / "hexag_perovs_wdiscards.db"
db : SQLite3Database = connect(from_db)

# Create the list of elements to be queried. 
symbols = []
for entry in db.select(selection='barrier'):
    for symbol in entry.symbols:
        symbols.append(symbol)

# Remove duplicates and sort the list.
symbols = list(set(symbols))
symbols.sort()

# Connecting to the new references database.
ref_db = RunConfiguration.structures_dir / "references.db"
to_db: SQLite3Database = connect(ref_db)

def update_or_write(db: SQLite3Database, atoms: Atoms, name: str, **kwargs):
    if db.count(name=name) > 0:
        ID = next(db.select(name=name)).id
        db.update(ID, atoms, name=name, **kwargs)
        return ID

    return db.write(atoms, name=name, **kwargs)


with MPRester(api_key=API_KEY) as mpr:
    for symbol in symbols:
        # Query the MP database for the structure.
        computed_entries = mpr.get_entries(chemsys_formula_mpids=symbol, conventional_unit_cell=False,
                            additional_criteria={'is_stable': True, "energy_above_hull": (0.0, 0.0)})
        # Order the entries by number of sites.
        computed_entries = sorted(computed_entries, key=lambda x: x.structure.num_sites)
        print(symbol,len(computed_entries))
        try:            
            entry = computed_entries[0]
            struct = entry.structure
            # Get the structure and convert it to an ASE Atoms object.
            atoms = AseAtomsAdaptor.get_atoms(struct, msonable=False)
            print(atoms)
            ID = update_or_write(to_db, atoms, name=symbol)
                        
        except IndexError:
            # This is to handle the case where there are no entries for the symbol.
            # The entry has been added manually to the database. retrieve it and run the task.
            try:
                entry = to_db.get(selection=f"name={symbol}")
                atoms = entry.toatoms()
                ID = entry.id
            except KeyError:
                print(f"Skipping {symbol} due to no entries")
                continue
                
        task = Task(
            code="codes/calc_refs.py",
            args={'element': symbol, 'id': ID},
            resources="24:1:xeon24el8_test:30m",
            name=f"Ref_{symbol}")

        with PersistentQueue() as pq:
            pq.submit(task)        
            sleep(0.5)
        
