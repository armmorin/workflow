from perqueue import PersistentQueue
from perqueue.selection import Selection
from ase.db import connect
from pathlib import Path

with PersistentQueue() as pq:
    entries = pq.get_entries()

s = Selection(states="s")
targets = s.filter(entries)
db = connect('structures/hexag_perovs_wdiscards.db')

print(f"Number of entries in db: {len(targets)}")
def rm_wavecar(id: int):
    entry = db.get(selection=f'id={id}')
    name = entry.name
    jobdir = Path(entry.dir)
    direc_list = [direc.parent for direc in jobdir.rglob("vasprun.xml")]
    for direc in direc_list:
        wavecar = direc / "WAVECAR"
        # chg = direc / "CHG"
        # procar = direc / "PROCAR"
        # chgcar = direc / "CHGCAR"
        files = [wavecar]
        for file in files:
            if file.exists():
                file.unlink()
                print(f"Removing {file.name} for {name}")

en_data = [en.data for en in targets]
keyslist = ['initial_id', 'final_id', 'neb_id']

# Get the values contained in the keyslist
items = [item[key] for item in en_data for key in keyslist if key in item]
# Get the unique values
items_list = list(set(items))
# Remove None from the list
items_list = [item for item in items_list if item is not None]

for item in items_list:
    try:
        rm_wavecar(item)
    except:
        print(f"Could not remove files for {item}")
        continue
