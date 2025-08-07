from pathlib import Path
from typing import Tuple
from ase.db import connect
from ase.io import read
from emmet.core.thermo import ThermoType
from mp_api.client import MPRester
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from sqlite3 import OperationalError
from herculestools.dft import RunConfiguration
from time import sleep

db_path = RunConfiguration.structures_dir / "hexag_perovs_wdiscards.db"
def main(**kwargs) -> Tuple[bool, None]:
    
    db_id = kwargs['db_id']
    with connect(db_path) as db:
        direc = Path(db.get(db_id).dir)
    
    atoms = read(direc / 'CONTCAR')
    entry = Vasprun(direc / 'vasprun.xml').get_computed_entry()
    compat = MaterialsProject2020Compatibility()
    entry.parameters['run_type'] = 'GGA'
    entry = compat.process_entry(entry) or entry  # function can return None    
#    end_point = MPRester("0gn6w1hTmwhMsHFGj") # <--- Legacy API key
    end_point = MPRester("FP32VD2CnVamajyKu5HU28pMBS7f3p1p")  # <--- New API key

    chemsys = list({*atoms.symbols})
    entries = end_point.get_entries_in_chemsys(chemsys, additional_criteria = {"thermo_types": [ThermoType.GGA_GGA_U]})
    entries.append(entry)

    pd = PhaseDiagram(entries)
    e_above_hull = pd.get_e_above_hull(entry)
    
    for attempt in range(5):
        with connect(db_path) as db:
            try:
                db.update(db_id, e_above_hull=e_above_hull)
                break
            except OperationalError as e:
                print(f"Attempt {attempt + 1} failed. Reason:{e}")
                sleep(1)

    print(f"Energy above the hull is: {e_above_hull:.3f}")
    #if e_above_hull < hull_tol:
    #    return True, {'db_id': db_id}
    return True, {'db_id': db_id}

from sys import argv
if __name__ == "__main__":
    main(int(argv[1]))
