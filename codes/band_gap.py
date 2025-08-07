from pathlib import Path
from typing import Tuple

from ase.db import connect
from ase.dft.bandgap import bandgap
from herculestools.dft import RunConfiguration
from herculestools.dft import load_run

def main(db_id:int = 1, gap_tol: float=1.0) -> Tuple[bool, None]:
    
    db_path = RunConfiguration.structures_dir / "hexag_perovs_wdiscards.db"
    with connect(db_path) as db:
        direc = Path(db.get(db_id).dir)

    atoms = load_run(direc)
    bg = bandgap(atoms.calc)

    with connect(db_path) as db:
        db.update(db_id, bandgap=str(bg), gap=bg[0])

    print(f"The band gap is: {bg}")
    if bg[0] > gap_tol:
        return True, None
    return True, None

if __name__ == "__main__":
    from sys import argv
    main(int(argv[1]))
