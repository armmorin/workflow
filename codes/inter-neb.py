# Prepare an interpolation between two images, and perform a singe point calculation of the created images.

from pathlib import Path
from os import environ
from rich.console import Console
from ase.io import read, write
import shutil
from ase.db import connect
from herculestools.dft import RunConfiguration as RC, create_Vasp_calc, set_magnetic_moments
from ase.neb import NEB

structures_dir = RC.structures_dir
db_path = structures_dir / "hexag_perovs_wdiscards.db"

here = Path(__file__).resolve().parent.parent

c = Console()
nnodes = int(environ['SLURM_NNODES'])
ncore = int(environ['NCORE'])

def main(image, direc, vasp={}, **kwargs):
    
    # Create the VASP calculation
    calc = create_Vasp_calc(structure=image, xc='PBEsol', directory=direc)
    calc.set(
        ibrion=-1,
        nsw=1,
        kpar=nnodes,
        ncore=ncore,
        **vasp)
    set_magnetic_moments(image)
    energy = image.get_potential_energy()
    return image, energy 

# Extract the system from the database
with connect(db_path) as db:
    row = db.get(id=3461) # La7Ge4ScO20_p1
    row_dir = Path(row.dir)
    name = row_dir.name
    # Get the initial and final images for the interpolation
    initial = read(row_dir / "iimg03/vasprun.xml")
    final = read(row_dir / "final/vasprun.xml")
    
    # Prepare the interpolated imges for a single point calculation
    images = [initial]
    images += [initial.copy() for i in range(5)]
    images += [final]
    
        # Interpolate between the two images
    neb = NEB(images)
    neb.interpolate(mic=True)
    
    sp = Path(f'single_point')
    sp.mkdir(parents=True, exist_ok=True)
    write(sp/f'{name}_interpolated.traj', images)
    
    # Perform the single point calculation
    energies = []
    for i, image in enumerate(images):
        direc = sp / f'{name}_{i}'
        direc.mkdir(parents=True, exist_ok=True)
        image, energy = main(image, direc)
        energies.append(energy)
        
        c.log(f"Image {i} has energy {energy}")
        # Move the important files from scratch to home
        for file in direc.glob("*"):
            if file.is_file() and file.name != 'WAVECAR':
                # Take the current path and move it to the home directory
                file_parts = list(file.parts)
                new_path = file_parts[:2] + ['energy'] + file_parts[3:]
                new_file = Path(*new_path)
                new_file.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy(file, new_file)
        
            
        
        