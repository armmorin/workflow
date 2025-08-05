from sys import argv
from os import environ
from rich.console import Console
from pathlib import Path
from ase import Atoms
import shutil
from ase.db import connect
from ase.db.sqlite import SQLite3Database
from pymatgen.io.vasp import Vasprun
import numpy as np
from pymatgen.electronic_structure.core import OrbitalType, Spin
from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.electronic_structure.dos import CompleteDos
import matplotlib.pyplot as plt
from ase.data import chemical_symbols
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

# Define a function to calculate the p-band center of an element in the DOS
def calculate_p_band_center(complete_dos: CompleteDos, element="O"):
    pdos = complete_dos.get_element_spd_dos(element)
    p_dos = pdos[OrbitalType.p]
    energies = p_dos.energies
    densities = p_dos.densities[Spin.up]  # Assuming non-spin-polarized calculation

    p_band_width = np.sum(energies * densities) / np.sum(densities)
    
    return p_band_width - complete_dos.efermi  # Reference to Fermi level

# For good measure, generate a plot of the PDOS and visualize the distribution of states.
# Plot the projected density of states from a given id in the database
def plot_dos(atoms: Atoms, name: str, complete_dos: CompleteDos):
    
    dos_dict = {}
    set_of_atoms = set(sorted(atoms.get_atomic_numbers())) # Find another way to sort the atoms, by numbers first
    set_of_atoms = [chemical_symbols[i] for i in set_of_atoms]
    cdos = complete_dos
    plotter = DosPlotter()
    gap = cdos.get_gap()
    for atom in set_of_atoms:
        element_dict = cdos.get_element_spd_dos(atom)
        plotter.add_dos_dict(element_dict)
        dos_dict[atom] = plotter.get_dos_dict()    
    
    # Plot the DOS with matplotlib
    plot_name = f"pdos_{name}.png"
    fig, ax = plt.subplots()
    cmap = plt.get_cmap('rainbow_r')
    # Get the total combinations of elements and orbitals from the dictionary
    tuples = [(i, j) for i in dos_dict for j in dos_dict[i]]
    permutations = len(tuples)
    #print(permutations)
    #ax.set_prop_cycle(color=[cmap(i) for i in np.linspace(0, 1, 8)])
    colors = cmap(np.linspace(0, 1, permutations))
    # Plot the DOS from the dictionary by element and its orbitals
    for i,element in enumerate(dos_dict):
        for j,orbital in enumerate(dos_dict[element]):
            # Give i,j a unique value from 0 to permutations
            k = i * len(dos_dict[element]) + j

            energies = np.asarray(dos_dict[element][orbital]['energies']) #+ fermi_energy # The energies are already shifted to the Fermi level
            densities_up = np.asarray(dos_dict[element][orbital]['densities']['1'])
            # If the values are all zero for the densities, skip the plot
            if np.allclose(densities_up, 0):
                continue
            densities_dn = np.asarray(dos_dict[element][orbital]['densities']['-1']) * -1
            ax.plot(energies, densities_up, c=colors[k], linewidth=1.25, label=f"{element} {orbital}")
            ax.plot(energies, densities_dn, c=colors[k], linestyle='--', linewidth=0.75)

    # Set the legend outside the plot
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_xlim(-15, 10)
    ax.set_ylim(-80, 80)
    ax.annotate(f'Band gap = {gap:.2f} eV', xy=(0.4, 0.05), xycoords='axes fraction')
    ax.grid(color='gray', linestyle='--', linewidth=0.25, which='both')
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel('Density of states')
    ax.set_title(f'Projected DOS for pristine {name}')
    #fig.show()
    fig.savefig(f'neb_geometries/{name}/{plot_name}', bbox_inches='tight', dpi=600)
    return c.log(plot_name)

# --- SCRIPT BODY ---

# Get the specified structure from the database
def main(db_id: int, vasp: dict= {}, **kwargs):
    db: SQLite3Database = connect(db_path)
    with db:
        row = db.get(db_id)
        direc = Path(row.dir)
        name = direc.name
        structure = row.toatoms()

    extent = [-15, 15]
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

    # POST-PROCESSING - Calculate the p-band center for the oxygen atoms and plot the pDOS

    if (vaspxml := dos_dir / 'vasprun.xml').is_file():
        vasprun = Vasprun(vaspxml.as_posix(), parse_dos=True, parse_eigen=False, parse_projected_eigen=False, parse_potcar_file=False)
        c.log(f'Processing {name} with Fermi level at {vasprun.efermi:.2f} eV')
        cdos = vasprun.complete_dos
        
        p_band_center = calculate_p_band_center(cdos, element='O')
        c.log(f"Oxygen p-band center: {p_band_center:.2f} eV")#relative to Fermi level")
        
        plot_dos(structure, name, cdos)
    else:
        c.log(f"VASP run failed for {name}")
        return False, {'db_id': db_id}
    
    # Save the result to the database
    with db:
        db.update(db_id, dos='done', p_band_center=p_band_center, e_fermi=vasprun.efermi)

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
    