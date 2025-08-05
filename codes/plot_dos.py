from ase.io import read
from ase.db import connect
import matplotlib.pyplot as plt
import numpy as np
from pymatgen.electronic_structure.plotter import DosPlotter
#from pymatgen.electronic_structure.dos import CompleteDos
from pymatgen.io.vasp import Vasprun
from ase.db.sqlite import SQLite3Database
from pathlib import Path
import matplotlib
from ase.data import chemical_symbols
matplotlib.use("Cairo")

# Connect to the database
db_path = "structures/hafnia_alloys.db"
db = connect(db_path)

# Plot the projected density of states from a given id in the database
def main(db_id: int, db: SQLite3Database, **kwargs):
    
    # Get thte structure from the database
    row = db.get(db_id)
    name = row.name
    if row.dos != 'done':
        print(f"DOS calculation for {name} not done")
        return None
    row_dir = Path(row.dir)
    
    # Read the DOSCAR file
    dos_path = row_dir / 'dos'
    xml_path = dos_path / 'vasprun.xml'
    if not xml_path.exists():
        print(f"vasprun.xml not found in {xml_path}")
        return None
    dos_dict = {}
    atoms = read(xml_path, index=0)
    set_of_atoms = set(sorted(atoms.get_atomic_numbers())) # Find another way to sort the atoms, by numbers first
    set_of_atoms = [chemical_symbols[i] for i in set_of_atoms]
    vasprun = Vasprun(dos_path / 'vasprun.xml')
    cdos = vasprun.complete_dos
    plotter = DosPlotter()
    
    for atom in set_of_atoms:
        element_dict = cdos.get_element_spd_dos(atom)
        plotter.add_dos_dict(element_dict)
        dos_dict[atom] = plotter.get_dos_dict()    
    
    # Store the contents of the dos_dict in a file
    with open(dos_path / 'pdos.json', 'w') as f:
        f.write(str(dos_dict))
    
    # Plot the DOS with matplotlib
    plot_name = f"pdos_{name}.png"
    fig, ax = plt.subplots()
    cmap = plt.get_cmap('rainbow')
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
            #print(i, j, k)
            #fermi_energy = dos_dict[element][orbital]['efermi']
            energies = np.asarray(dos_dict[element][orbital]['energies']) #+ fermi_energy # The energies are already shifted to the Fermi level
            densities_up = np.asarray(dos_dict[element][orbital]['densities']['1'])
            # If the values are all zero for the densities, skip the plot
            if np.allclose(densities_up, 0):
                continue
            densities_dn = np.asarray(dos_dict[element][orbital]['densities']['-1']) * -1
            ax.plot(energies, densities_up, c=colors[k], label=f"{element} {orbital}")
            ax.plot(energies, densities_dn, c=colors[k], linestyle='--')
#            break
#        break
    # Set the legend outside the plot
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_xlim(-10, 10)
    ax.set_ylim(-75, 75)
    ax.annotate(f'Band gap = {row.gap:.2f} eV', xy=(0, -70), xytext=(-40, 0), textcoords='offset points')
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel('Density of states')
    ax.set_title(f'Projected DOS for {name}')
    fig.savefig(plot_name, bbox_inches='tight', dpi=300)
    
    return print(plot_name)

if __name__ == "__main__":
    from sys import argv
    args = [int(i) for i in argv[1:]]
    if len(args) == 1:
        print(main(*args, db))
    else:
        for arg in args:
            print(main(arg, db))
        
    