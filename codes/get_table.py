from ase.db import connect
import pandas as pd
from pathlib import Path
from collections import Counter
from ase import data
from ase.db.sqlite import SQLite3Database
import mendeleev as md
from decimal import Decimal
from ase import data
from collections import Counter

cwd = Path(__file__).parent.parent
db_path = cwd / 'structures/hexag_perovs_wdiscards.db'
#print(db_path)
db = connect(db_path.as_posix())
df = pd.DataFrame()

def get_oxy_vac_energy(entry: str, db: SQLite3Database):
    # Get the oxygen vacancy formation energy
    relax_row = db.get(name=f"{entry}_r")
    energy_r = relax_row.energy
    # The energy of the structure with the oxygen vacancies
    energy_vi = db.get(name=f"{entry}_vi").energy
    energy_vf = db.get(name=f"{entry}_vf").energy
    
    # The energies of the hydrogen and the water molecules
    mol_db = connect('structures/molecules.db')
    h2_energy = mol_db.get(name='H2').energy
    h2o_energy = mol_db.get(name='H2O').energy
    
    # The deltaE of water formation is an empirical value
    DE_H2O = -2.97
    
    # The formation energy of the oxygen vacancy at the different sites
    o_energy_vi = energy_vi - (4 * energy_r) + h2o_energy - h2_energy - DE_H2O
    o_energy_vf = energy_vf - (4 * energy_r) + h2o_energy - h2_energy - DE_H2O
    
    return o_energy_vi, o_energy_vf
    
def get_magmom(entry: str, db: SQLite3Database):
    # Get the magnetic moment of the structure with the oxygen vacancies
    try:
        magmom_vi = db.get(name=f"{entry}_vi").magmom
        magmom_vf = db.get(name=f"{entry}_vf").magmom
    except AttributeError:
        magmom_vi = 0
        magmom_vf = 0
    return magmom_vi, magmom_vf

def format_e(n):
    a = '%.3E' % Decimal(n)
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]

def calculate_density(entry):
    # Calculate the theoretical density of the structure
    volume = entry.volume * 1E-24
    natoms = entry.natoms
    atoms = entry.toatoms()
    atoms_count = Counter(atoms.numbers)        
    normalized_count = {key: value / natoms for key, value in atoms_count.items()}
    mass = sum(data.atomic_masses[key] * value for key, value in normalized_count.items())
    navog = 6.022e23
    
    rho = natoms * mass / (volume * navog)
    return rho

def get_pettifor_number(symbol):
    element = md.element(symbol)
    return element.pettifor_number
    
# Start from the entries that have the 'barrier' key in the database.
entries = db.select(selection='barrier') 
dic_list = []
for entry in entries:
    properties = {}
    properties['neb_id'] = entry.id
    name_parts = entry.name.split('_')
    basename = '_'.join(name_parts[0:2])
    dops = int(name_parts[1][-1])
    atomsrow = entry.toatoms()
    numbers = [atom for atom in atomsrow.numbers if atom != 8]
    sort_numbers = [x[0] for x in Counter(numbers).most_common()]
    symbols = [data.chemical_symbols[i] for i in sort_numbers]
    properties['Pettifor A'] = get_pettifor_number(symbols[0])
    properties['Pettifor B'] = get_pettifor_number(symbols[1])
    properties['Pettifor dopant'] = get_pettifor_number(symbols[-1])    
    
    properties['system name'] = f"{symbols[0]}{symbols[1]}{symbols[-1]}O_p{dops}"
    properties['site A'] = symbols[0]
    properties['site B'] = symbols[1]
    properties['dopant'] = symbols[-1]
    properties['dopant site'] = dops
    properties['electroneg A'] = md.element(symbols[0]).electronegativity_pauling()
    properties['electroneg B'] = md.element(symbols[1]).electronegativity_pauling()
    properties['electroneg dopant'] = md.element(symbols[-1]).electronegativity_pauling()
    # If Yb is the dopant, the electronegativity is 1.1
    if symbols[-1] == 'Yb':
        properties['electroneg dopant'] = '1.1'
    properties['site A cov rad, pm'] = data.covalent_radii[sort_numbers[0]]
    properties['site B cov rad, pm'] = data.covalent_radii[sort_numbers[1]]
    properties['dopant cov rad, pm'] = data.covalent_radii[sort_numbers[-1]]
    
    # Add the density of the structure
    properties['density, g/cm^3'] = f"{calculate_density(entry):.3f}"
    
    # Get the relative energies with respect of the other dopants in the same system.
    energies = sorted([db.get(name=f"{name_parts[0]}_p{i}_r").energy] for i in range(1, 4))
    energy = db.get(name=f"{basename}_r").energy
    min_e = min(energies)[0]
    rel_energy_atom = (energy - min_e) / len(numbers)
    properties['rel. energy, eV/atom'] = f"{rel_energy_atom:.3f}"
    gap = db.get(name=f"{basename}_r").gap
    properties['gap, eV'] = f"{gap:.3f}"
    properties['E above hull, eV/atom'] = f"{db.get(name=f'{basename}_r').e_above_hull:.3f}"
    # Getting the oxygen vacancy formation energies
    vi_energy, vf_energy = get_oxy_vac_energy(basename, db)
    properties['O vac formation ini, eV'] = f"{vi_energy:.3f}"
    properties['O vac formation fin, eV'] = f"{vf_energy:.3f}"
    properties['Avg O vac formation, eV'] = f"{(vi_energy + vf_energy) / 2:.3f}"

    # Magnetic properties of the structures with a vacancy
    magmom_vi, magmom_vf = get_magmom(basename, db)
    properties['magmom ini'] = f"{magmom_vi:.3f}"
    properties['magmom fin'] = f"{magmom_vf:.3f}"
    properties['magmom avg'] = f"{(magmom_vi + magmom_vf) / 2:.3f}"
    
    # NEB related properties
    delta_e = entry.delta_e
    properties['dE endpoints, eV'] = f"{delta_e:.3f}"
    barrier = entry.barrier
    properties['barrier, eV'] = f"{barrier:.3f}"
    dic_list.append(properties)
    
    #break
df = pd.DataFrame(dic_list)

# reorder the entries by the system name followed by the dopant site. Reset the index.
df = df.sort_values(by=['Pettifor A', 'Pettifor B', 'Pettifor dopant', 'dopant site']).reset_index(drop=True)

df.to_csv(cwd/'hexag_perovs_schema.csv', index=False)
