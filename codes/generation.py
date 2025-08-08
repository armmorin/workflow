from functools import lru_cache
from itertools import product
from json import load
from typing import List, Tuple, Union

from ase import Atoms
from ase.db import connect
import numpy as np

from herculestools.dft import RunConfiguration


def indices(idx: int, Nb: int) -> Tuple[int, int, int]:
    i = idx // (Nb * (Nb-1))
    j = (idx // (Nb-1)) % Nb
    k = idx % (Nb-1)
    return i, j, k


def get_combination(index: int) -> Tuple[str, str, str]:
    A = ['Na', 'Mg', 'K',  'Sr', 'Ba', 'Li', 'Ca', 'Y', 'La']
    B = ['Ti', 'V',  'Zr', 'Nb', 'Mo', 'Hf', 'Ta', 'W', 'Ce', 'Sm',
         'Gd', 'Sn', 'Er', 'Sb', 'Ga', 'Ge', 'Bi', 'Al',
         'Sc', 'Mn', 'Fe', 'Pr', 'Nd', 'Zn', 'Yb', 'Si', 'Dy']

    i, j, k = indices(index, len(B))
    A_atom = A[i]
    B_atom = B[j]
    _ = B.pop(j)
    Bp_atom = B[k]
    return A_atom, B_atom, Bp_atom


@lru_cache
def load_data() -> dict:
    with open(RunConfiguration.get().home/"data/periodic_table.json") as fd:
        return load(fd)


def electrons(atoms: Atoms) -> int:
    return atoms.numbers.sum()


def get_ox(states: list, amount: int) -> list:
    states = [s for s in states if s > 0]
    return np.unique([p for p in product(states, repeat=amount)])

def get_B_ox_states(states: list, amount: int) -> List[list]:
    states = [s for s in states if s > 0]
    return np.unique(
        np.sort(
            list(product(states, repeat=amount)), axis=1
        ), axis=0
    )


def electroneutral_ox_states(A_atom: str, B_atom: str, Bp_atom: str
                         ) -> Tuple[bool, Union[dict, None]]:
    data = load_data()
    ox_states = {e: [int(i) for i in data[e]['Shannon radii'].keys()]
                 for e in [A_atom, B_atom, Bp_atom]}

    ox_A = get_ox(ox_states[A_atom], 7)[-1]
    ox_Bs = get_B_ox_states(ox_states[B_atom], 4)
    ox_Bp = get_ox(ox_states[Bp_atom], 1)
    B_site_ox_states = [
        [*b, bp] for b in ox_Bs for bp in ox_Bp
    ]
    combs = [(sum(p) + ox_A*7 - 2*20, p) for p in B_site_ox_states
             if 0 == sum(p) + ox_A*7 - 2*20]

    if not combs:
        return False, None
    if len(combs) == 1:
        c = combs[0][1]
        return True, {A_atom: ox_A, B_atom: c[:-1], Bp_atom: c[-1]}
    
    X_B = data[B_atom]['X']
    X_Bp = data[Bp_atom]['X']
    states = reversed([p for _, p in combs])

    if X_B/X_Bp <= 0.9:
        c = max(states, key=lambda p: p[-1] - np.mean(p[:-1]))
    elif X_B/X_Bp < 1.1:
        c = min(states, key=lambda p: abs(np.mean(p[:-1]) - p[-1]))
    else:
        c = max(states, key=lambda p: np.mean(p[:-1]) - p[-1])
    
    return True, {A_atom: ox_A, B_atom: c[:-1], Bp_atom: c[-1]}


def Shannon_radius(atom: str, oxidation_state: int, coordination: str = "VI"):
    data = load_data()
    dct = data[atom]['Shannon radii'][str(oxidation_state)]
    # Get the required coordination, or just any
    dct = dct.get(coordination, dct[list(dct.keys())[0]])
    return dct.get('', dct.get('High Spin'))['ionic_radius']


def goldschmidt_tolerance(oxidation_states: dict) -> float:
    A_atom, B_atom, Bp_atom, *_ = oxidation_states.keys()

    r_A = Shannon_radius(A_atom, oxidation_states[A_atom])
    r_O = Shannon_radius('O', '-2')
    r_Bp = Shannon_radius(Bp_atom, oxidation_states[Bp_atom])
    r_B = np.mean([
        Shannon_radius(B_atom, ox) for ox in oxidation_states[B_atom]
    ] + [r_Bp])
    return (r_A + r_O) / (np.sqrt(2) * (r_B + r_O) )


def decorate(template: Atoms, A: str, B: str, Bp: str, Bp_index: int
             ) -> Atoms:
    decorated = template.copy()
    decorated.symbols[decorated.symbols == 'Ba'] = A
    decorated.symbols[decorated.symbols == 'Nb'] = B
    decorated.symbols[Bp_index] = Bp
    return decorated


def _main(A: str, B: str, Bp: str):
    RunConfiguration.load()
    hex_per = Atoms(A*7 + B*4 + Bp + "O"*20)
    print("Composition:", hex_per.symbols)

    even_electrons = (electrons(hex_per) % 2 == 0)
    print("Has even number of electrons:", even_electrons)
    is_neutral, ox_states = electroneutral_ox_states(A, B, Bp)
    print("Is electroneutral:", is_neutral)
    if not is_neutral:
        return False, None
    print(ox_states)
    tol = goldschmidt_tolerance(ox_states)
    print("Has tolerance factor:", tol)
    is_tolerant = (tol > 0.8)


    if not (even_electrons and is_neutral and is_tolerant):
        return False, None

    db_path = RunConfiguration.get().structures_dir / "hexagonal_perovskites.db"
    from ase.db.sqlite import SQLite3Database
    db: SQLite3Database
    with connect(db_path) as db:
        template = db.get('Ba7Nb5O20').toatoms()
    
    p1, p2, p3 = (11, 8, 10)
    dop_p1 = decorate(template, A, B, Bp, p1)
    dop_p2 = decorate(template, A, B, Bp, p2)
    dop_p3 = decorate(template, A, B, Bp, p3)

    ids = {"ids": []}
    with connect(db_path) as db:
        def insert(atoms: Atoms, idx: int):
            name = f"{hex_per.symbols}_p{idx+1}"
            if db.count(name=name):
                ids['ids'].append(next(db.select(name=name)).id)
            else:
                ids["ids"].append(db.write(atoms, name=name))
        
        for i, a in enumerate([dop_p1, dop_p2, dop_p3]):
            insert(a, i)
    
    return True, ids

def main(index: int=1):
    a, b, bp = get_combination(index)
    return _main(a, b, bp)

if __name__ == '__main__':
    from sys import argv
    if len(argv) == 2:
        a, b, bp = get_combination(int(argv[1]))
    elif len(argv) == 4:
        a, b, bp = argv[1:]
    
    _main(a, b, bp)

# Get the constituents - A, B, B'
# Check rules
# Determine 'fitness' - either fail, or continue - EXIT CLAUSE

# Load the base structure
# Decorate it in the 3 positions
# Finish