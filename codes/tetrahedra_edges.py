# tetrahedra_edges.py

# Import necessary libraries
from ase.io import read as ase_read
from pymatgen.io.cif import CifParser
from pymatgen.io.ase import AseAtomsAdaptor
from ase.data import atomic_numbers
from pymatgen.analysis.local_env import CrystalNN
import numpy as np
import warnings

def read_structure(filename):
    """
    Read a structure file (POSCAR or CIF) and return an ASE Atoms object.
    
    Args:
    filename (str): Path to the structure file.
    
    Returns:
    ase.Atoms: ASE Atoms object representing the structure.
    """
    file_extension = filename.split('.')[-1].lower()
    
    if file_extension in ['poscar', 'contcar'] or 'poscar' in filename.lower():
        return ase_read(filename)
    elif file_extension == 'cif':
        structure = CifParser(filename).parse_structures(primitive=True)[0]
        return AseAtomsAdaptor.get_atoms(structure)
    else:
        raise ValueError(f"Unsupported file format: {file_extension}. Please use POSCAR or CIF.")

def identify_transition_metals(structure):
    """
    Identify transition metal atoms in the given structure.
    
    Args:
    structure (ase.Atoms): ASE Atoms object representing the structure.
    
    Returns:
    list: Indices of transition metal atoms in the structure.
    """
    transition_metals = []
    for atom in structure:
        if 21 <= atomic_numbers[atom.symbol] <= 30 or \
           39 <= atomic_numbers[atom.symbol] <= 48 or \
           57 <= atomic_numbers[atom.symbol] <= 80 or \
           89 <= atomic_numbers[atom.symbol] <= 112:
            transition_metals.append(atom.index)
    return transition_metals

def find_nearest_neighbors(structure, atom_index):
    """
    Find nearest neighbors of a given atom using pymatgen's CrystalNN.
    
    Args:
    structure (ase.Atoms): ASE Atoms object representing the structure.
    atom_index (int): Index of the atom to find neighbors for.
    
    Returns:
    list: Indices of nearest neighbor atoms.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pymatgen_structure = AseAtomsAdaptor.get_structure(structure)
        cnn = CrystalNN()  # CrystalNN inherently considers PBC for crystalline structures
        neighbors = cnn.get_nn_info(pymatgen_structure, atom_index)
        neighbor_indices = [neighbor['site_index'] for neighbor in neighbors]
        return neighbor_indices

def calculate_alignment_score(angle):
    """
    Calculate the alignment score for a given angle.
    
    Args:
    angle (float): Angle in degrees.
    
    Returns:
    float: Alignment score (lower is better aligned).
    """
    return min(angle, abs(180 - angle))

def analyze_coordination(structure, atom_indices):
    """
    Analyze the coordination environment of specified atoms.
    
    Args:
    structure (ase.Atoms): ASE Atoms object representing the structure.
    atom_indices (list): List of atom indices to analyze.
    
    Returns:
    dict: Coordination environments for each specified atom.
    """
    coordination_environments = {}
    for atom_index in atom_indices:
        neighbors = find_nearest_neighbors(structure, atom_index)
        coordination_environments[atom_index] = neighbors
    return coordination_environments

def check_symmetry(structure, atom_indices):
    """
    Check symmetry of coordination polyhedra for specified atoms.
    
    Args:
    structure (ase.Atoms): ASE Atoms object representing the structure.
    atom_indices (list): List of atom indices to analyze.
    
    Returns:
    dict: Symmetry information for each specified atom.
    """
    symmetries = {}
    pymatgen_structure = AseAtomsAdaptor.get_structure(structure)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cnn = CrystalNN()
        for atom_index in atom_indices:
            site = pymatgen_structure[atom_index]
            neighbors = cnn.get_nn_info(pymatgen_structure, atom_index)
            symmetries[atom_index] = {
                'type': structure[atom_index].symbol,
                'coordination_number': cnn.get_cn(pymatgen_structure, atom_index)
            }
    return symmetries

def get_atom_info(structure, atom_index):
    """
    Get information about a specific atom in the structure.
    
    Args:
    structure (ase.Atoms): ASE Atoms object representing the structure.
    atom_index (int): Index of the atom to get information for.
    
    Returns:
    dict: Information about the specified atom.
    """
    atom = structure[atom_index]
    return {
        'index': atom_index,
        'symbol': atom.symbol,
        'coordinates': atom.position,
    }

def calculate_angle(v1, v2):
    """
    Calculate the angle between two vectors.
    
    Args:
    v1, v2 (numpy.array): Vectors to calculate angle between.
    
    Returns:
    float: Angle in degrees.
    """
    cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    angle_rad = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    angle_deg = np.degrees(angle_rad)
    return angle_deg

def adjust_position_with_pbc(pos, ref_pos, cell):
    """
    Adjust position vector considering periodic boundary conditions.
    
    Args:
    pos (numpy.array): Position to adjust.
    ref_pos (numpy.array): Reference position.
    cell (numpy.array): Unit cell of the structure.
    
    Returns:
    numpy.array: Adjusted position vector.
    """
    diff_vector = pos - ref_pos
    diff_vector -= np.dot(np.round(np.dot(diff_vector, np.linalg.inv(cell))), cell)
    adjusted_pos = ref_pos + diff_vector
    return adjusted_pos

def calculate_diagonals(structure, atom_indices):
    """
    Calculate diagonal lengths of octahedra with periodic boundary conditions.
    
    Args:
    structure (ase.Atoms): ASE Atoms object representing the structure.
    atom_indices (list): List of atom indices to analyze.
    
    Returns:
    dict: Diagonal lengths for each specified atom.
    """
    diagonals = {}
    cell = structure.get_cell()
    lattice_vectors = [cell[0], cell[1], cell[2]]
    
    for atom_index in atom_indices:
        neighbors = find_nearest_neighbors(structure, atom_index)
        if len(neighbors) == 6:  # Ensure it's an octahedron
            central_position = structure[atom_index].position
            neighbor_positions = [adjust_position_with_pbc(structure[neighbor].position, central_position, cell) for neighbor in neighbors]
            diagonals[atom_index] = [None, None, None]
            
            for i in range(3):  # For each lattice vector a, b, c
                best_angle_diff = float('inf')
                best_diagonal_length = None
                
                for j in range(6):
                    for k in range(j+1, 6):
                        diagonal_length = np.linalg.norm(neighbor_positions[j] - neighbor_positions[k])
                        diff_vector = neighbor_positions[j] - neighbor_positions[k]
                        angle_deg = calculate_angle(diff_vector, lattice_vectors[i])
                        angle_diff_to_0_or_180 = min(abs(angle_deg), abs(180 - angle_deg))
                        
                        print(f"Diagonal between atoms {j} and {k}: length = {diagonal_length}, angle with lattice vector {i} = {angle_deg} degrees")
                        
                        if angle_diff_to_0_or_180 < best_angle_diff:
                            best_angle_diff = angle_diff_to_0_or_180
                            best_diagonal_length = diagonal_length
                
                print(f"Selected diagonal for lattice vector {i}: length = {best_diagonal_length}, angle difference to 0 or 180 degrees = {best_angle_diff}")
                diagonals[atom_index][i] = best_diagonal_length
    
    return diagonals

def interpret_vector_input(structure, vector_input):
    """
    Interpret vector input as either a lattice vector or a custom vector.
    
    Args:
    structure (ase.Atoms): ASE Atoms object representing the structure.
    vector_input (str or list): Input vector as either 'a', 'b', 'c' or a list of 3 numbers.
    
    Returns:
    numpy.array: Interpreted vector.
    """
    if isinstance(vector_input, str):
        vector_input = vector_input.lower()
        if vector_input == 'a':
            return structure.cell[0]
        elif vector_input == 'b':
            return structure.cell[1]
        elif vector_input == 'c':
            return structure.cell[2]
        elif vector_input == 'd':
            return structure.cell[0] + structure.cell[1]           
            
        else:
            raise ValueError("Invalid string input. Use 'a', 'b', or 'c' for lattice vectors.")
    elif isinstance(vector_input, (list, np.ndarray)):
        return np.array(vector_input)
    else:
        raise ValueError("Vector input must be a string ('a', 'b', 'c') or a list/array of 3 numbers.")

def analyze_coordination_and_edges(structure, atom_index, alignment_vector):
    """
    Analyze the coordination environment and edge alignments for a given atom.
    
    Args:
    structure (ase.Atoms): ASE Atoms object representing the structure.
    atom_index (int): Index of the atom to analyze.
    alignment_vector (numpy.array): Vector to align edges with.
    
    Returns:
    tuple: (neighbors, edges) where neighbors is a list of neighbor indices and
           edges is a list of tuples (neighbor_index, edge_length, angle_with_alignment, alignment_score)
           sorted by alignment_score (lower is better aligned).
    """
    neighbors = find_nearest_neighbors(structure, atom_index)
    central_position = structure[atom_index].position
    cell = structure.get_cell()
    
    edges = []
    for neighbor in neighbors:
        neighbor_position = adjust_position_with_pbc(structure[neighbor].position, central_position, cell)
        edge_vector = neighbor_position - central_position
        edge_length = np.linalg.norm(edge_vector)
        angle_with_alignment = calculate_angle(edge_vector, alignment_vector)
        alignment_score = calculate_alignment_score(angle_with_alignment)
        edges.append((neighbor, edge_length, angle_with_alignment, alignment_score))
    
    # Sort edges by alignment score (smaller is better)
    edges.sort(key=lambda x: x[3])
    
    return neighbors, edges

def calculate_angle_between_edges(structure1, index1, vector1, structure2, index2, vector2):
    """
    Calculate the angle between the most aligned edges of two different structures.

    Args:
    structure1, structure2 (ase.Atoms): ASE Atoms objects representing the crystal structures.
    index1, index2 (int): Indices of the central atoms in each structure.
    vector1, vector2 (str or list): Alignment vectors for each structure (as strings or numpy arrays).

    Returns:
    tuple: (angle, edge1, edge2) where angle is the angle (in degrees) between the two most aligned edges,
           and edge1, edge2 are the most aligned edge information for each structure.
    """
    # Process the first structure
    alignment_vector1 = interpret_vector_input(structure1, vector1)
    neighbors1, edges1 = analyze_coordination_and_edges(structure1, index1, alignment_vector1)
    edge1 = edges1[0]  # Most aligned edge for structure1

    # Process the second structure
    alignment_vector2 = interpret_vector_input(structure2, vector2)
    neighbors2, edges2 = analyze_coordination_and_edges(structure2, index2, alignment_vector2)
    edge2 = edges2[0]  # Most aligned edge for structure2

    # Calculate edge vectors
    central_position1 = structure1[index1].position
    neighbor_position1 = adjust_position_with_pbc(structure1[edge1[0]].position, central_position1, structure1.get_cell())
    edge_vector1 = neighbor_position1 - central_position1

    central_position2 = structure2[index2].position
    neighbor_position2 = adjust_position_with_pbc(structure2[edge2[0]].position, central_position2, structure2.get_cell())
    edge_vector2 = neighbor_position2 - central_position2

    # Calculate angle between edge vectors
    angle = calculate_angle(edge_vector1, edge_vector2)

    return angle, edge1, edge2, edge_vector1, edge_vector2, neighbor_position1, neighbor_position2

def main(filename, metal_index, alignment_vector_input, print_edges=True):
    """
    Main function to analyze a structure file and calculate edge alignments.
    
    Args:
    filename (str): Path to the structure file (POSCAR or CIF).
    metal_index (int): Index of the metal atom to analyze.
    alignment_vector_input (str or list): Alignment vector input.
    print_edges (bool): Whether to print detailed edge information.
    
    Returns:
    tuple: Information about the most aligned edge.
    """
    structure = read_structure(filename)
    
    atom_info = get_atom_info(structure, metal_index)
    print(f"\nAtom Information:")
    print(f"Index: {atom_info['index']}")
    print(f"Element: {atom_info['symbol']}")
    print(f"Coordinates: {atom_info['coordinates']}")
    
    alignment_vector = interpret_vector_input(structure, alignment_vector_input)
    
    neighbors, edges = analyze_coordination_and_edges(structure, metal_index, alignment_vector)
    
    print(f"\nFile: {filename}")
    print(f"Coordination number for atom {metal_index}: {len(neighbors)}")
    print(f"Alignment vector: {alignment_vector}")
    
    if print_edges:
        print("\nEdges:")
        for edge in edges:
            print(f"Edge to atom {edge[0]}: length = {edge[1]:.3f}, angle with alignment vector = {edge[2]:.2f} degrees, alignment score = {edge[3]:.2f}")
    
    return edges[0]

# Example usage
if __name__ == "__main__":
    # Example with a POSCAR file
    #poscar_file = 'Ba7Nb4MoO20_p1_ts_translated.cif' #'POSCAR_ts'
    #metal_index_poscar = 44
    # poscar_file = 'Ba7Nb4MoO20_p1_r_translated.cif' #'POSCAR_r'
    # metal_index_poscar = 11
    # result_poscar = main(poscar_file, metal_index_poscar, 'c', print_edges=True)
    # print(f"\nPOSCAR - Most aligned edge with 'c' vector: to atom {result_poscar[0]}, length = {result_poscar[1]:.3f}, angle = {result_poscar[2]:.2f} degrees, alignment score = {result_poscar[3]:.2f}")

    # Example usage of the angle calculation between edges from different structures
    structure1 = read_structure('Ba7V4MoO20_p3_r.cif') # 'POSCAR_r'
    structure2 = read_structure('Ba7V4MoO20_p3_ts.cif') # 'POSCAR_ts'
    index1 = 7 #11  # Example atom index for structure1
    index2 = 28 #44  # Example atom index for structure2
    vector1 = 'c'  # Alignment vector for structure1
    vector2 = 'c'  # Alignment vector for structure2

    angle, edge1, edge2 = calculate_angle_between_edges(structure1, index1, vector1, structure2, index2, vector2)

    print(f"\nAngle between most aligned edges: {angle:.2f} degrees")
    print(f"Edge 1: to atom {edge1[0]}, length = {edge1[1]:.3f}, angle with alignment vector = {edge1[2]:.2f} degrees, alignment score = {edge1[3]:.2f}")
    print(f"Edge 2: to atom {edge2[0]}, length = {edge2[1]:.3f}, angle with alignment vector = {edge2[2]:.2f} degrees, alignment score = {edge2[3]:.2f}")
    
    index1 = 10 #11  # Example atom index for structure1
    index2 = 30 #44  # Example atom index for structure2
    vector1 = 'c'  # Alignment vector for structure1
    vector2 = 'c'  # Alignment vector for structure2

    angle, edge1, edge2 = calculate_angle_between_edges(structure1, index1, vector1, structure2, index2, vector2)

    print(f"\nAngle between most aligned edges: {angle:.2f} degrees")
    print(f"Edge 1: to atom {edge1[0]}, length = {edge1[1]:.3f}, angle with alignment vector = {edge1[2]:.2f} degrees, alignment score = {edge1[3]:.2f}")
    print(f"Edge 2: to atom {edge2[0]}, length = {edge2[1]:.3f}, angle with alignment vector = {edge2[2]:.2f} degrees, alignment score = {edge2[3]:.2f}")