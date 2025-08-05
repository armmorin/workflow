# © Copyright [2023] Technical University of Denmark, Steen Lysgaard

from PIL import Image
import os, sys
from pathlib import Path
import matplotlib.pyplot as plt
from ase.io import read, write
from ase.io.pov import write_pov
from ase import Atoms
from ase.data import covalent_radii
from ase.neb import NEBTools
import numpy as np
from typing import List

def texture_dict(texture_key:str = 'steen', texture_version:str = 'old'):
    """ 
    Returns the texture dictionary for the povray image

    Args:
        texture_key (str): Key for the texture dictionary
        texture_version (str): Version of the texture dictionary
    
    Returns:
        str: Texture string

    """
     # These new styles were an attempt to port the old styles o the correct
    # gamma, many had or still have unphysical light properties inorder to
    # acheive a certain look.
    material_styles_dict = dict(
        simple='finish {phong 0.7 ambient 0.4 diffuse 0.55}',
        # In general, 'pale' doesn't conserve energy and can look
        # strange in many cases.
        pale=('finish {ambient 0.9 diffuse 0.30 roughness 0.001 '
              'specular 0.2 }'),
        intermediate=('finish {ambient 0.4 diffuse 0.6 specular 0.1 '
                      'roughness 0.04}'),
        vmd=(
            'finish {ambient 0.2 diffuse 0.80 phong 0.25 phong_size 10.0 '
            'specular 0.2 roughness 0.1}'),
        jmol=('finish {ambient 0.4 diffuse 0.6 specular 1 roughness 0.001 '
              'metallic}'),
        ase2=('finish {ambient 0.2 brilliance 3 diffuse 0.6 metallic '
              'specular 0.7 roughness 0.04 reflection 0.15}'),
        ase3=('finish {ambient 0.4 brilliance 2 diffuse 0.6 metallic '
              'specular 1.0 roughness 0.001 reflection 0.0}'),
        glass=('finish {ambient 0.4 diffuse 0.35 specular 1.0 '
               'roughness 0.001}'),
        glass2=('finish {ambient 0.3 diffuse 0.3 specular 1.0 '
                'reflection 0.25 roughness 0.001}'),
        )

    # These styles were made when assumed_gamma was 1.0 which gives poor color
    # reproduction, the correct gamma is 2.2 for the sRGB standard.
    material_styles_dict_old = dict(
        steen = "finish {ambient 0.0 diffuse 0.65 " "specular 0.1 }",
        simple='finish {phong 0.7}',
        pale=('finish {ambient 0.5 diffuse 0.85 roughness 0.001 '
                'specular 0.200 }'),
        intermediate=('finish {ambient 0.3 diffuse 0.6 specular 0.1 '
                        'roughness 0.04}'),
        vmd=('finish {ambient 0.0 diffuse 0.65 phong 0.1 phong_size 40.0 '
                'specular 0.5 }'),
        jmol=('finish {ambient 0.2 diffuse 0.6 specular 1 roughness 0.001 '
                'metallic}'),
        ase2=('finish {ambient 0.05 brilliance 3 diffuse 0.6 metallic '
                'specular 0.7 roughness 0.04 reflection 0.15}'),
        ase3=('finish {ambient 0.15 brilliance 2 diffuse 0.6 metallic '
                'specular 1.0 roughness 0.001 reflection 0.0}'),
        glass=('finish {ambient 0.05 diffuse 0.3 specular 1.0 '
                'roughness 0.001}'),
        glass2=('finish {ambient 0.01 diffuse 0.3 specular 1.0 '
                'reflection 0.25 roughness 0.001}'),
        )
    if texture_version == 'old':
        return material_styles_dict_old[texture_key]
    else:
        return material_styles_dict[texture_key]


def write_povray(fname, atoms, bbox, rotation, transmittances, radii):
    """
    Create a povray image of the atoms object

    Args:
        fname (str): Name of the povray file
        atoms (Atoms): Atoms object
        bbox (list): Bounding box
        rotation (str): Rotation of the image
        transmittances (list): List of transmittances for each atom
        radii (list): List of radii for each atom
    """
    cwd = os.getcwd()
    os.chdir(image_path)
    pov_obj = write_pov(
        fname,
        atoms,
        bbox=bbox,
        show_unit_cell=2,
        rotation=rotation,
        povray_settings={
            "textures": [texture_dict(texture_key = 'intermediate', texture_version= 'new')]* len(atoms), #Steen seetings            
            "transmittances": transmittances,
            "transparent": True,
            "canvas_width": 800,
            "display": False,
        },
        radii=radii,
    )
    # Run povray to create the png
    png_path = pov_obj.render()
    os.chdir(cwd)


def get_atoms_transmittances_and_radii(images: List[Atoms],moving_index: int = None, reduce_radii_muliplier:float = 0.5,alpha:float = 0.1, reduce_other_atoms_radii:List[str] = None): 
    """
    Get the atoms object, transmittances and radii for the povray image

    Args:
        images (list): List of NEB images
        moving_index (int): Index of the moving atom if None the last atom is considered to be the moving atom
        reduce_radii_muliplier (float): Multiplier to reduce the radii of the adatoms
        alpha (float): Transmittance of the adatoms

    Returns:
        Atoms: Atoms object
        list: List of transmittances
        list: List of radii
    """
    # if moving_index is None, the last atom is considered to be the moving atom
    # For bulk, take the average position of the Na atoms in all the images
    Nbulk = len(images[0]) - 1
    bulk_positions = np.zeros((Nbulk, 3))

    if moving_index is None:
        index = [a.index for a in images[0]]
        moving_index = index[-1]
        index = index[:-1]
        
    else:
        index = [a.index for a in images[0] if a.index != moving_index]
    for image in images:
        bulk_positions += image.positions[index]
    bulk_positions /= len(images)
    atoms = Atoms(
        np.array(images[0].get_chemical_symbols())[index], positions=bulk_positions, cell=images[0].cell, pbc=images[0].pbc
    )
    # Add another version of the bulk atoms, shifted by one cell in the -x direction
    #ac = atoms.copy()
    #ac.translate(atoms.cell[0] * -1)
    #atoms += ac
    # Add another version of the bulk atoms, shifted by one cell in the +y direction
    #ac = atoms.copy()
    #ac.translate(atoms.cell[1])
    #atoms += ac
    transmittances = [0] * Nbulk #* 4
    
    # Reduce radii of some specfic atoms
    if reduce_other_atoms_radii:
        radii = []
        for a in atoms:
            # If the reduced atom is mentioned by atom type
            if isinstance(reduce_other_atoms_radii[0], str):        
                if a.symbol in reduce_other_atoms_radii:
                    radii.append(covalent_radii[a.number] * reduce_radii_muliplier)

                else:
                    radii.append(covalent_radii[a.number])
            # If the reduced atom is mentioned by index
            elif isinstance(reduce_other_atoms_radii[0], int):
                if a.index in reduce_other_atoms_radii:
                    radii.append(covalent_radii[a.number] * reduce_radii_muliplier)
                else:
                    radii.append(covalent_radii[a.number])
            # No reduction
            else: 
                radii.append(covalent_radii[a.number])
        radii = np.array(radii)
    else:
        radii = covalent_radii[atoms.numbers]

    # Add the Cl atoms
    for i, image in enumerate(images):
        adatoms = image[moving_index]
        # Lower the radii of all adatoms by 50%
        adradii = np.array([covalent_radii[adatoms.number] * reduce_radii_muliplier])
        atoms += adatoms
        radii = np.concatenate((radii, adradii))
        if i not in [0, len(images) - 1]:
            transmittances.extend([alpha] * 1)#len(adatoms)) # only one atom moves
        else:
            transmittances.extend([0] * 1)#len(adatoms)) # only one atom moves

    return atoms, transmittances, radii


def add_periodic_image(atoms, direction):
    ac = atoms.copy()
    ac.translate(np.dot(direction, atoms.cell))
    atoms += ac

# base_path = Path('/home/energy/stly/calcs/Salbage/Al111/barriers-with-solvation')
base_path = Path(".")
image_path = base_path / "image"

image_path.mkdir(exist_ok=True)
#calc_path = "/home/energy/mahpe/Playground/Kyoto/NEB/2_Fe24Na23O96P24/neb_last.traj"
calc_path = "/home/energy/armoma/phd/hex_perovs/discarded/NEB/Sr7V4MoO20_p1/Sr7V4MoO20_p1_5.traj"
#calc_path = "/home/energy/mahpe/Playground/Kyoto/NEB_pbe/3_Fe24NaO96P24/neb_last.traj"
#calc_path = "/home/energy/mahpe/Playground/Kyoto/NEB_pbe/0_Fe24Na16O96P24/neb_last.traj"
images = read(calc_path , index="-5:")
Ion_moving = 65
#Na_moving = 144 
#Na_moving = 22

# Remove all atoms above y=5Å
#y_limit = 5.205
#ref_image = images[3]

# finding atom index above y=y_limit
#index = []
#for i, atom in enumerate(ref_image):
#    if atom.position[1] > y_limit:
#        index.append(i)
# removing the atoms above y=5Å
#for image in images:
#    del image[index]

print('Creating povray images')
atoms, transmittances, radii = get_atoms_transmittances_and_radii(images, moving_index=Ion_moving, reduce_radii_muliplier=0.85, alpha=0.5,
                                                                  reduce_other_atoms_radii=[])
bbox = None #[-2, 23, 5, 21]
name = 'Sr7V4MoO20_p1'
#name = 'FePO4'
#name = 'Na0.66FePO4'

rotation = "-180x,0y,90z"
#rotation = "90x,0y,0z"
print(f'Writing povray images {name}.pov')
write_povray(f"{name}.pov", atoms, bbox, rotation, transmittances, radii)

