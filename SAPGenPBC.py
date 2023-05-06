"""
Self-Avoiding Protein Chain Generator in Periodic Boundary Conditions (SAPGenPBC)

Author: Tobias Materzok https://github.com/TobiasMaterzok/
        Annabelle Canestraight

This script generates a number of non-overlapping protein chains in a periodic box based on the provided
protein sequence, the box dimensions and the target density value. The protein chains are 
generated with random coordinates and dihedral angles, and their positions are adjusted
according to the periodic boundary conditions to not overlap with each other.

Inputs:

    Protein sequence (sequence): a string containing the protein sequence (one-letter code)
    Box dimensions (xmax, ymax, zmax): float values for the box dimensions in Angstroms
    Density (dens): float value for the protein density in g/cm^3

Outputs:

    Protein chains in PDB format files named ChainX.pdb, where X is the chain number

Usage:

    Ensure that you have the required Python packages installed: numpy, BioPython, PeptideBuilder
    Run the script with the following command line arguments:
    python <script_name> <sequence> <xmax> <ymax> <zmax> <dens>
    Example: python protein_chain_generator.py MSCCPPSCA 50 50 50 1.3

Note: This script may take some time to run, depending on the input parameters.
For parallelization and automation, use the accompanying script
run_generators_and_concatenate_pdbs.sh, which starts multiple generator jobs, monitors progress,
and concatenates PDB files upon completion.
"""

import sys 
import random
import numpy as np
import Bio.PDB
import PeptideBuilder
from PeptideBuilder import Geometry

#sequence= 'MSCCPPSCATPSCPKPCCSPCCSPCGYPTGGLGSLGCCPCPCGPSSCCGSSTSARCLGITSGASVSCINQIPASCEPLRVGGYTACGGCPPCGRIC'

sequence=str(sys.argv[1])
xmax=float(sys.argv[2]) # Angstrom
ymax=float(sys.argv[3]) # Angstrom
zmax=float(sys.argv[4]) # Angstrom
dens=float(sys.argv[5])


inter_cutoff = 1.9
intra_cutoff = 1.85
NA=6.02214076*10**23 # mol-1

# Define residue lengths and molar mass of amino acids
residue_length_dict={'M': 8,'S': 6,'C': 6,'P': 7,'A': 5,'T': 7,'K': 9,'G': 4,'Y': 12,
 'L': 8,'R': 11,'I': 8,'V': 7,'N': 8,'Q': 9,'E': 9,'D': 8}
molarmass={'M': 149.2,'S': 105.1,'C': 121.2,'P': 115.1,'A': 89.1,'T': 119.1,'K': 146.2,'G': 75.1,
 'Y': 181.2,'L': 131.2,'R': 174.2,'I': 131.2,'V': 117.1,'N': 132.1,'Q': 146.2,'E': 147.1,'D': 133.1}


# Calculate the total molar mass of the protein sequence and convert it to mass in grams
mass = sum(molarmass[aa] for aa in sequence) / NA

# Estimate the *natural* number of protein chains needed to achieve the desired mass density
while True:
    volume = xmax * ymax * zmax
    # Convert the volume to cubic centimeters
    volume_in_cm = volume * (10**(-10))**3 * (10**(2))**3

    # Calculate the target mass of protein chains within the box using the supplied density
    totalmass = dens * volume_in_cm

    # Estimate the number of protein chains needed to achieve the desired mass density
    ngenchains = totalmass / mass

    # Calculate the mean squared deviation between the actual number of chains
    # and the rounded integer value of the estimated number of chains
    msd = (1 - round(ngenchains) / ngenchains) ** 2

    # The mean squared deviation (msd) is used to quantify the accuracy of the
    # approximation of the number of protein chains required to achieve the
    # desired mass density. A smaller msd indicates that the rounded integer
    # value of the estimated number of chains is closer to the actual value.
    # This ensures that the generated protein chains in the simulation box are
    # close to the desired mass density.

    # Check if the mean squared deviation is smaller than or equal to the threshold value
    # The threshold value of 3e-07 is chosen as an acceptable level of accuracy for the
    # approximation. This value can be adjusted based on the desired level of accuracy
    # for the initial model creation. A smaller value would result in a more accurate representation
    # of the desired mass density but may require more iterations to reach the target.

    print(ngenchains, msd)

    # Check if the mean squared deviation is smaller than or equal to the threshold value
    if msd <= 3e-07:
        # Set the final number of protein chains and print the results
        ngenchains = int(round(ngenchains))
        print(ngenchains, msd, xmax, ymax)

        # Break since the desired mass density has been achieved
        break

    # If the mean squared deviation is larger than the threshold value,
    # increase the dimensions of the simulation box by a small factor (1.00001)
    xmax *= 1.00001
    ymax *= 1.00001

def Periodic(coordz):
    box_dimensions = np.array([xmax, ymax, zmax])
    pbc_coordz = np.mod(coordz, box_dimensions)
    return pbc_coordz

def get_rnd_coord(box_vals):
    xlow, xhigh, ylow, yhigh, zlow, zhigh = box_vals
    return np.array([np.random.uniform(xlow, xhigh), np.random.uniform(ylow, yhigh), np.random.uniform(zlow, zhigh)])

def get_rnd_degree():
    return random.uniform(-180,180)

# Function to create a partial protein chain from already existing angles and positions
def rebuild(i, angles_list, rand_coord):
    def create_geo(res_idx, angles):
        geo = Geometry.geometry(sequence[res_idx])
        geo.psi_im1, geo.omega, geo.phi = angles
        return geo

    # Initialize the first residue of the protein chain with its random coordinates
    geo = Geometry.geometry(sequence[0])
    first_angles = angles_list[0]
    geo.psi_im1, geo.omega, geo.phi = first_angles
    Structure = PeptideBuilder.initialize_res(geo)
    chain = Structure[0]['A']
    residue = chain[1]
    active_chain = []

    for atom in residue:
        new_coord = atom.get_coord() + rand_coord
        atom.set_coord(new_coord)

    # Add the remaining residues in the partial protein chain
    for k in range(1, i):
        geo = create_geo(k, angles_list[k])
        Structure = PeptideBuilder.add_residue(Structure, geo)

    # Apply periodic boundary conditions to the coordinates
    for res in chain:
        for atom in res:
            cor = atom.get_coord()
            ncor = np.array(cor)
            active_chain.append(Periodic(ncor))

    return Structure, active_chain

# Function to generate a new protein chain within the box
def Add_chain(sequence, box_vals):
    angles_list, active_chain = [], []
    rand_coord = get_rnd_coord(box_vals)

    def generate_random_angles(geo):
        geo.psi_im1, geo.omega, geo.phi = get_rnd_degree(), get_rnd_degree(), get_rnd_degree()

    # Initialize the first residue of the protein chain with random angles and random coordinates
    geo = Geometry.geometry(sequence[0])
    generate_random_angles(geo)
    angles_list.append([geo.psi_im1, geo.omega, geo.phi])
    Structure = PeptideBuilder.initialize_res(geo)

    chain = Structure[0]['A']
    residue = chain[1]
    for atom in residue:
        new_coord = atom.get_coord() + rand_coord
        atom.set_coord(new_coord)
        active_chain.append(Periodic(np.array(new_coord)))

    i, tries, iteration, backtrack = 1, 0, 0, 0
    maxiter = 1000000
    # Iterate through the remaining residues in the sequence
    while i < len(sequence) and iteration < maxiter:
        geo = Geometry.geometry(sequence[i])
        generate_random_angles(geo)

        # Add the residue with random angles
        Structure = PeptideBuilder.add_residue(Structure, geo)
        chain = Structure[0]['A']
        residue, ver_list = chain[i+1], []
        length = len(active_chain)
        iteration += 1

        pbc = []
        for atom in residue:
            # Apply periodic boundary conditions to the coordinates of the newly added residue
            coords = np.array(atom.get_coord())
            pbc_coords = Periodic(coords)
            pbc.append(pbc_coords)

        for pbc_coords in np.array(pbc):
            # Check if any atoms in the new residue are too close to atoms in earlier residues 
            # within the same chain or in other chains
            overlap = 0
            if overlap == 0:
                intra_chain_coords = np.array(active_chain[-(residue_length_dict[sequence[i-1]]-3)])
                intra_chain_distances = np.linalg.norm(pbc_coords - intra_chain_coords[:, np.newaxis], axis=1)
                if np.any(intra_chain_distances < intra_cutoff):
                        ver_list.append(False)
                        overlap = 1
            if overlap == 0:
                early_intra_chain_coords = np.array(active_chain[:(length-residue_length_dict[sequence[i-1]]+2)])
                early_intra_chain_distances = np.linalg.norm(pbc_coords - early_intra_chain_coords[:, np.newaxis], axis=-1)
                if np.any(early_intra_chain_distances < intra_cutoff):
                    ver_list.append(False)
                    overlap = 1
            if overlap == 0:
                inter_chain_coords = np.array(all_chains_position_list)
                inter_chain_distances = np.linalg.norm(pbc_coords - inter_chain_coords[:, np.newaxis], axis=-1)
                if np.any(inter_chain_distances < inter_cutoff):
                    ver_list.append(False)
                    overlap = 1

        max_backtrack = 10
        def backtrack_n_residues(n):
            nonlocal Structure, active_chain, angles_list, i
            Structure, active_chain = rebuild(i - n, angles_list, rand_coord)
            i -= n
            angles_list = angles_list[:i]        

        # If the new residue clashes with other atoms, try different positions or backtrack to a previous residue
        if False in ver_list:
            tries += 1
            if tries <= 25:
                Structure, active_chain = rebuild(i, angles_list, rand_coord)
                angles_list=angles_list[:i]
            elif tries > 25 and tries < 100 and i <= 5:
                rand_coord = get_rnd_coord(box_vals)
                Structure, active_chain = rebuild(0, angles_list, rand_coord)
                tries=0
                i=1
                angles_list=angles_list[:i] 
            elif tries > 25 and tries <= 100 and i > backtrack and backtrack <= max_backtrack:
                # Backtracking: systematically go back different numbers of residues
                backtrack_n_residues(backtrack)
                tries = 0
                backtrack += 1
            else:
                random.seed(m-2*i)
                rand_coord = get_rnd_coord(box_vals)
                Structure, active_chain = rebuild(0, angles_list, rand_coord)
                tries = 0
                i = 1
                angles_list = angles_list[:i]
                backtrack = 0
        else:
            # If the new residue does not clash with other atoms, move on to the next residue in the sequence
            active_chain.extend(pbc)
            i += 1
            angles_list.append([geo.psi_im1, geo.omega, geo.phi])
            backtrack = 0
    return active_chain, Structure

all_chains_position_list=[]
for m in range(1, ngenchains + 1):
    # Generate a new protein chain and save its structure to a PDB file
    pdb_output = Bio.PDB.PDBIO()
    pdb_file_name = str("Chain" + str(m) + ".pdb")
    boxvals = [0,xmax,0,ymax,0,zmax]

    new_peptide_coordinates, peptide_structure_new_chain = Add_chain(sequence, boxvals)

    # Add the coordinates of the new protein chain to the list of all chain coordinates
    pdb_output.set_structure(peptide_structure_new_chain)
    pdb_output.save(pdb_file_name)
    for b in new_peptide_coordinates:
        all_chains_position_list.append(b)
    print("Done", m)
