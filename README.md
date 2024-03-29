# Protein Chain Generator for Periodic Boundary Conditions

![SAPGen](SAPGen.jpg?raw=true "SAPGen")

## Technology Highlights:
SAPGenPBC is a Python-based project that employs backtracking, vectorization, and periodic boundary conditions to generate non-overlapping protein chains within a periodic box up to a target density - Outputing in the PDB format. Backtracking systematically explores alternative conformations when clashes between atoms are detected. By going back to previous residues and trying different angles, the code is able to balance between global and local search. It effectively leverages NumPy's capabilities for efficient computation and offers a modular design for easy maintenance and adaptability. 
The automation via bash/shell scripts accelerates convergence via running N generator jobs in parallel, continuously monitoring their progress, and stopping the remaining jobs when one finishes. 

### If you use SAPGenPBC in your work, please cite the following paper:

Materzok, T.; Canestraight, A.; Gorb, S.; Müller-Plathe, F. ["How Does Gecko Keratin Stick to Hydrophilic and Hydrophobic Surfaces in the Presence and Absence of Water? An Atomistic Molecular Dynamics Investigation"](https://doi.org/10.1021/acsnano.2c08627). ACS Nano 2022, 16 (11), 19261–19270.


## Why a Random Walker Can Be Used for Some Protein Sequences

In the context of the gecko keratin study, the Ge-cprp-9 protein, particularly its head and tail regions, contain large intrinsically disordered regions (IDRs) that fold into random coils. IDRs, as well as intrinsically disordered proteins (IDPs), have a flat energy landscape in contrast to the typical, strongly funneled energy landscapes of functional proteins. The result of a flat energy landscape in IDPs is that they do not adopt a favored three-dimensional (3D) conformation. Since IDPS are coil-like and disordered, their initial three-dimensional structure can be created using a self-avoiding random walker.

## Overview

### SAPGenPBC.py

This python tool generates a specified number of protein chains in a periodic box based on the provided protein sequence and density value. The protein chains are generated with random coordinates and dihedral angles, and without overlaps (self-avoiding) according to the periodic boundary conditions.

### run_generators_and_concatenate_pdbs.sh

This script starts 48 generator jobs in parallel, continuously checks for any generator that has finished, and stops the others when one finishes. It then concatenates the PDB files in the completed generator directory to create a final PDB file.

### concatenate_pdbs.sh

This script concatenates multiple PDB files into a single PDB file with the given box dimensions.

## Inputs

- Number of cores to run the generator in parallel on
- Protein sequence: a string containing the protein sequence (one-letter code)
- Box dimensions: float values for the box dimensions in Angstroms (xmax, ymax, zmax)
- Density: float value for the protein density in g/cm^3

## Outputs

- Protein chains in PDB format files named ChainX.pdb, where X is the chain number
- Concatenated PDB files using the provided shell scripts

## Dependencies

- Python 3.x
- NumPy
- BioPython
- PeptideBuilder

## Usage

1. Download the repository: The files assume that this command is run in the home directory

```
cd ~
git clone https://github.com/TobiasMaterzok/PDB-Protein-Chain-Generator-in-Periodic-Boundary-Conditions
```

2. Ensure that you have the required Python packages installed:

```
pip install numpy biopython PeptideBuilder
```

3. To run the generator embarrassingly parallel and find a solution fast, navigate to an empty directory and run

```
~/PDB-Protein-Chain-Generator-in-Periodic-Boundary-Conditions/run_generators_and_concatenate_pdbs.sh CORE_NUMBER SEQUENCE BOX_X BOX_Y BOX_Z DENSITY
```

Example (7x7x7 nm at 1.3 g cm^-3 density on 48 cores):
```
~/PDB-Protein-Chain-Generator-in-Periodic-Boundary-Conditions/run_generators_and_concatenate_pdbs.sh 48 MSCCPPSCA 70 70 70 1.3
```

3.B. You can also run the generator with the following command line arguments:

```
python SAPGenPBC.py <core_number> <sequence> <xmax> <ymax> <zmax> <dens>
```

Example:

```
python SAPGenPBC.py 48 MSCCPPSCA 70 70 70 1.3
```

Note: This script may take some time to run, depending on the input parameters. It is recommended to use the run_generators_and_concatenate_pdbs.sh wrapper.

## License

This project is licensed under the GNU General Public License v3.0 - see the LICENSE file for details.    
