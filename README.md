# Protein Chain Generator for Periodic Boundary Conditions

![SAPGen](SAPGen.jpg?raw=true "SAPGen")

## Technology Highlights:
SAPGenPBC is a Python-based project that employs backtracking, vectorization, and periodic boundary conditions to generate non-overlapping protein chains within a periodic box up to a target density using the PDB format. Backtracking systematically explores alternative conformations when clashes between atoms are detected. By going back to previous residues and trying different angles, the script balances between global and local search. It effectively leverages NumPy's capabilities for efficient computation and offers a modular design for easy maintenance and adaptability. 
The automation via bash/shell scripts accelerates convergence via running N generator jobs in parallel, continuously monitoring their progress, and stopping the remaining jobs when one finishes. 

This project was necessary to achieve the interdisciplinary approach in ["How Does Gecko Keratin Stick to Hydrophilic and Hydrophobic Surfaces in the Presence and Absence of Water? An Atomistic Molecular Dynamics Investigation"](https://pubs.acs.org/doi/full/10.1021/acsnano.2c08627).

## Why a Random Walker Can Be Used for Some Protein Sequences

In the context of the gecko keratin study, the Ge-cprp-9 protein, particularly its head and tail regions, contain large intrinsically disordered regions (IDRs) that fold into random coils. IDRs, as well as intrinsically disordered proteins (IDPs), have a flat energy landscape in contrast to the typical, strongly funneled energy landscapes of functional proteins. The result of a flat energy landscape in IDPs is that they do not adopt a favored three-dimensional (3D) conformation. Since IDPS are coil-like and disordered, their initial three-dimensional structure can be created using a self-avoiding random walker.

## Overview

### SAPGenPBC.py

This python tool generates a specified number of protein chains in a periodic box based on the provided protein sequence and density value. The protein chains are generated with random coordinates and dihedral angles, and their positions are adjusted according to the periodic boundary conditions.

### run_generators_and_concatenate_pdbs.sh

This script starts 48 generator jobs in parallel, continuously checks for any generator that has finished, and stops the others when one finishes. It then concatenates the PDB files in the completed generator directory to create a final PDB file.

### concatenate_pdbs.sh

This script concatenates multiple PDB files into a single PDB file with the given box dimensions.

## Inputs

- Protein sequence: a string containing the protein sequence (one-letter code)
- Box dimensions: float values for the box dimensions in Angstroms (xmax, ymax, zmax)
- Density: float value for the protein density in g/cm^3

## Outputs

- Protein chains in PDB format files named ChainX.pdb, where X is the chain number
- Concatenated PDB files using the provided shell scripts

## Dependencies

- Python 3.x
- numpy
- BioPython
- PeptideBuilder

## Usage

1. Ensure that you have the required Python packages installed:

```
pip install numpy biopython PeptideBuilder
```

2. Copy the files into ~/tools_ua_gecko/

3. To run the generator embarrassingly parallel, navigate to an empty directory and run

```
~/tools_ua_gecko/run_generators_and_concatenate_pdbs.sh SEQUENCE BOX_X BOX_Y BOX_Z DENSITY
```

Example (7x7x7 nm at 1.3 g cm^-3 density):
```
~/tools_ua_gecko/run_generators_and_concatenate_pdbs.sh MSCCPPSCA 70 70 70 1.3
```

1.B. Run the generator with the following command line arguments:

```
python SAPGenPBC.py <sequence> <xmax> <ymax> <zmax> <dens>
```

Example:

```
python SAPGenPBC.py MSCCPPSCA 70 70 70 1.3
```

Note: This script may take some time to run, depending on the input parameters.

## License

This project is licensed under the MIT License - see the LICENSE file for details.    
