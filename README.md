# Protein Chain Generator for Periodic Boundary Conditions

## Technology Highlights:
SAPGenPBC is a Python-based project that employs backtracking, vectorization, and periodic boundary conditions to generate non-overlapping protein chains within a periodic box up to a target density using the PDB format. It effectively leverages NumPy's capabilities for efficient computation and offers a modular design for easy maintenance and adaptability. This project was necessary to achieve the interdisciplinary approach in ["How Does Gecko Keratin Stick to Hydrophilic and Hydrophobic Surfaces in the Presence and Absence of Water? An Atomistic Molecular Dynamics Investigation"](https://pubs.acs.org/doi/full/10.1021/acsnano.2c08627).

Backtracking systematically explores alternative conformations when clashes between atoms are detected. By going back to previous residues and trying different angles, the script balances between global and local search.

## Why a Random Walker Can Be Used for Some Protein Sequences

In the context of the gecko keratin study, the Ge-cprp-9 protein, particularly its head and tail regions, contain large intrinsically disordered regions (IDRs) that fold into random coils. IDRs, as well as intrinsically disordered proteins (IDPs), have a flat energy landscape in contrast to the typical, strongly funneled energy landscapes of functional proteins. The result of a flat energy landscape in IDPs is that they do not adopt a favored three-dimensional (3D) conformation. Since IDPS are coil-like and disordered, their initial three-dimensional structure can be created using a self-avoiding random walker.

## Overview

The script generates a specified number of protein chains in a periodic box based on the provided protein sequence and density value. The protein chains are generated with random coordinates and dihedral angles, and their positions are adjusted according to the periodic boundary conditions.

## Inputs

- Protein sequence: a string containing the protein sequence (one-letter code)
- Box dimensions: float values for the box dimensions in Angstroms (xmax, ymax, zmax)
- Density: float value for the protein density in g/cm^3

## Outputs

- Protein chains in PDB format files named ChainX.pdb, where X is the chain number

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

2. Run the script with the following command line arguments:

```
python <script_name> <sequence> <xmax> <ymax> <zmax> <dens>
```

Example:

```
python SAPGenPBC.py MSCCPPSCA 50 50 50 1.3
```

Note: This script may take some time to run, depending on the input parameters.

## License

This project is licensed under the MIT License - see the LICENSE file for details.    
