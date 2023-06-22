#!/bin/bash

#
# Author: Tobias Materzok https://github.com/TobiasMaterzok/
#
# This script concatenates multiple PDB files into a single PDB file with the 
# given box dimensions. The resulting file has a new model number and an 
# "ENDMDL" statement at the end to signify the end of the model.
#
# Usage: ./concatenate_pdb.sh NUMBER_OF_CHAINS BOX_SIZE_X BOX_SIZE_Y BOX_SIZE_Z
#
# Input: 
# - NUMBER_OF_CHAINS: the number of PDB files to concatenate
# - BOX_SIZE_X, BOX_SIZE_Y, BOX_SIZE_Z: the dimensions of the box in which the 
#  concatenated structure will be placed
#
# Output:
# - all.pdb: the concatenated PDB file
#

# Set the number of chains to concatenate
N=$1

# Set the prefix for the PDB file names
cname="Chain"

# Set the dimensions of the box in which the concatenated structure will be placed
boxlenx=$2
boxleny=$3
boxlenz=$4

# Create the header for the concatenated PDB file and write it to the all.pdb file
echo "CRYST1  "$boxlenx".000  "$boxleny".000  "$boxlenz".000  90.00  90.00  90.00 P 1           1" > all.pdb

# Create a new model for the concatenated PDB file
echo "MODEL        1" >> all.pdb

# Iterate over the PDB files to be concatenated
for i in `seq 1 1 $N`
do
	# Remove the last line (which is usually "END") from each PDB file and append the remaining lines to the all.pdb file
	head -n -1 "$cname""$i".pdb >> all.pdb
done

# Add an "ENDMDL" statement to the end of the concatenated PDB file to signify the end of the model
echo "ENDMDL" >> all.pdb
