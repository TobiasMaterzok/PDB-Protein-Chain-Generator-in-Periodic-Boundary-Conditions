#!/bin/bash

#
# Author: Tobias Materzok https://github.com/TobiasMaterzok/
#
# Description:
# This script starts 48 generator jobs in parallel. It continuously checks for
# any generator that has finished. When one generator finishes, it stops the
# others, deletes their directories, and concatenates the PDB files in the
# completed generator directory to create a final PDB file.
#
# Usage: ./run_generators_and_concatenate_pdbs.sh SEQUENCE BOX_X BOX_Y BOX_Z DENSITY
#
# Input:
# - SEQUENCE: The protein sequence which will be filled into the pbc box
# - BOX_X, BOX_Y, BOX_Z: The dimensions of the box for the SAPGenPBC.py script
# - DENSITY: The target density to which SAPGenPBC.py fills the box with chains
#
# Output: all.pdb - the final concatenated PDB file
#

tools=~/tools_ua_gecko/

sequence=$1
box_x=$2
box_y=$3
box_z=$4
density=$5

# Start 48 generator jobs in parallel
for i in $(seq 1 1 48)
do
    mkdir uniform_"$i"
    cd uniform_"$i"

    # Run the generator script in the background and store its process ID
    nohup python $tools/SAPGenPBC.py $sequence $box_x $box_y $box_z $density &
    pid[$i]=$!
    cd ..
done

# Check continuously if any generator has finished
finished=0
generator_id=0
while [ $finished -eq 0 ]
do
    for i in $(seq 1 1 48)
    do
        # If the generator job with process ID ${pid[$i]} is not running, it has finished
        if ! kill -0 ${pid[$i]} 2>/dev/null; then
            finished=1
            generator_id=$i
            break
        fi
        sleep 60
    done
done

# Stop other generators and delete their directories
for j in $(seq 1 1 48)
do
    if [ $j -ne $generator_id ]; then
    echo "kill $j"
        # Stop the generator job with the process ID ${pid[$j]}
        kill ${pid[$j]} 2>/dev/null
        wait ${pid[$j]} 2>/dev/null
        # Remove the directory for the stopped generator job
        rm -r uniform_"$j"
    fi
done

# Concatenate PDB files in the completed generator directory
echo $generator_id
cd uniform_"$generator_id"
N=$(ls Chain*.pdb | wc -l)
# Run the concatenate_pdb.sh script to concatenate the PDB files
$tools/concatenate_pdbs.sh $N $box_x $box_y $box_z
cp all.pdb ../all.pdb
cp Chain* ../
cd ..
rm -r uniform_"$generator_id"
