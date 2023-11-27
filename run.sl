#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:01:00
#SBATCH --output=test.out
#SBATCH -A anakano_429

module purge
module load gcc/11.3.0

# Update the path to your executable and input file
EXECUTABLE="./raytracer"
INPUT_FILE="./InputFiles/spheres.scene"
OUTPUT_FILE="test.png"

# Run the executable
$EXECUTABLE $INPUT_FILE $OUTPUT_FILE

scp $OUTPUT_FILE aamarath@discovery.usc.edu:home1/aamarath/CSCI596_Final_Project
