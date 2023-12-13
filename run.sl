#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=00:20:00
#SBATCH --output=carcRun.out
#SBATCH -A anakano_429

module purge
module load gcc/11.3.0


# Run the executable
INPUT_SCENE="./InputFiles/spheres.scene"
OUTPUT_FILE="carcRun.out"

./raytracer 1   $INPUT_SCENE localRun.png > $OUTPUT_FILE
./raytracer 2   $INPUT_SCENE localRun.png >> $OUTPUT_FILE
./raytracer 4   $INPUT_SCENE localRun.png >> $OUTPUT_FILE
./raytracer 8   $INPUT_SCENE localRun.png >> $OUTPUT_FILE
./raytracer 16  $INPUT_SCENE localRun.png >> $OUTPUT_FILE
./raytracer 32  $INPUT_SCENE localRun.png >> $OUTPUT_FILE
./raytracer 64  $INPUT_SCENE localRun.png >> $OUTPUT_FILE