#!/bin/bash
INPUT_SCENE="./InputFiles/spheres.scene"
OUTPUT_FILE="localRun.out"

./raytracer 1   $INPUT_SCENE localRun.png > $OUTPUT_FILE
./raytracer 2   $INPUT_SCENE localRun.png >> $OUTPUT_FILE
./raytracer 4   $INPUT_SCENE localRun.png >> $OUTPUT_FILE
./raytracer 8   $INPUT_SCENE localRun.png >> $OUTPUT_FILE
./raytracer 16  $INPUT_SCENE localRun.png >> $OUTPUT_FILE
./raytracer 32  $INPUT_SCENE localRun.png >> $OUTPUT_FILE
./raytracer 64  $INPUT_SCENE localRun.png >> $OUTPUT_FILE