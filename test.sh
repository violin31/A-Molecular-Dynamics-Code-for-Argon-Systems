#!/bin/bash

# This script is testing that the main program does not return an error when running. In the absence of references, the outputs can't be tested.
# Create directory lib in case it doesn't exist.
mkdir -p lib

# Checking that the code compiles
gcc -o lib/test_main src/dynamics.c -lm

# Check if compilation was successful
if [ $? -ne 0 ]; then
    echo "Compilation failed."
    exit 1
fi

# Create a test case that will be used as input
echo "Preparing test cases..."
echo "5" > test_input.txt
echo "0.0    0.0    0.0    39.948" >> test_input.txt
echo "0.0    0.0    0.5    39.948" >> test_input.txt
echo "0.1    0.2   -0.5    39.948" >> test_input.txt
echo "0.3    0.0   -0.8    39.948" >> test_input.txt
echo "-0.4    0.0   0.2    39.948" >> test_input.txt

# Run the program, using the compile_and_run.sh file
echo "Running the program..."
./compile_and_run.sh ./test_input.txt 0.01 1000   # three arguments to run: `compile_and_run.sh <file_path> <delta_time> <max_steps>` where `file_path` is the input file containing the initial state of the molecule, `delta_time` is the time step used for the Verlet Algorithm and `max_steps` is the maximum number of iterations of the Verlet Algorithm.

# Check the output
if [ $? -ne 0 ]; then
    echo "Program execution failed. Check your code."
    exit 1
fi

echo "Test passed, the program ran without errors. Warning: this does not guarantee that the outputs are right. Output can be found in the \"trajectory.txt\" file."
