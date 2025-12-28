# Architecture of the repository

In `lib` you will find the executable.  This directory is created when `compile_and_run.sh` is run.  
In `src` you will find the source files for the project.

# Sources

The main source is `dynamics.c`. It contains the main function and private functions to compute the trajectory of the atoms in a molecule with Verlet algorithm. A documentation for those functions is available at the beginning of the file.

# Output

The output of the algorithm is a file `trajectory.txt` stored in the same directory as the `compile_and_run.sh` file.

# Testing

In order to test the code, a `test.sh` file is provided, which runs the code with a given output and verifies that it did not return an error.
