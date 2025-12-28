// All the codes here will be write in C
#include <stdio.h> // contains definition of useful fonctions (fprintf, fscanf, fopen ...)
#include <stdlib.h>// required for malloc
#include <math.h> // for math functions
#include <string.h>


double** malloc_2d(size_t m, size_t n);
void free_2d(double** a);


// +--------------------------------------------------------------------
// |  Allocation functions
// +--------------------------------------------------------------------


// Atom coordinates will be provided in nanometers.
// You will need to allocate two-dimensional arrays. 
// The following functions should help you allocate and free 2-dimensional arrays:
double** malloc_2d(size_t m, size_t n)
{
    double** a = malloc(m*sizeof(double*)); // memory allocation in 1 dimension 
    if (a == NULL)
    {
        return NULL;
    }
    
    a[0] = malloc(n*m*sizeof(double));
    
    if (a[0] == NULL)
    {
        free(a);
        return NULL;
    }
    
    for (size_t i=1 ; i<m ; i++)
    {
        a[i] = a[i-1]+n;
    }

    return a;
}


void free_2d(double** a)
{
    free(a[0]);
    a[0] = NULL;
    free(a);
}


// +--------------------------------------------------------------------
// |  Functions headers
// +--------------------------------------------------------------------

// Throughout this tutorial, we will use Argon atoms (mass = 39.948 g/mol) with the following parameters:

// ϵ : 0.0661 j/mol
static const double EPSILON = 0.0661; // j/mol

// σ : 0.3345 nm
static const double SIGMA = 0.3345; // nm

size_t read_Natoms(FILE* input_file);

void read_molecule( FILE* input_file,
                    size_t Natoms,
                    double** coord,
                    double* mass);

void compute_distances( size_t Natoms,
                        double** coord,
                        double** distance);

double VLJ(double epsilon,
           double sigma,
           double r);

double V(double epsilon,
         double sigma,
         size_t Natoms,
         double** distance);

void init_to_zero_velocity(size_t Natoms,
                           double** velocity);

double T(size_t Natoms,
         double** velocity,
         double* mass);

double E(double epsilon,
         double sigma,
         size_t Natoms,
         double** distance,
         double** velocity,
         double* mass);

double U(double epsilon,
         double sigma,
         double r);

void compute_acc(size_t Natoms,
                 double** coord,
                 double* mass,
                 double** distance,
                 double** acceleration);

void Verlet_algorithm(double time_step,
                 size_t Natoms,
                 double** coord,
                 double* mass,
                 double** distance,
                 double** velocity,
                 double** acceleration);

void compute_trajectory(double time_step,
                        size_t max_steps,
                        FILE* input_file,
                        FILE* output_file);


// +--------------------------------------------------------------------
// |  Functions implementation
// +--------------------------------------------------------------------


// 1 Describing the Atoms
// We will create functions to read the atomic data from an input file. 
//
// The input file follows this format:
// The first line contains the number of atoms (Natoms).
// Each subsequent line contains the x, y, and z coordinates followed by the mass of an atom.



 
// Create a function with the following prototype which reads the number of atoms from an opened input file:
size_t read_Natoms(FILE* input_file)
{
    size_t Natoms = 0;
    // The first line contains the number of atoms (Natoms).
    // size_t n = read(fd , &number_of_atoms , 1 * sizeof(int));
    //  size_t n = fread(&Natoms, sizeof(int), 1, input_file);
    fscanf(input_file, "%zu", &Natoms);
    return Natoms;
}


// Then, write a function with the following prototype:
void read_molecule( FILE* input_file, 
                    size_t Natoms,    
                    double** coord,   
                    double* mass)     
{
    // The first line contains the number of atoms (Natoms).
    // Each subsequent line contains the x, y, and z coordinates followed by the mass of an atom.
    for (size_t i = 0; i < Natoms; i++)
    {   
        // fread(&(mass[i]), sizeof(double), 1, input_file);
        fscanf(input_file, "%lf %lf %lf %lf", &coord[i][0], &coord[i][1], &coord[i][2], &mass[i]);
    }
}
// This function reads the atomic coordinates and masses, and stores the data in the arrays coord and mass provided as parameters. 
// coord is a two dimensional array allocated such that coord[i][2] returns the z coordinate of atom i.

// Write another function which takes as input the array of coordinates
// and returns internuclear distances between each pair in a two-dimensional
// double array of size (Natoms × Natoms)
void compute_distances(size_t Natoms,
            double** coord,
            double** distance)
{

    for (size_t i = 0; i < Natoms; i++)
        for (size_t j = 0; j < Natoms; j++)
        {
            double dx = coord[i][0] - coord[j][0];
            double dy = coord[i][1] - coord[j][1];
            double dz = coord[i][2] - coord[j][2];
            distance[i][j] = sqrt(dx * dx + dy * dy + dz * dz); 

        }
}

// To use the sqrt function in C, you need to include <math.h>, 
// and you also need to compile using the -lm option to link with the libm.so math library.
// You might need to put the -lm option at the end of the command line (after your files).


// 2 The Lennard-Jones potential

// Write a function which takes as input ϵ, σ, the number of atoms and the array of distances, 
// and computes the total potential energy where VLJ(r) is the Lennard-Jones potential VLJ(r): 
double VLJ(double epsilon,
           double sigma,
           double r) 
{
    return 4*epsilon*(pow(sigma/r,12)-pow(sigma/r,6));
}
// To compute x^n in C, you will need to call the power function located in the libm.so library, defined in math.h.
 
// We will now write a function to calculate the total potential energy of the system using the Lennard-Jones potential. 
// This function will take the following inputs: ϵ, σ, the number of atoms (Natoms), 
// and the array containing the distances between each pair of atoms.
double V(double epsilon,
         double sigma,
         size_t Natoms,
         double** distance)
{
    double LJ = 0;

    for (size_t i = 0; i < Natoms; i++)
        for (size_t j = i+1; j < Natoms; j++)
        {
            LJ += VLJ(epsilon, sigma, distance[i][j]);                
        }

    return LJ;
}


// 3 Computing the total energy

// We will now write functions to compute the total kinetic energy and the total energy of the system.
// First, we will write a function to calculate the total kinetic energy of the system T. 
// This function will take as input the number of atoms, the array of velocities (3×Natoms) and masses (Natoms). 
// The velocities are initialized to zero. 
void init_to_zero_velocity(size_t Natoms,
                           double** velocity)
{
    for (size_t i = 0; i < Natoms; i++)
    {
        velocity[i][0] = 0.0;
        velocity[i][1] = 0.0;
        velocity[i][2] = 0.0;
    }
}

// The total kinetic energy T is given by T.
double T(size_t Natoms,
         double** velocity,
         double* mass)
{
    double K = 0;

    for (size_t i = 0; i < Natoms; i++)
    {
    // velocity norm is sqrt(vix² + viy² + viz²)
    // velocity norm squared is vix² + viy² + viz²
        double velocity_squared_i = velocity[i][0] * velocity[i][0] + velocity[i][1] * velocity[i][1] + velocity[i][2] * velocity[i][2];
        K += mass[i]*velocity_squared_i;
    }

    return (1.0/2.0)*K;
}

// Next, we will write a function to compute the total energy of the system,
// which is the sum of the total kinetic energy T and the total potential energy V. 
// The total energy E is given by: E = T + V.

double E(double epsilon,
         double sigma,
         size_t Natoms,
         double** distance,
         double** velocity,
         double* mass)
{
    return T(Natoms, velocity, mass) +  V(epsilon, sigma, Natoms, distance);
}


// 4 Computing the acceleration

// We need to write a function that computes the acceleration vector for each atom and stores it in a double precision array.
double U(double epsilon,
         double sigma,
         double r)
{
    return 24*(epsilon/r)*(pow(sigma/r,6)-2*pow(sigma/r,12));
}

void compute_acc(size_t Natoms,
                 double** coord,
                 double* mass,
                 double** distance,
                 double** acceleration)
{
        for (size_t i = 0; i < Natoms; i++)
            for (size_t k = 0; k < 3; k++)
			{
				double sum = 0;
				for (size_t j = 0; j < Natoms; j++)
				{
                    if (i!=j)
                    {
                        // fprintf(stdout, "distance[%zu][%zu]=%lf\n", i, j, distance[i][j]); // nan = Not a number if distance == 0
                        sum += U(EPSILON, SIGMA, distance[i][j]) * (coord[i][k] - coord[j][k]) / distance[i][j];
                    }
				}
				acceleration[i][k] = -(1.0/mass[i]) * sum;
			}
}

void Verlet_algorithm(double time_step, size_t Natoms, double** coord, double* mass, double** distance, double** velocity, double** acceleration)
{
    size_t iter_atom = 0;
    size_t iter_coordinate = 0;

    // 1. Update position vector for each atom
    for (iter_atom=0; iter_atom<Natoms; ++iter_atom)
    {
        for (iter_coordinate=0; iter_coordinate<3; ++iter_coordinate)
        {
            coord[iter_atom][iter_coordinate] = coord[iter_atom][iter_coordinate] + velocity[iter_atom][iter_coordinate]*time_step
                                                + acceleration[iter_atom][iter_coordinate]*(time_step*time_step)/2.0;
        }
    }

    // 2. Update the part of the velocity depending on previous acceleration
    for (iter_atom=0; iter_atom<Natoms; ++iter_atom)
    {
        for (iter_coordinate=0; iter_coordinate<3; ++iter_coordinate)
        {
            velocity[iter_atom][iter_coordinate] = velocity[iter_atom][iter_coordinate] + acceleration[iter_atom][iter_coordinate]*time_step/2.0;
        }
    }

    // 3. Update acceleration
    // 3.1 Start by computing new distances
    compute_distances(Natoms, coord, distance);

    // 3.2 Compute acclerations
    compute_acc(Natoms, coord, mass, distance, acceleration);

    // 4. Update the part of the velocity depending on new acceleration
    for (iter_atom=0; iter_atom<Natoms; ++iter_atom)
    {
        for (iter_coordinate=0; iter_coordinate<3; ++iter_coordinate)
        {
            velocity[iter_atom][iter_coordinate] = velocity[iter_atom][iter_coordinate] + acceleration[iter_atom][iter_coordinate]*time_step/2.0;
        }
    }
}

void compute_trajectory(double time_step, size_t max_steps, FILE* input_file, FILE* output_file)
{
    size_t step = 0;
    size_t iter_atom = 0;

    // 1. Describing the Atoms
    size_t Natoms = read_Natoms(input_file);

    // Allocation
    double** coord = malloc_2d(Natoms, 3);
    double*  mass  = malloc(Natoms * sizeof(double));
    double** distance = malloc_2d(Natoms, Natoms);
    double** velocity = malloc_2d(Natoms, 3);
    double** acceleration = malloc_2d(Natoms, 3);
        
    // Read input file
    read_molecule(input_file, Natoms, coord, mass);

    // Fill distance
    compute_distances(Natoms, coord, distance);

    // 2. The Lennard-Jones potential
    double Vr = V(EPSILON, SIGMA, Natoms, distance);

    // 3. Computing the total energy
    // Init all velocities to zero
    init_to_zero_velocity(Natoms, velocity);

    // Compute kinetic energy
    double Tr = T(Natoms, velocity, mass);

    double Total_Energy = E(EPSILON, SIGMA, Natoms, distance, velocity, mass);
    
    // 4. Computing the acceleration
    compute_acc(Natoms, coord, mass, distance, acceleration);

    // Print initial state
    if (output_file != NULL)
    {
        fprintf(output_file, "%zu\n", Natoms);
        fprintf(output_file, "# T:%lf V:%lf E:%lf \n", Tr, Vr, Total_Energy);
        for (iter_atom=0; iter_atom < Natoms; ++iter_atom)
        {
           fprintf(output_file, "Ar %lf %lf %lf \n", coord[iter_atom][0], coord[iter_atom][1], coord[iter_atom][2]); 
        //    fprintf(output_file, "Ar vx:%lf vy:%lf vz:%lf \n", velocity[iter_atom][0], velocity[iter_atom][1], velocity[iter_atom][2]); 
        //    fprintf(output_file, "Ar ax:%lf ay:%lf az:%lf \n", acceleration[iter_atom][0], acceleration[iter_atom][1], acceleration[iter_atom][2]); 
        }
    }

    // Compute each step
    for (step=0; step<max_steps; ++step)
    {
        // Propagate coordinates and velocities
        Verlet_algorithm(time_step, Natoms, coord, mass, distance, velocity, acceleration);

        // Print new state
        if (output_file != NULL)
        {
            // Compute Energies
            Vr = V(EPSILON, SIGMA, Natoms, distance);
            Tr = T(Natoms, velocity, mass);
            Total_Energy = E(EPSILON, SIGMA, Natoms, distance, velocity, mass);

            // Print trajectory
            fprintf(output_file, "%zu\n", Natoms);
            fprintf(output_file, "# T:%lf V:%lf E:%lf \n", Tr, Vr, Total_Energy);
            for (iter_atom=0; iter_atom < Natoms; ++iter_atom)
            {
                fprintf(output_file, "Ar %lf %lf %lf \n", coord[iter_atom][0], coord[iter_atom][1], coord[iter_atom][2]); 
                // fprintf(output_file, "Ar vx:%lf vy:%lf vz:%lf \n", velocity[iter_atom][0], velocity[iter_atom][1], velocity[iter_atom][2]); 
                // fprintf(output_file, "Ar ax:%lf ay:%lf az:%lf \n", acceleration[iter_atom][0], acceleration[iter_atom][1], acceleration[iter_atom][2]); 
            }
        }

    }
    
    // Free memory
    free_2d(acceleration);
    free_2d(velocity);
    free_2d(distance);
    free(mass);
    free_2d(coord);

}


// +--------------------------------------------------------------------
// |  Test main
// +--------------------------------------------------------------------


// Try the functions
int main(int argc, char *argv[])
{
    const char * fichier_molecules = argv[1];
    double time_step = atof(argv[2]);
    size_t max_steps = atol(argv[3]);

    const char * trajectory_file = "./trajectory.txt";
    
    // Open file
    FILE * input_file = fopen(fichier_molecules, "r");
    FILE * output_file = fopen(trajectory_file, "w");
    
    if ((input_file == NULL) || (output_file == NULL))
    {
        return -1;
    }
    else
    {
        compute_trajectory(time_step, max_steps, input_file, output_file);

        // Close file
        fclose(input_file);
        fclose(output_file);
    }        
    
    
    return 0;
}





