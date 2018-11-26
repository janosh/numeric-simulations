/*
* This program carries out a gravitational force computation in a 2D periodic space using the particle-mesh (PM) method
* with clouds-in-cell (CIC) charge/mass assignment. We wish to measure the effective force law the method delivers.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>


// this function inserts a particle of unit mass into the grid according to the CIC mass assignment
void insert_particle(fftw_complex *grid, int N_grid, double *pos) {

    srand48(clock());  // initializes the random number generation

    for (int i = 0; i < 2; i++) {
        pos[i] = drand48();   // generates random doubles in the range [0,1]
    }

    double cell_f[2];  // floating-point cell index
    for (int i = 0; i < 2; i++) {
        cell_f[i] = pos[i] * N_grid - 0.5;
    }

    int cell[2];    // cell index
    for (int i = 0; i < 2; i++) {
        cell[i] = floor(cell_f[i]);
    }

    grid[cell[0] * N_grid + cell[1]][0] = (1 - (cell_f[0] - cell[0])) * (1 - (cell_f[1] - cell[1]));
    grid[N_grid * ((cell[0] + 1) % N_grid) + cell[1]][0] = (cell_f[0] - cell[0]) * (1 - (cell_f[1] - cell[1]));
    grid[N_grid * cell[0] + (cell[1] + 1) % N_grid][0] = (1 - (cell_f[0] - cell[0])) * (cell_f[1] - cell[1]);
    grid[N_grid * ((cell[0] + 1) % N_grid) + (cell[1] + 1) % N_grid][0] = (cell_f[0] - cell[0]) * (cell_f[1] - cell[1]);
}


void get_potential(fftw_complex *grid, int N_grid, fftw_complex *potential) {
    // allocate memory to hold the Fourier transformed grid
    fftw_complex *grid_kspace = calloc(N_grid * N_grid, sizeof(fftw_complex));

    // setup up Fourier transformation by passing FFTW the necesasary array pointers and those array's dimensions
    fftw_plan plan_density = fftw_plan_dft_2d(N_grid, N_grid, grid, grid_kspace, FFTW_FORWARD, FFTW_ESTIMATE);

    // now do the actual transform
    fftw_execute(plan_density);

    fftw_complex *potential_kspace = calloc(N_grid * N_grid, sizeof(fftw_complex));

    potential_kspace[0][0] = 0; // here we treat the special case of k = 0
    potential_kspace[0][1] = 0; // by setting the potential's real and imaginary part to zero

    for (int i = 1; i < N_grid * N_grid; i++) {
        int l1 = (int) i/N_grid, l2 = i % N_grid;

        if (2 * l1 < N_grid && 2 * l2 < N_grid) {
            potential_kspace[i][0] = -1/(M_PI * (l1 * l1 + l2 * l2)) * grid_kspace[i][0];
            potential_kspace[i][1] = -1/(M_PI * (l1 * l1 + l2 * l2)) * grid_kspace[i][1];
        }

        if (2 * l1 >= N_grid && 2 * l2 < N_grid) {
            potential_kspace[i][0] = -1/(M_PI * ((l1 - N_grid) * (l1 - N_grid) + l2 * l2)) * grid_kspace[i][0];
            potential_kspace[i][1] = -1/(M_PI * ((l1 - N_grid) * (l1 - N_grid) + l2 * l2)) * grid_kspace[i][1];
        }

        if (2 * l1 < N_grid && 2 * l2 >= N_grid) {
            potential_kspace[i][0] = -1/(M_PI * (l1 * l1 + (l2 - N_grid) * (l2 - N_grid))) * grid_kspace[i][0];
            potential_kspace[i][1] = -1/(M_PI * (l1 * l1 + (l2 - N_grid) * (l2 - N_grid))) * grid_kspace[i][1];
        }

        if (2 * l1 >= N_grid && 2 * l2 >= N_grid) {
            potential_kspace[i][0] = -1/(M_PI * ((l1 - N_grid) * (l1 - N_grid) + (l2 - N_grid) * (l2 - N_grid))) * grid_kspace[i][0];
            potential_kspace[i][1] = -1/(M_PI * ((l1 - N_grid) * (l1 - N_grid) + (l2 - N_grid) * (l2 - N_grid))) * grid_kspace[i][1];
        }
    }

    fftw_plan plan_potential = fftw_plan_dft_2d(N_grid, N_grid, potential_kspace, potential, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(plan_potential);

}


void get_forces(fftw_complex *potential, int N_grid, double force_field[N_grid][N_grid][2]) {
    for (int i = 0; i < N_grid; i++) {
        for (int j = 0; j < N_grid; j++) {
            force_field[i][j][0] = -(potential[N_grid * ((i + 1) % N_grid) + j][0] - potential[((i == 0) ? (N_grid - 1) : (i - 1)) * N_grid + j][0]) * N_grid/2;
            force_field[i][j][1] = -(potential[i * N_grid + (j + 1) % N_grid][0] - potential[i * N_grid + ((j == 0) ? (N_grid - 1) : (j - 1))][0]) * N_grid/2;
        }
    }
}


void get_100_samples(double *pos, int N_grid, double force_field[N_grid][N_grid][2], FILE *radii, FILE *forces) {

    srand48(clock());  // initializes the random number generation

    double r_min = 0.3/N_grid, r_max = 0.5;
    double forces_at_points[100][2] = {0}, point_grid[N_grid][N_grid];

    for (int i = 0; i < 100; i++) {
        // create 100 random points around the particle created earlier using the following equations and store them in random_points
        double p = drand48(), q = drand48();
        double deltax = r_min * pow(r_max/r_min, p) * cos(2 * M_PI * q);
        double deltay = r_min * pow(r_max/r_min, p) * sin(2 * M_PI * q);
        double rad = sqrt(deltax * deltax + deltay * deltay);

        fwrite(&rad, sizeof(double), 1, radii);

        double random_point[2] = {fmod(pos[0] + deltax, 1), fmod(pos[1] + deltay, 1)};  // fmod(x, y) is the modulus function for floats (defined in math.h)

        // now place these points into the force grid according to the CIC method
        double cell_f[2];  // floating-point cell index
        for (int j = 0; j < 2; j++) {
            cell_f[j] = random_point[j] * N_grid - 0.5;
        }

        int cell[2];    // cell index
        for (int j = 0; j < 2; j++) {
            cell[j] = floor(cell_f[j]);
        }

        point_grid[cell[0]][cell[1]] = (1 - (cell_f[0] - cell[0])) * (1 - (cell_f[1] - cell[1]));
        point_grid[(cell[0] + 1) % N_grid][cell[1]] = (cell_f[0] - cell[0]) * (1 - (cell_f[1] - cell[1]));
        point_grid[cell[0]][(cell[1] + 1) % N_grid] = (1 - (cell_f[0] - cell[0])) * (cell_f[1] - cell[1]);
        point_grid[(cell[0] + 1) % N_grid][(cell[1] + 1) % N_grid] = (cell_f[0] - cell[0]) * (cell_f[1] - cell[1]);

        for (int j = 0; j < 2; j++) {
            forces_at_points[i][j] = point_grid[cell[0]][cell[1]] * force_field[cell[0]][cell[1]][j] + point_grid[(cell[0] + 1) % N_grid][cell[1]] * force_field[(cell[0] + 1) % N_grid][cell[1]][j] + point_grid[cell[0]][(cell[1] + 1) % N_grid] * force_field[cell[0]][(cell[1] + 1) % N_grid][j] + point_grid[(cell[0] + 1) % N_grid][(cell[1] + 1) % N_grid] * force_field[(cell[0] + 1) % N_grid][(cell[1] + 1) % N_grid][j];
        }

        double force_magn = sqrt(forces_at_points[i][0] * forces_at_points[i][0] + forces_at_points[i][1] * forces_at_points[i][1]);
        fwrite(&force_magn, sizeof(double), 1, forces);


        for (int j = 0; j < N_grid; j++) {
            for (int k = 0; k < N_grid; k++) {
                point_grid[j][k] = 0;
            }
        }
    }

}


int main(void) {

    const int N_grid = 256;    // number of grid cells

    // data file pointers to store the radii and force values
    FILE *radii = fopen("radii.dat", "w"), *forces = fopen("forces.dat", "w");

    for (int j = 0; j < 10; j++) {

        // create and initialize the grid array
        fftw_complex *grid = calloc(N_grid * N_grid, sizeof(fftw_complex));

        double pos[2];  // the particle's position in the 2D periodic space of unit length

        // insert a particle at a random position into the grid
        insert_particle(grid, N_grid, pos);

        // create and initialize the potential array
        fftw_complex *potential = calloc(N_grid * N_grid, sizeof(fftw_complex));
        for(int row = 0; row < N_grid; row++)
        for(int col = 0; col < N_grid; col++)
        {
            potential[row * N_grid + col][0] = 0;  // real part
            potential[row * N_grid + col][1] = 0;  // imaginary part
        }

        // calculate the potential based on the current mass distribution in the grid array
        get_potential(grid, N_grid, potential);

        // create and initialize the force array
        double force_field[N_grid][N_grid][2] = {0};

        // calculate and store the forces in the force array
        get_forces(potential, N_grid, force_field);

        // execute the function to generate
        get_100_samples(pos, N_grid, force_field, radii, forces);
    }
}
