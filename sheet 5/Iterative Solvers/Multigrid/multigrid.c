#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


// modulo function to implement periodic boundaries on the grid
int mod(int x, int N) {
    x = x % N;
    return (x < 0) ? (x + N) : x;
}


// This function improves the solution Ax=b with one Gauss-Seidel step. The result is used to overwrite x. The matrix A is hardcoded and represents the Poisson equation.
void gaussseidel_step(double *phi, double *rho, int N, double h, double G)
{
    for(int i = 0; i < N; i++)
    for(int j = 0; j < N; j++)
    {
        phi[i * N + j] = (phi[mod(i + 1, N) * N + j] + phi[i * N + mod(j + 1, N)] + phi[mod(i - 1, N) * N + j] + phi[i * N + mod(j - 1, N)] - 4 * M_PI * G * h * h * rho[i * N + mod(j - 1, N)])/4;
    }
}


void calc_residual(double *phi, double *rho, double *res, int N, double h, double G)
{
    for(int i = 0; i < N; i++)
    for(int j = 0; j < N; j++)
    {
        res[i * N + j] = 4 * M_PI * G * h * h * rho[i * N + j] - (phi[mod(i + 1, N) * N + j] + phi[mod(i - 1, N) * N + j] + phi[i * N + mod(j + 1, N)] + phi[i * N + mod(j - 1, N)] - 4 * phi[i * N + j]);
    }
}


// This function calculates the norm of the vector of length N, defined as the usual quadratic vector norm.
double norm_of_residual(double *res, int N)
{
    double res_sum = 0;

    for(int i = 0; i < N; i++)
    for(int j = 0; j < N; j++)
    res_sum += res[i * N + j] * res[i * N + j];

    return sqrt(res_sum);
}


// This function restricts the NxN mesh stored in 'fine[ ]' to the N/2xN/2 mesh stored in 'coarse[ ]'.
void do_restrict(int N, double *fine, double *coarse)
{
    for(int i = 0; i < N; i += 2)
    for(int j = 0; j < N; j += 2)
    {
        coarse[i/2 * N/2 + j/2] = fine[mod(i-1, N) * N + mod(j-1, N)]/16 + fine[mod(i-1, N) * N + j]/8 + fine[mod(i-1, N) * N + mod(j+1, N)]/16 + fine[i * N + mod(j-1, N)]/8 + fine[i * N + j]/4 + fine[i * N + mod(j+1, N)]/8 + fine[mod(i+1, N) * N + mod(j-1, N)]/16 + fine[mod(i+1, N) * N + j]/8 + fine[mod(i+1, N) * N + mod(j+1, N)]/16;
    }
}


// This function interpolates the the NxN mesh stored in 'coarse[ ]' to the 2Nx2N mesh stored in 'fine[ ]'
void do_prolong(int N, double *coarse, double *fine)
{
    for(int i = 0; i < 2 * N; i += 2)
    for(int j = 0; j < 2 * N; j += 2)
    {
        fine[i * 2 * N + j] = coarse[i/2 * N + j/2];
    }
    for(int i = 0; i < 2 * N; i += 2)
    for(int j = 1; j < 2 * N; j += 2)
    {
        fine[i * 2 * N + j] = coarse[i/2 * N + j/2]/2 + coarse[i/2 * N + mod((j+1)/2, N)]/2;  // intended integer divison with remainder j/2
    }
    for(int i = 1; i < 2 * N; i += 2)
    for(int j = 0; j < 2 * N; j += 2)
    {
        fine[i * 2 * N + j] = coarse[i/2 * N + j/2]/2 + coarse[mod((i+1)/2, N) * N + j/2]/2;
    }
    for(int i = 1; i < 2 * N; i += 2)
    for(int j = 1; j < 2 * N; j += 2)
    {
        fine[i * 2 * N + j] = (coarse[i/2 * N + j/2] + coarse[mod((i+1)/2, N) * N + j/2] + coarse[i/2 * N + mod((j+1)/2, N)] + coarse[mod((i+1)/2, N) * N + mod((j+1)/2, N)])/4;
    }
}


// This function carries out a V-cycle using the multigrid approach, down to N = 4.
void do_v_cycle(double *phi, double *rho, int N, double h, double G)
{
    gaussseidel_step(phi, rho, N, h, G);

    if(N > 4)
    {
        // allocate some storage for the residual and error on the current and a coarsened mesh
        double *res = calloc(N * N, sizeof(double));
        double *err = calloc(N * N, sizeof(double));
        double *coarse_res = calloc(N/2 * N/2, sizeof(double));
        double *coarse_err = calloc(N/2 * N/2, sizeof(double));

        // calculate the residual
        calc_residual(phi, rho, res, N, h, G);

        // restrict the residual
        do_restrict(N, res, coarse_res);

        // now multiply the residual on the coarser mesh (our b) with a factor of 4 to take care of the (h) -> (2h) change, and the fact that we do not change A
        for(int i = 0; i < N/2; i++)
        for(int j = 0; j < N/2; j++)
        coarse_res[i * N/2 + j] *= 4;

        // call the V-cycle again, recursively
        do_v_cycle(coarse_err, coarse_res, N/2, h, G);

        // now interpolate the error to the finer mesh
        do_prolong(N/2, coarse_err, err);

        // finally add the error to the current solution
        for(int i = 0; i < N; i++)
        for(int j = 0; j < N; j++)
        {
            phi[i * N + j] += err[i * N + j];
        }

        // free up temporary storage
        free(coarse_err); free(coarse_res); free(err); free(res);
    }

    gaussseidel_step(phi, rho, N, h, G);
}


int main(int argc, char **argv)
{
    const int N = 256, Nsteps = 2000;
    const double L = 1.0, h = L / N, G = 1.0, eta = L/10, rho0 = 10;


    // allocate some storage for our fields
    double *phi = calloc(N * N, sizeof(double));
    double *rho = calloc(N * N, sizeof(double));
    double *res = calloc(N * N, sizeof(double));


    // now set-up the density field
    for(int i = 0; i < N; i++)
    for(int j = 0; j < N; j++)
    {
        double dx = (i - N / 2 + 0.5) * h, dy = (j - N / 2 + 0.5) * h;
        double r2 = dx * dx + dy * dy;

        rho[i * N + j] = rho0 * exp(-r2 / (2 * eta * eta));
    }

    // open a file for outputting the residual values, and then do 2000 Jacobi steps
    FILE *datafile = fopen("res_multigrid.txt", "w");
    printf("\nMultigrid\n");
    for(int i = 0; i <= Nsteps; i++)
    {
        do_v_cycle(phi, rho, N, h, G);

        calc_residual(phi, rho, res, N, h, G);

        double S = norm_of_residual(res, N);

        if (10 * i % Nsteps == 0) {
            printf("step: %4d\t\tresidual = %g\n", i, S);
        }
        fprintf(datafile, "%g\n", S);
    }

    double phi_sum = 0;
    FILE *phi_data = fopen("phi_mg.dat", "w");
    for(int i = 0; i < N; i++)
    for(int j = 0; j < N; j++)
    {
        phi_sum += phi[i * N + j];
        fwrite(&phi[i * N + j], sizeof(double), 1, phi_data);
    }
    printf("phi_sum = %g\n", phi_sum);
}
