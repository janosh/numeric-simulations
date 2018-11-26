#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


/* this is our data type for storing particle properties */
typedef struct {
    double pos[3], vel[3], acc[3], pot;
} part;


/* This function initializes our particle set. */
void initialize_system(part *particle, int N1d, double mean_spacing, double sigma)
{
    srand48(42);

    int n = 0;
    for(int i = 0; i < N1d; i++)
    for(int j = 0; j < N1d; j++)
    for(int k = 0; k < N1d; k++, n++)
    {
        particle[n].pos[0] = mean_spacing * i;
        particle[n].pos[1] = mean_spacing * j;
        particle[n].pos[2] = mean_spacing * k;

        for(int l = 0; l < 3; l++)
        {
            double u1, u2;
            do {
                u1 = drand48(), u2 = drand48();
            } while(u1 == 0);   /* u1 may not take the value 0 to ensure that the logarithm is well-defined */

            particle[n].vel[l] = sigma * sqrt(-2.0 * log(u1)) * cos(2 * M_PI * u2);
        }
    }
}


/* This function updates the velocities by applying the accelerations for the given time interval. */
void kick(part *particle, int Ntot, double dt)
{
    for(int i = 0; i < Ntot; i++)
    for(int k = 0; k < 3; k++)
    particle[i].vel[k] += particle[i].acc[k] * dt;
}


/* This function drifts the particles with their velocities for the given time interval. The particles are mapped periodically back to the box if needed. */
void drift(part *particle, int Ntot, double boxsize, double dt)
{
    for(int i = 0; i < Ntot; i++)
    for(int k = 0; k < 3; k++)
    {
        particle[i].pos[k] += particle[i].vel[k] * dt;

        while(particle[i].pos[k] >= boxsize)
        particle[i].pos[k] -= boxsize;
        while(particle[i].pos[k] < 0)
        particle[i].pos[k] += boxsize;
    }
}


/* This function calculates the potentials and forces for all particles. For simplicity, we do this by going through all particle pairs i-j, and then adding the contributions both to i and j. */
void calc_forces(part *particle, int Ntot, double boxsize, double rcut)
{
    /* first, set all the accelerations and potentials to zero */
    for(int i = 0; i < Ntot; i++)
    {
        particle[i].pot = 0;
        for(int k = 0; k < 3; k++)
        particle[i].acc[k] = 0;
    }

    /* sum over all distinct pairs */
    for(int i = 0; i < Ntot; i++)
    for(int j = i + 1; j < Ntot; j++)
    {
        double dist_sqr = 0, dr[3];
        /* we calculate the Euclidean distance between two particles */
        for(int k = 0; k < 3; k++)
        {
            dr[k] = particle[i].pos[k] - particle[j].pos[k];

            /* this is to heed boundary conditions which render particles closer to the correspoding particle in the previous/next box if their separation is larger than half/smaller than minus half the box size */
            if(dr[k] > 0.5 * boxsize)
            dr[k] -= boxsize;
            if(dr[k] < -0.5 * boxsize)
            dr[k] += boxsize;

            dist_sqr += dr[k] * dr[k];
        }

        double r6 = dist_sqr * dist_sqr * dist_sqr, r12 = r6 * r6;
        if (dist_sqr < rcut * rcut) {
            double pot = 4 * (1/r12 - 1/r6);
            particle[i].pot += pot;
            particle[j].pot += pot;


            for(int k = 0; k < 3; k++)
            {
                double acc = 24 * dr[k] * (2/(r12 * dist_sqr) - 1/(r6 * dist_sqr));
                particle[i].acc[k] += acc;
                particle[j].acc[k] -= acc;
            }
        }
    }
}


/* This function calculates the total kinetic and total potential energy, averaged per particle. It also returns the instantaneous kinetic temperature. */
void calc_energies(part *particle, int Ntot, double epsilon, double *Ekin, double *Epot, double *temp)
{
    *Ekin = 0; *Epot = 0; /* we ensure that no contribution from the previous step enters our current calculation */

    for (int i = 0; i < Ntot; i++) {
        for (int j = 0; j < 3; j++)
        {
            *Ekin += 0.5 * particle[i].vel[j] * particle[i].vel[j];
        }
        *Epot += 0.5 * particle[i].pot;    /* since we assigned potential energies resulting from mutual interactions to both particles, we need to divide by 2 in order to avoid overcounting */
    }
    *temp = 2 * *Ekin * epsilon / (3 * Ntot);
}


/* This function rescales the velocities by a prescribed factor vscale */
void rescale_velocities(part *particle, int Ntot, double target_Ekin)
{
    double Ekin = 0;
    for (int i = 0; i < Ntot; i++)
    for (int k = 0; k < 3; k++)
    Ekin += 0.5 * particle[i].vel[k] * particle[i].vel[k];

    double vscale = sqrt(target_Ekin/Ekin);
    for (int i = 0; i < Ntot; i++)
    for (int k = 0; k < 3; k++)
    particle[i].vel[k] *= vscale;
}


/* main driver routine */
int main(int argc, char **argv)
{
    /* declaration of some constants to characterize our simulation */
    const int
    N1d = 8,
    Ntot = N1d * N1d * N1d,
    Nstep = 60000;
    const double
    dt = 0.01;

    /* declaration of constants that characterize the system */
    const double
    target_temp = 30,
    epsilon = 120,
    sigma = sqrt(target_temp/epsilon),
    target_Ekin = 3.0/2 * Ntot * target_temp/epsilon,
    mean_spacing = 5,
    boxsize = N1d * mean_spacing,
    rcut = 10;

    /* finally, we introduce the dynamic variables of our simulation */
    double Ekin, Epot, temp;

    /* allocate storage for our particles */
    part *particle = calloc(Ntot, sizeof(part));

    /* let's initialize the particles */
    initialize_system(particle, N1d, mean_spacing, sigma);
    /* calculate the forces at t=0 */
    calc_forces(particle, Ntot, boxsize, rcut);

    /* create an output file */
    char fname[20];
    sprintf(fname, "%dK.txt", (int) target_temp);
    FILE *fd = fopen(fname, "w");

    /* now we carry out time integration steps using the leapfrog */
    for(int step = 0; step <= Nstep; step++)
    {
        /* measure energies and output this to the file and screen */
        calc_energies(particle, Ntot, epsilon, &Ekin, &Epot, &temp);
        fprintf(fd, "%6d\t%10g\t%10g\t%10g\t%10g\n", step, Ekin, Epot, Ekin + Epot, temp);
        if(100 * step % Nstep == 0)
        printf("%6d\t%10g\t%10g\t%10g\t%10g\n", step, Ekin, Epot, Ekin + Epot, temp);

        kick(particle, Ntot, 0.5 * dt);

        drift(particle, Ntot, boxsize, dt);

        calc_forces(particle, Ntot, boxsize, rcut);

        kick(particle, Ntot, 0.5 * dt);

        /* every 100-th step, we call the rescaling routine */
        if((step % 100) == 0)
        rescale_velocities(particle, Ntot, target_Ekin);
    }
}
