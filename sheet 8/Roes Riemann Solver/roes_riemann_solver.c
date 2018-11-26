#include <math.h>
#include <stdio.h>
#include <stdlib.h>


void riemann_roe_isothermal(double rhoL, double uL, double vL, double rhoR, double uR, double vR, double cs, double *flux)
{
	double
	uavr = (sqrt(rhoL) * uL + sqrt(rhoR) * uR)/(sqrt(rhoL) + sqrt(rhoR)),
	vavr = (sqrt(rhoL) * vL + sqrt(rhoR) * vR)/(sqrt(rhoL) + sqrt(rhoR)),

	lambda1 = uavr - cs,
	lambda2 = uavr + cs,
	lambda3 = uavr,

	u1 = rhoR - rhoL,
	u2 = rhoR * uR - rhoL * uL,
	u3 = rhoR * vR - rhoL * vL,

	alpha1 = ((uavr + cs) * u1 - u2)/(2 * cs),
	alpha2 = (-(uavr - cs) * u1 + u2)/(2 * cs),
	alpha3 = u3 - vavr * u1,

	fluxL[3] = {rhoL * uL, rhoL * (uL * uL + cs * cs), rhoL * uL * vL},
	fluxR[3] = {rhoR * uR, rhoR * (uR * uR + cs * cs), rhoR * uR * vR},

	K1[3] = {1, uavr - cs, vavr},
	K2[3] = {1, uavr + cs, vavr},
	K3[3] = {0, 0, 1};

	for (int i = 0; i < 3; i++) {
		flux[i] = 0.5 * (fluxL[i] + fluxR[i] - (alpha1 * fabs(lambda1) * K1[i] + alpha2 * fabs(lambda2) * K2[i] + alpha3 * fabs(lambda3) * K3[i]));
	}


}


void sweep(double *rho, double *u, double *v, double dt, double dx, double cs, int N, int M, int *offset1, int *offset2)
{
	double Un[N][M][3], Unp1[N][M][3], fluxL[3], fluxR[3];

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			Un[i][j][0] = rho[M * i + j];
			Un[i][j][1] = rho[M * i + j] * u[M * i + j];
			Un[i][j][2] = rho[M * i + j] * v[M * i + j];
		}
	}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			int
			lcix = (i + offset1[0] + N) % N,	/* left cell index x */
			lciy = (j + offset1[1] + M) % M;	/* left cell index y */
			riemann_roe_isothermal(rho[M * lcix + lciy], u[M * lcix + lciy], v[M * lcix + lciy], rho[M * i + j], u[M * i + j], v[M * i + j], cs, fluxL);

			int
			rcix = (i + offset2[0] + N) % N,	/* right cell index x */
			rciy = (j + offset2[1] + M) % M;	/* right cell index y */
			riemann_roe_isothermal(rho[M * i + j], u[M * i + j], v[M * i + j], rho[M * rcix + rciy], u[M * rcix + rciy], v[M * rcix + rciy], cs, fluxR);

			for (int k = 0; k < 3; k++) {
				Unp1[i][j][k] = Un[i][j][k] + dt/dx * (fluxL[k] - fluxR[k]);
			}
		}
	}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			rho[M * i + j] = Unp1[i][j][0];
			u[M * i + j] = Unp1[i][j][1]/rho[M * i + j];
			v[M * i + j] = Unp1[i][j][2]/rho[M * i + j];
		}
	}
}


double max(double *array, int length)
{
	double max = array[0];

	for(int i = 1; i < length; i++)
	if(max < array[i]) max = array[i];

	return max;
}


int main(void) {

	/* Testing Roe’s approximate Riemann solver for the isothermal problem */
	double
	rhoL[3] = {1.0, 2.5, 2.0},
	uL[3] = {1.0, 2.0, -1.0},
	vL[3] = {2.0, 3.0, -2.0},

	rhoR[3] = {3.0, 1.0, 1.0},
	uR[3] = {1.0, -3.0, -1.0},
	vR[3] = {0.0, -2.0, 2.0},

	cs = 2.0,
	flux[3];

	for (int i = 0; i < 3; i++) {
		riemann_roe_isothermal(rhoL[i], uL[i], vL[i], rhoR[i], uR[i], vR[i], cs, flux);
		for (int j = 0; j < 3; j++) {
			printf("f%d = %g\n", j, flux[j]);
		}
	}

	/* Carrying out a multidimensional simulation */
	int
	N = 60, M = 30,
	offsetx1[2] = {-1,0}, offsetx2[2] = {1,0},
	offsety1[2] = {0,-1}, offsety2[2] = {0,1};

	double
	Lx = 3, Ly = 1.5,
	dx = Lx/N, dy = Ly/M,
	Tmax = 1.5,
	CCFL = 0.4,
	dt = CCFL * ((dx <= dy) ? dx : dy)/cs,
	*u = calloc(N * M, sizeof(double)),
	*v = calloc(N * M, sizeof(double)),
	*rho = calloc(N * M, sizeof(double));

	/* Initiializing the density */
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			if (fabs(dx * i - Lx/2) < Lx/4 && fabs(dy * j - Ly/2) < Ly/4) {
				rho[M * i + j] = 4;
			}
			else {
				rho[M * i + j] = 1;
			}
		}
	}

	double time = 0;
	for (int i = 0; time < Tmax; i++) {
		if (i % 100 == 0) printf("time = %g\t(dt = %.6g)\n", time, dt);

		sweep(rho, u, v, dt, dx, cs, N, M, offsetx1, offsetx2);
		sweep(rho, v, u, dt, dy, cs, N, M, offsety1, offsety2);

		time += dt;
		double maxu = max(u, N * M), maxv = max(v, N * M);
		dt = CCFL * ((dx <= dy) ? dx : dy)/(cs + ((maxu > maxv) ? maxu : maxv));
	}

	FILE *rho_data = fopen("rho.txt", "w");
	for (int j = 0; j < N; j++) {
		for (int k = 0; k < M; k++)
		fprintf(rho_data, "%g\t", rho[M * j + k]);
		fprintf(rho_data, "\n");
	}
}
