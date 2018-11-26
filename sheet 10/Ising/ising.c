#include <stdlib.h>
#include <stdio.h>
#include <math.h>


double calc_energy(int *lattice, int N, int x, int y, double beta)
{
	double energy = 0;

	for (int i = -1; i <= 1; i += 2) {
		int nn_xpos = (x + i + N) % N;	// implements periodic boundary for nearest-neighbor x-position
		energy += beta/2 * (1 - lattice[N * x + y] * lattice[N * nn_xpos + y]);
	}

	for (int j = -1; j <= 1; j += 2) {
		int nn_ypos = (y + j + N) % N;	// implements periodic boundary for nearest-neighbor y-position
		energy += beta/2 * (1 - lattice[N * x + y] * lattice[N * x + nn_ypos]);
	}

	return energy;
}


void initialize_lattice(int *lattice, int N)
{
	for (int i = 0; i < N * N; i++) {
		if (drand48() < 0.5) lattice[i] = 1;
		else lattice[i] = -1;
	}
}


void sweep_lattice(int *lattice, int N, double beta)
{
	// chess-board iteration scheme
	for(int j = 0; j < 2; j++)
	for(int x = 0; x < N; x++)
	for(int y = (((x & 1) == j) ? 1 : 0); y < N; y += 2)
	{
		lattice[N * x + y] = 1;	// prepare current lattice site for calculation of Eplus
		double Eplus = calc_energy(lattice, N, x, y, beta);	// interaction energy at current lattice site if spin up
		lattice[N * x + y] = -1;	// prepare current lattice site for calculation of Eminus
		double Eminus = calc_energy(lattice, N, x, y, beta);	// interaction energy at current lattice site if spin down

		double pplus = exp(-Eplus)/(exp(-Eplus) + exp(-Eminus));
		if (pplus > drand48()) lattice[N * x + y] = 1;	// reset lattice site to spin up with probability pplus, otherwise leave it at spin down
	}
}


void thermalize_lattice(int *lattice, int N, int relax_steps, double beta)
{
	for(int i = 0; i < relax_steps; i++)	// thermalization loop
	sweep_lattice(lattice, N, beta);
}


void writeout_lattice(int *lattice, int N, double beta)
{
	char filename[100];	// C-style string to hold dynamic filename
	snprintf(filename, 100 * sizeof(char), "lattice%02d.txt", (int) (10 * beta));
	FILE *lattice_file = fopen(filename, "w");
	for (int i = 0; i < N * N; i++) {
		if(i % N == 0) fprintf(lattice_file, "\n");
		fprintf(lattice_file, "%d\t", lattice[i]);
	}
}


void calc_magnetization(int *lattice, int N, int avr_steps, double beta, FILE *Mdata)
{
	double avr_magnetization = 0, crt_magnetization;

	for(int i = 0; i < avr_steps; i++) {	// averaging loop to measure magnetization

		crt_magnetization = 0;

		for (int i = 0; i < N * N; i++)
		crt_magnetization += (double) lattice[i]/(N * N);

		avr_magnetization += fabs(crt_magnetization)/avr_steps;

		sweep_lattice(lattice, N, beta);
	}

	printf("M_avr(beta = %.2f) = %g\n", beta, avr_magnetization);
	fprintf(Mdata, "%-10g\t%-10g\n", beta, avr_magnetization);	// write out average magnetizations and correspoding temperature to file for plotting (%-10g for left justification in a 10 space-wide column)
}


int main()
{
	int N = 32, *lattice = malloc(N * N * sizeof(int));	// length of lattice and lattice itself
	double beta[] = {100.0, 10.0, 1.8, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 0.95, 0.9, 0.875, 0.85, 0.825, 0.8, 0.775, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.1};	// values of beta, i.e. temperature
	int relax_steps = 3000;	// number of times the lattice is sweeped in chessboard pattern to achieve thermalization
	int avr_steps = 1000;	// number of times the lattice is sweeped to measure and the mgneization and then average over all measurements
	FILE *Mdata = fopen("magnetization.txt", "w");

	for(int i = 0; i < sizeof(beta)/sizeof(double); i++)	// beta loop
	{
		initialize_lattice(lattice, N);

		thermalize_lattice(lattice, N, relax_steps, beta[i]);

		writeout_lattice(lattice, N, beta[i]);

		calc_magnetization(lattice, N, avr_steps, beta[i], Mdata);
	}
}
