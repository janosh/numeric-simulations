#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

double mp_integral = 0;


double func(double x) {
	double f = 1.5 * (1 - x * x);
	return f;
}


void midpoint(int dim, const int dim_stor, const int intervals, double *x) {
    if(dim == 0) {
		double f = 1;
		for (int i = 0; i < dim_stor; i++) f *= func(x[i]);
		mp_integral += f/pow(intervals, dim_stor);
	}
    else {
        for (int i = 0; i < intervals; i++) {
			x[dim - 1] = (i + 0.5)/intervals;
			midpoint(dim - 1, dim_stor, intervals, x);
		}
    }
}


void montecarlo(int dim, int samples) {
	clock_t time = clock();
	double result = 0, x, f;
	for (int j = 0; j < samples; j++) {
		f = 1;
		for (int k = 0; k < dim; k++) {
			x = drand48();
			f *= func(x);
		}
		result += f/samples;
	}
	time = clock() - time;
	printf("I_mc(dim = %d) = %10.7g\t(took %g sec)\n", dim, result, (double) time/CLOCKS_PER_SEC);
}


int main() {
	const int dim = 12, intervals = 6, samples = 2e4;

	for (int i = 1; i <= dim; i++) {
		clock_t time = clock();
		double x[i];
		midpoint(i, i, intervals, x);
		time = clock() - time;
		printf("I_mp(dim = %d) = %10.7g\t(took %g sec)\n", i, mp_integral, (double) time/CLOCKS_PER_SEC);
		mp_integral = 0;
		montecarlo(i, samples);
	}
}
