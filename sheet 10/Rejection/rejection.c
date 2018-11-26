#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


double p(double x)
{
	double p = 1/(pow(x - 2, 4) + pow(sin(x - 3), 8));
	return p;
}


double f(double x)
{
	if(x < 0 || x > 5) return 0;

	int i;
	double xn[] = {0.0, 1.8, 2.35, 3.0, 5.0}, yn[] = {0.01, 0.15, 2.5, 0.1, 0.002};

	for (int j = 0; j < 4; j++)
	if (x >= xn[j] && x <= xn[j+1]) i = j;

	double f = (yn[i+1] - yn[i])/(xn[i+1] - xn[i]) * (x - xn[i]) + yn[i];
	return f;
}


double cdf()
{
	double x, z = 1.81975 * drand48();    // 1/1.81975 = normalization factor of PDF f(x)

	if (z < 0.144)
	return x = ( - 0.01 + sqrt(0.0001 + 4 * 0.0388889 * z))/(2 * 0.0388889);

	else if (z < 0.87275)
	return x = (0.00454545 * 1659 + sqrt(0.00454545 * 1659 * 0.00454545 * 1659 - 4 * 0.00454545 * 470 * (6.79582 - z)))/(2 * 0.00454545 * 470);

	else if (z < 1.71775)
	return x = -( - 3.69231 * 3.02708 + sqrt(3.69231 * 3.02708 * 3.69231 * 3.02708 - 4 * 0.5 * 3.69231 * (15.1976 + z)))/3.69231;

	else
	return x = -( - 0.001 * 247 + sqrt(0.001 * 247 * 0.001 * 247 + 4 * 0.001 * 24.5 * (1.19725 - z)))/(2 * 0.001 * 24.5);
}


double calc_uniform_rejecrate(double xlength, double ylength, int sample_length)
{
	int acc_count = 0, tot_count = 0;
	FILE *sample_file = fopen("rejection_sample.txt", "w");

	while (acc_count < sample_length) {
		double x = xlength * drand48(), y = ylength * drand48();
		if (y <= p(x)) {
			fprintf(sample_file, "%g\n", x);
			acc_count++;
		}
		tot_count++;
	}

	double rejec_rate = 1 - (double) acc_count/tot_count, p_0 = xlength * ylength * (1 - rejec_rate);
	printf("\nRegular rejection method\n\trejection rate r = %g\n\tnormalisation p_0 = %g\n\n", rejec_rate, p_0);
	return p_0;
}


void calc_nonuniform_rejecrate(int sample_length, double p_0)
{
	int acc_count = 0, tot_count = 0;
	FILE *sample_file = fopen("rejection_aux_sample.txt", "w");

	while (acc_count < sample_length) {
		double x = cdf(), y = f(x) * drand48();
		if (y <= p(x)/p_0) {
			fprintf(sample_file, "%g\n", x);
			acc_count++;
		}
		tot_count++;
	}

	double rejec_rate = 1 - (double) acc_count/tot_count;
	printf("Rejection method with auxiliary function\n\trejection rate r = %g\n\n", rejec_rate);
}


int main()
{
	double xmin = 0, xmax = 5, xlength = xmax - xmin;
	double ymin = 0, ymax = 30.5, ylength = ymax - ymin;
	int sample_length = 1e6;

	double p_0 = calc_uniform_rejecrate(xlength, ylength, sample_length);

	calc_nonuniform_rejecrate(sample_length, p_0);
}
