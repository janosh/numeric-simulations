#include <stdlib.h>		// included for drand48
#include <stdio.h>		// included for printf
#include <math.h>		// included for fmod

// range of the RANDU random number generator (2^31)
const double range = 2147483648;

// IBM's 'horrible' RANDU random number generator
double RANDU(double previous){
	double next = fmod(65539 * previous, range);
	return next;
}

int main(){
	// array's to hold x,y coordinates of 1000 pseudo-random points
	int npoints = 1000;
	double x[npoints], y[npoints];

	// seed variable
	const double seed = 1.0;
	// storage variable to hold the current value of RANDU between cycles of the for-loop
	double stor = RANDU(seed);

	// assign odd numbers in the RANDU sequence of random numbers to x-coordinates of points and the immediately following even number as the corresponding y-coordinate
	for (int i = 0; i < npoints; i++) {
		x[i] = stor/range;
		double temp = RANDU(stor);
		y[i] = temp/range;
		stor = RANDU(temp);
	}

	// write RANDU points to file
	FILE *randu = fopen("randu.txt", "w");
	for (int i = 0; i < npoints; i++) {
		fprintf(randu, "%-15g%-15g\n", x[i], y[i]);
	}


	// overwrite coordinates with drand48() values
	for (int i = 0; i < npoints; i++) {
		x[i] = drand48();
		y[i] = drand48();
	}

	// write drand48() points to file
	FILE *drand = fopen("drand.txt", "w");
	for (int i = 0; i < npoints; i++) {
		fprintf(drand, "%-15g%-15g\n", x[i], y[i]);
	}



	int points_in_zoom = 0;
	// variables to hold two successive RANDU numbers
	double rand1 = RANDU(seed), rand2 = RANDU(rand1);

	while (points_in_zoom < npoints) {
		if (0.2 <= rand1/range && rand1/range <= 0.2005 && 0.3 <= rand2/range && rand2/range <= 0.3005) {
			x[points_in_zoom] = rand1/range;
			y[points_in_zoom] = rand2/range;
			points_in_zoom++;
			printf("number of points in RANDU zoom: %d\n", points_in_zoom);
		}
		rand1 = RANDU(rand2), rand2 = RANDU(rand1);
	}

	FILE *randu_zoom = fopen("randu_zoom.txt", "w");
	for (int i = 0; i < npoints; i++) {
		fprintf(randu_zoom, "%-15g%-15g\n", x[i], y[i]);
	}



	points_in_zoom = 0;

	while (points_in_zoom < npoints) {
		// variables to hold two successive drand48 numbers
		rand1 = drand48(), rand2 = drand48();
		if (0.2 <= rand1 && rand1 <= 0.2005 && 0.3 <= rand2 && rand2 <= 0.3005) {
			x[points_in_zoom] = rand1;
			y[points_in_zoom] = rand2;
			points_in_zoom++;
			printf("number of points in drand zoom: %d\n", points_in_zoom);
		}
	}

	FILE *drand_zoom = fopen("drand_zoom.txt", "w");
	for (int i = 0; i < npoints; i++) {
		fprintf(drand_zoom, "%-15g%-15g\n", x[i], y[i]);
	}
}
