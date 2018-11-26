#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gmp.h>

int compare (const void *a, const void *b) {
	const double *da = (const double *) a;
	const double *db = (const double *) b;
	return (fabs(*da) > fabs(*db)) - (fabs(*da) < fabs(*db));
}

void forwardSum(double *doubleArray, int nDoubles) {
	double forwardSum = 0;
	for(int i = 0; i < nDoubles; ++i) {
		forwardSum += doubleArray[i];
	}
	printf("Forward sum of numbers.dat:\t\t\t%f\n", forwardSum);
}

void backwardSum(double *doubleArray, int nDoubles) {
	double backwardSum = 0;
	for(int i = nDoubles; i > 0; --i){
		backwardSum += doubleArray[i];
	}
	printf("Backward sum of numbers.dat:\t\t\t%f\n", backwardSum);
}

void sortedSum(double *doubleArray, int nDoubles) {
	double sortedSum = 0;
	for(int i = 0; i < nDoubles; ++i) {
		sortedSum += doubleArray[i];
	}
	printf("Sum sorted by magnitude:\t\t\t%f\n", sortedSum);
}

void longSortedSum(double *doubleArray, int nDoubles) {
	long double longSortedSum = 0;
	for(int i = 0; i < nDoubles; ++i) {
		longSortedSum += doubleArray[i];
	}
	printf("Long double precision sum sorted by magnitude:\t%Lf\n", longSortedSum);
}

void gmpSortedSum(double *doubleArray, int nDoubles) {
	mpf_t gmpSum;
	mpf_init2 (gmpSum, 512);
	for(int i = 0; i < nDoubles; ++i) {
		mpf_t arrayEntry;
		mpf_init2 (arrayEntry, 512);
		mpf_set_d (arrayEntry,doubleArray[i]);
		mpf_add(gmpSum,gmpSum,arrayEntry);
	}
	double dGmpSum = mpf_get_d(gmpSum);
	// mpf_clear (gmpSum); and mpf_clear (arrayEntry); not required becuase the program is about to exit anyway.
	printf("512 bit precision sum sorted by magnitude:\t%f\n", dGmpSum);
}

int main() {

	FILE *ptr_to_datafile;
	ptr_to_datafile = fopen("numbers.dat", "r");
	int nDoubles;
	fread(&nDoubles, sizeof(int), 1, ptr_to_datafile);
	printf("Number of doubles in numbers.dat: %d\n", nDoubles);

	double doubleArray[nDoubles];
	fread(doubleArray, sizeof(double), nDoubles, ptr_to_datafile);

	forwardSum(doubleArray, nDoubles);
	backwardSum(doubleArray, nDoubles);

	qsort(doubleArray, nDoubles, sizeof(double), compare);

	sortedSum(doubleArray, nDoubles);
	longSortedSum(doubleArray, nDoubles);
	gmpSortedSum(doubleArray, nDoubles);
}
