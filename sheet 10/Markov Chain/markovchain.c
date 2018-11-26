#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double p(double x) {
	return exp(-pow(x + 2 * cos(x) * cos(x), 2));
}


void generate_markovchain(int chain_length, double *chain) {
	int chain_pos = 0; chain[0] = 0;

	while (chain_pos < chain_length) {
		double
		x = chain[chain_pos] + (2 * drand48() - 1),
		p_ratio = p(x)/p(chain[chain_pos]),
		r = (1 > p_ratio) ? p_ratio : 1;
		chain_pos++;
		if (drand48() < r) chain[chain_pos] = x;
		else chain[chain_pos] = chain[chain_pos - 1];
	}

	FILE *markovchain_file = fopen("markovchain.txt", "w");
	for (int i = 0; i < chain_length; i++) {
		fprintf(markovchain_file, "%g\n", chain[i]);
	}
}

int compare(const void *a, const void *b) {
	if (*(double*) a <  *(double*) b) return -1;
	if (*(double*) a == *(double*) b) return 0;
	else return 1;
}


void count_distinct_entries(int chain_length, double *chain) {
	qsort(chain, chain_length, sizeof(double), compare);
	int num_dist_entr = 1;
	for (int i = 1; i < chain_length; i++) {
		if (chain[i] != chain[i-1]) num_dist_entr++;
	}
	printf("number of distinct entries in the markov chain: %d\n", num_dist_entr);
}


int main() {
	int chain_length = 1e6;
	double *chain = malloc(chain_length * sizeof(double));
	generate_markovchain(chain_length, chain);
	count_distinct_entries(chain_length, chain);
}
