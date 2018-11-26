#include <stdlib.h>
#include <math.h>
#include <stdio.h>

const double k_B = 1.38e-16,
Lambda_0 = 1e-22,
T_0 = 2e4,
alpha = 10.0,
beta = -0.5,
n_H = 1.0,
T_init = 1e7,
T_fin = 6e3,
Delta_t = 8e10,
prefactor = -2/(3 * k_B) * n_H * Lambda_0;

double fcalc(const double temp) {
    if (temp <= T_0) {
        return prefactor * pow(temp/T_0, alpha);
    }
    else {
        return prefactor * pow(temp/T_0, beta);
    }
}

void stepper(double *temp, const double Delta_t) {
    double k_1 = fcalc(*temp) * Delta_t;
    double k_2 = fcalc(*temp + k_1/2) * Delta_t;
    *temp = *temp + k_2;
}

void adaptive_stepper(double *temp, double *adap_Delta_t, double *time_passed) {
    // one full step
    double k_1 = fcalc(*temp) * *adap_Delta_t;
    double k_2 = fcalc(*temp + k_1/2) * *adap_Delta_t;
    double one_step_temp = *temp + k_2;
    // two half steps
    double k_3 = fcalc(*temp + k_1/4) * *adap_Delta_t/2;
    double two_step_temp = *temp + k_3;
    double k_4 = fcalc(two_step_temp) * *adap_Delta_t;
    double k_5 = fcalc(two_step_temp + k_4/4) * *adap_Delta_t/2;
    two_step_temp = two_step_temp + k_5;
    double Delta_T = fabs(two_step_temp - one_step_temp);
    if (Delta_T > 50) {
        *adap_Delta_t /= 2;
        return;
    }
    if (Delta_T < 0.01) {
        *adap_Delta_t *= 2;
    }
    *temp = two_step_temp;
    *time_passed += *adap_Delta_t;
}

int main() {
    FILE *data = fopen("temp.dat", "w");
    double temp = T_init;
    int step = 0;
    while (temp > T_fin && step < 1e6) {
        stepper(&temp, Delta_t);
        double time_passed = step * Delta_t;
        fwrite(&time_passed, sizeof(double), 1, data);
        fwrite(&temp, sizeof(double), 1, data);
        ++step;
    }
    data = fopen("adaptive_temp.dat", "w");
    temp = T_init; step = 0;
    double time_passed = 0.0, adap_Delta_t = Delta_t;
    while (temp > T_fin && step < 1e6) {
        adaptive_stepper(&temp,&adap_Delta_t,&time_passed);
        fwrite(&time_passed, sizeof(double), 1, data);
        fwrite(&temp, sizeof(double), 1, data);
        ++step;
    }
}
