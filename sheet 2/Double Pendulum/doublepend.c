#include <stdlib.h>
#include <math.h>
#include <stdio.h>

const double    delta_t = 0.05, period = 100,
                g = 1.0, m1 = 0.5, m2 = 1.0, l1 = 2.0, l2 = 1.0,
                phi1_init = 50.0 / 180.0 * M_PI, phi2_init = -119.0 / 180.0 * M_PI, phi1_dot_init = 0, phi2_dot_init = 0,
                q1_init = 0, q2_init = 0;

double y[4] = {phi1_init, phi2_init, q1_init, q2_init}, f[4], M_inv[4];

double Ecalc() {
    double hamilton = y[2] * f[0] + y[3] * f[1] - (m1 + m2)/2 * l1 * l1 * f[0] * f[0] - m2/2 * l2 * l2 * f[1] * f[1] - m2 * l1 * l2 * f[0] * f[1] * cos(y[0] - y[1]) + (m1 + m2) * g * l1 * (1 - cos(y[0])) + l2 * (1 - cos(y[1]));
    return hamilton;
}

void fcalc() {
    // matrix inverter
    double det_M = m2 * l1 * l1 * l2 * l2 * (m1 + m2 * pow(sin(y[0] - y[1]), 2));
    M_inv[0] = m2 * l2 * l2 / det_M;
    M_inv[1] = -m2 * l1 * l2 * cos(y[0] - y[1])/ det_M;
    M_inv[2] = M_inv[1];
    M_inv[3] = (m1 + m2) * l1 * l1 / det_M;
    // ODE RHS calculator
    f[0] = M_inv[0] * y[2] + M_inv[2] * y[3];
    f[1] = M_inv[2] * y[2] + M_inv[3] * y[3];
    f[2] = -m2 * l1 * l2 * f[0] * f[1] * sin(y[0] - y[1]) - (m1 + m2) * g * l1 * sin(y[0]);
    f[3] = m2 * l1 * l2  * f[0] * f[1] * sin(y[0] - y[1]) - m2 * g * l2 * sin(y[1]);
}

int main() {
    const double E0 = Ecalc();
    double E_err, time_passed = 0;
    FILE *error = fopen("error.dat", "w");
    FILE *movie = fopen("movie.dat", "w");
    // time-loop
    for (int step = 0; step * delta_t < period; step++) {
        double stor[4], k1[4], k2[4], k3[4], k4[4];
        /*
        // Runge-Kutta 2 routines
        fcalc();
        for (int i = 0; i < 4; i++) {
            k1[i] = delta_t * f[i];
            stor[i] = y[i];
            y[i] = y[i] + k1[i]/2;
        }
        fcalc();
        for (int i = 0; i < 4; i++) {
            k2[i] = delta_t * f[i];
            y[i] = stor[i] + k2[i];
        }
        */
        // Runge-Kutta 4 routines
        fcalc();
        for (int i = 0; i < 4; i++) {
            k1[i] = delta_t * f[i];
            stor[i] = y[i];
            y[i] = y[i] + k1[i]/2;
        }
        fcalc();
        for (int i = 0; i < 4; i++) {
            k2[i] = delta_t * f[i];
            y[i] = stor[i] + k2[i]/2;
        }
        fcalc();
        for (int i = 0; i < 4; i++) {
            k3[i] = delta_t * f[i];
            y[i] = stor[i] + k3[i];
        }
        fcalc();
        for (int i = 0; i < 4; i++) {
            k4[i] = delta_t * f[i];
            y[i] = stor[i] + (k1[i] + k4[i])/6 + (k2[i] + k3[i])/3;
        }
        // phi1 and phi2 readout
        fwrite(&y[0], sizeof(double), 1, movie);
        fwrite(&y[1], sizeof(double), 1, movie);
        // energy error and passed time readout
        time_passed = step * delta_t;
        fwrite(&time_passed, sizeof(double), 1, error);
        E_err = (Ecalc() - E0)/E0;
        fwrite(&E_err, sizeof(double), 1, error);
    }
}
