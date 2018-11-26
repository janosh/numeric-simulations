#include <stdio.h>

int main() {
    float  x = 0.01;
    double y = x;
    double z = 0.01;
    int i = 10000 * x, j = 10000 * y, k = 10000 * z;

    printf("\nx = %.30f\ny = %.30f\nz = %.50f\n\n", x, y, z);
    printf("10000*x = %.30f\n10000*y = %.30f\n10000*z = %.30f\n\n", 10000*x, 10000*y, 10000*z);
    printf("i = %d\nj = %d\nk = %d\n\n", i, j, k);
}
