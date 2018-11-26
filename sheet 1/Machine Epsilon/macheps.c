#include <stdio.h>
#include <float.h>

int main() {
    // approximate calculations:
    float f_epsilon = 1;
    while (1 + f_epsilon/2 != 1) {
        f_epsilon /= 2;
    }
    double d_epsilon = 1;
    while (1 + d_epsilon/2 != 1) {
        d_epsilon /= 2;
    }
    long double l_epsilon = 1;
    while (1 + l_epsilon/2 != 1) {
        l_epsilon /= 2;
    }

    printf("approximations:\n\tfloat epsilon = \t%e\n\tdouble epsilon = \t%e\n\tlong double epsilon = \t%Le\n", f_epsilon, d_epsilon, l_epsilon);

    // what would be done in actual programming:
    printf("according to <float.h>:\n\tfloat epsilon = \t%e\n\tdouble epsilon = \t%e\n\tond double epsilon = \t%Le\n", FLT_EPSILON, DBL_EPSILON, LDBL_EPSILON);
}
