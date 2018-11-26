#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>


int main (void)
{
  double x = 1, z = 0;
  mpf_t y;

  mpf_init2 (y, 512); /* initialize to 512 bit precision */
  mpf_set_d (y, x); /* set y = x */
  z = mpf_get_d(y); /* set z = y */

  printf("This should be 1: %g\n", z);
}
