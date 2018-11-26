#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <time.h>


void gaussseidel_step(double *phi, double *rho, double *res, int N, double h)
{
  double *phinew = malloc(N * N * sizeof(double));
  int i,j,k, im, jm, ip, jp;

  double phisum = 0;

  for(i=0; i < N; i++)
    for(j=0; j < N; j++)
      {
	phinew[i*N + j] = phi[i*N + j];

	im = i-1;
	jm = j-1;
	ip = i+1;
	jp = j+1;
	if(im < 0)
	  im += N;
	if(jm < 0)
	  jm += N;
	if(ip >= N)
	  ip -= N;
	if(jp >= N)
	  jp -= N;

  	res[i*N + j] = 4 * M_PI * h * h * rho[i*N + j] 
                                 - ( - 4 * phi[i*N + j]
				     + phi[ip*N + j] + phi[i*N + jp] + 
				       phi[im*N + j] + phi[i*N + jm]); 
      }

  for(k=0; k<2; k++)
    for(i=0; i < N; i++)
      for(j=(((i&1)==k)?1:0); j < N; j+=2)
  

  //    for(i=0; i < N; i++)
  //    for(j=0; j < N; j++)
	{
	  im = i-1;
	  jm = j-1;
	  ip = i+1;
	  jp = j+1;
	  if(im < 0)
	    im += N;
	  if(jm < 0)
	    jm += N;
	  if(ip >= N)
	    ip -= N;
	  if(jp >= N)
	    jp -= N;
	  
	  phinew[i*N + j] = 0.25 * (
				    phinew[ip*N + j] + phinew[i*N + jp] + 
				    phinew[im*N + j] + phinew[i*N + jm] - 
				    4 * M_PI * h * h * rho[i*N + jm]);

	  phisum += phinew[i*N + j];
	}
  
  phisum=0;
	
  for(i=0; i < N; i++)
    for(j=0; j < N; j++)
      phi[i*N + j] = phinew[i*N + j] - phisum / (N*N);
  
  free(phinew);
}


double norm_of_residual(double *res, int N)
{
  int i, j;
  double sum = 0;

  for(i=0; i < N; i++)
    for(j=0; j < N; j++)
      sum += res[i*N + j]*res[i*N + j];
}


int main(int argc, char **argv)
{
  int i, j;
  int N = 128;
  double L = 1.0;
  double h = L / N;

  /* allocate some storage for the image, and then read it */

  double *phi =   malloc(N * N * sizeof(double));
  double *rho =   malloc(N * N * sizeof(double));
  double *res =   malloc(N * N * sizeof(double));

  double sum = 0;

  for(i=0; i < N; i++)
    for(j=0; j < N; j++)
      {
	phi[i*N + j] = 0;  

	double dx =  (i - N/2 + 0.5)*h;
	double dy =  (j - N/2 + 0.5)*h;
	double r2 = dx*dx + dy*dy;
	double eta = 0.1*L;

        rho[i*N + j] = 10.0 * exp(-r2 / (2*eta*eta));  
	
	sum += rho[i*N + j];
      }

  for(i=0; i < N; i++)
    for(j=0; j < N; j++)
      {
	rho[i*N + j] -= sum / (N * N);	
      }


  FILE *fd = fopen("res_gausseidel_redblack.txt", "w");

  for(i=0; i<2000; i++)
    {
      gaussseidel_step(phi, rho, res, N, h);
      
      double r = norm_of_residual(res, N);

      printf("iter=%d:  residual=%g\n", i, r);

      fprintf(fd, "%g\n", r);
    }

  double phisum = 0;

  fprintf(fd, "%d\n", N);
  for(i=0; i<N; i++)
    for(j=0; j < N; j++)
      {
	fprintf(fd, "%g\n", phi[i*N + j]);
	phisum += phi[i*N + j];
      }
  printf("phisum=%g\n", phisum);

  fclose(fd);
}
