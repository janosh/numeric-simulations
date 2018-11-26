/*============================================================================*/
/*! \file kelvin.c
 *  \brief Problem generator for KH instability. 
/*============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"


/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
 
  /* Read problem parameters */

  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  double u, v, w, x1, x2, x3;
  double dens1 = par_getd("problem", "density1");
  double dens2 = par_getd("problem", "density2");
  double v1 = par_getd("problem", "vx1");
  double v2 = par_getd("problem", "vx2");
  double pressure = par_getd("problem", "pressure");
  double ampl =  par_getd("problem", "ampl");
  double kwave = 2 * 2.0 * M_PI / 1.0; 


  for (k=ks; k<=ke; k++) 
    for (j=js; j<=je; j++) 
      for (i=is; i<=ie; i++) 
	{
	  cc_pos(pGrid, i, j, k, &x1, &x2, &x3); 
          /* x1 and x2 now hold the x and y coordinates of the cell centre, respectively */

	  /*** HERE CALCULATE YOUR DATA FOR THE MESH CELL (x1, x2, x3)*/

	  double dens = .... ;
	  double vx =   .... ;
	  double vy =   .... ;

	  /************************************************************/

	  pGrid->U[k][j][i].d = dens;

	  pGrid->U[k][j][i].M1 = dens * vx;
	  pGrid->U[k][j][i].M2 = dens * vy;
	  pGrid->U[k][j][i].M3 = 0; 

	  pGrid->U[k][j][i].E = pressure / Gamma_1 + 0.5 * dens * (vx * vx + vy * vy);            
	}
}

 


/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name)
{
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
  return;
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}
