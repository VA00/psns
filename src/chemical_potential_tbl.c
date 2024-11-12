#include "thermal.h"
#include "thermal_tbl.h" 
#include <interp.h>
#include <math.h>
#include "chemical_potential_table.c"

double chemical_potential_tbl(double lg_T, double lg_rhoYe)
{

  double mu = 0.0, *tbl;

  tbl = (double *) chemical_potential_table;
  

  mu = bilinear_interp(lg_T, lg_rhoYe, tbl, 
    N_T_mu+1, N_rhoYe_mu+1, lg_T_1_mu, lg_T_2_mu, lg_rhoYe_1_mu, lg_rhoYe_2_mu, 1, 0.0);

  return mu;

}

