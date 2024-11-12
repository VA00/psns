#include <math.h>
#include "thermal.h"
#include "thermal_tbl.h" 
#include <interp.h>
#include "plasmaT_R_table.c"
#include "plasmaT_Q_table.c"
#include "plasmaL_R_table.c"
#include "plasmaL_Q_table.c"

double plasmaL_Q_tbl(double lg_T, double lg_rhoYe, int flavour)
{
  double lg_Q=-99.9, *tbl;

  tbl = (double *) plasmaLtable_Q[flavour-1];

  lg_Q =  bilinear_interp(lg_T, lg_rhoYe, tbl, 
    N_T_thermal+1, N_rhoYe_thermal+1, lg_T_1_thermal, lg_T_2_thermal, lg_rhoYe_1_thermal, lg_rhoYe_2_thermal, 0, -99.9);

  return pow(10.0,lg_Q);

}

double plasmaL_R_tbl(double lg_T, double lg_rhoYe, int flavour)
{
  double lg_R = -89.9, *tbl;

  tbl = (double *) plasmaLtable_R[flavour-1]; 

  lg_R =  bilinear_interp(lg_T, lg_rhoYe, tbl, 
    N_T_thermal+1, N_rhoYe_thermal+1, 
	lg_T_1_thermal, lg_T_2_thermal, lg_rhoYe_1_thermal, lg_rhoYe_2_thermal, 0, -89.9);

  return pow(10.0,lg_R);  

}



double plasmaT_Q_tbl(double lg_T, double lg_rhoYe, int flavour)
{

  double lg_Q=-99.9, *tbl;

  tbl = (double *) plasmaTtable_Q[flavour-1];

  lg_Q =  bilinear_interp(lg_T, lg_rhoYe, tbl, 
    N_T_thermal+1, N_rhoYe_thermal+1, lg_T_1_thermal, lg_T_2_thermal, lg_rhoYe_1_thermal, lg_rhoYe_2_thermal, 0, -99.9);

  return pow(10.0,lg_Q);

}

double plasmaT_R_tbl(double lg_T, double lg_rhoYe, int flavour)
{
  double lg_R=-89.9, *tbl;

  tbl = (double *) plasmaTtable_R[flavour-1]; 


  lg_R =  bilinear_interp(lg_T, lg_rhoYe, tbl, 
    N_T_thermal+1, N_rhoYe_thermal+1, lg_T_1_thermal, lg_T_2_thermal, lg_rhoYe_1_thermal, lg_rhoYe_2_thermal, 0, -89.9);

  return pow(10.0,lg_R);

}

