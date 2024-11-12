#include <math.h>
#include "thermal.h"
#include "thermal_tbl.h" 
#include <interp.h>
/*
This file has been part of pair_tbl.c; due to very long compile
time using icc pair_tbl.c has been split into three:
pair_tbl1.c (nue)
pair_tbl2.c      (nuebar)
pair_tbl3.c              (numu)
*/
#include "pair_table_nue.c" 



double pair_tbl_nue(double Enu, double lg_T, double lg_rhoYe)
{

/*
NOTE: this code below is working as expected with libtool
FIXME: restore original idea (passing pointer, rather tham multiple 
return(s))   and test if this is working with libtool
*/
//   if( flavour==1 ) tbl = &thermal_Enue[0][0][0];
//   if( flavour==2 ) tbl = &thermal_Enuebar[0][0][0];
//   if( flavour==3 ) tbl = &thermal_Enumu[0][0][0];
//   if( flavour==4 ) tbl = &thermal_Enumu[0][0][0];
// 
//   if( (flavour<1) || (flavour>4) ) return 0.0;

//   if( flavour==1 ) printf("flavour==%d",flavour); ;
//   if( flavour==2 ) printf("flavour==%d",flavour);
//   if( flavour==3 ) printf("flavour==%d",flavour);
//   if( flavour==4 ) printf("flavour==%d",flavour);

  return trilinear_interp(Enu, lg_T, lg_rhoYe, &thermal_Enue[0][0][0],
         N_Enu_thermal+1,  N_T_thermal+1, N_rhoYe_thermal+1, 
         0.0, 20.0,
         lg_T_1_thermal, lg_T_2_thermal, lg_rhoYe_1_thermal, lg_rhoYe_2_thermal,
	 0,0.0); //OUT_OF_DOMAIN_VALUE=0.0 (do not extrapolate)
}
