#include <math.h>
#include "thermal.h"
#include "thermal_tbl.h" 
#include <interp.h>
#include "pair_Q_table.c"
#include "pair_R_table.c"
/*
nue - electron neutrino
nuebar - electron antineutrino
numu - average of muon neutrino and muon antineutrino spectrum
NOTE: tau neutrino/antineutrino  spectrum = muon neutrino/antineutrino spectrum
NOTE: difference between numu and numubar is much smaller than
difference between nue and nuebar, therefore average is used.
*/


/* these functions provide convenient access to pre-calculated tables generated
using ,,normal'' functions (without _tbl suffix). Possible speed-up
is 100x-1000x. Much more memory is possible used, but no cache memory usage
has been tested or estimated. 

Main difference is:
pair_Q(kT, mu, flavour)      versus pair_Q_tbl(lg_T, lg_rhoYe, flavour)
i.e. log10 of the temperature T (Kelvins) and electron density rhoYe (g/cc) is used as an input
rather than temperature kT (MeV) and chemical potential mu (MeV). Additional spped-up
is achieved if mu is not available.

 flavour 1 - nue, 2 - nuebar, 3- numu, 4 -numubar */
double pair_Q_tbl(double lg_T, double lg_rhoYe, int flavour)
{

  double lg_Q = -99.9, *tbl;
  
  tbl = (double *) pairtable_Q[flavour-1];

  lg_Q = bilinear_interp(lg_T, lg_rhoYe, tbl, 
    N_T_thermal+1, N_rhoYe_thermal+1, lg_T_1_thermal, lg_T_2_thermal, lg_rhoYe_1_thermal, lg_rhoYe_2_thermal, 0, -99.9);

  return pow(10.0,lg_Q);

}


double pair_R_tbl(double lg_T, double lg_rhoYe, int flavour)
{

  double lg_R = -89.9, *tbl;

  tbl = (double *) pairtable_R[flavour-1];
  

 lg_R = bilinear_interp(lg_T, lg_rhoYe, tbl, 
    N_T_thermal+1, N_rhoYe_thermal+1, lg_T_1_thermal, lg_T_2_thermal, lg_rhoYe_1_thermal, lg_rhoYe_2_thermal, 0, -89.9);

  return pow(10.0,lg_R);
}

double pair_tbl(double Enu, double lg_T, double lg_rhoYe, int flavour)
{

/*
NOTE: this code below is working as expected with libtool
FIXME: restore original idea (passing pointer, rather than multiple 
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

  if( flavour==1 )   return pair_tbl_nue   (Enu,  lg_T, lg_rhoYe);

  if( flavour==2 )   return pair_tbl_nuebar(Enu,  lg_T, lg_rhoYe);

  if( (flavour==3) || (flavour==4) )   
                     return pair_tbl_numu   (Enu, lg_T, lg_rhoYe);



}

