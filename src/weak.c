#include <stdio.h>
#include <math.h>
#include <fermidirac.h>
#include <interp.h>
#include <nse.h>
#include "weak.h"
#include "thermal.h"
#include "../psns.h"
#include "ffn_table_nue.c"
#include "ffn_table_nuebar.c"
#include "ffn_data.c"
#define DEBUG 0

/*
TODO: more inteligent integration e.g. Gauss-Kronrod (Laguerre) or Tanh-Sinh instead of trapezoidal rule
TODO: more efficient find root, now bisection
TODO: various optimizations (f. inlining, profiling, compiler options)
TODO: in NSE_neutrino_spectrum_table() most of interpolations is used only
      to select between capture/decay. It might be efficient to implement this
      as e.g. sign of eff_rate.

FFN weak reates tbl_number=0 for combined rate, tbl_number=1 for mean
neutrino energy, 2 for decay rate, 3 for capture rate 
*/


double 
ffn_rate_nue(double T9, double lg_rhoYe, int nuclide, int tbl_number)
{
  return ffn_bilinear(T9,lg_rhoYe,iso_nue[tbl_number][nuclide]);
}

double 
ffn_rate_nuebar(double T9, double lg_rhoYe, int nuclide, int tbl_number)
{
  return ffn_bilinear(T9,lg_rhoYe,iso_nuebar[tbl_number][nuclide]);
}


double 
weak_Q(double T9, double lg_rhoYe, int nuclide, int flavour)
{

  /* no muon/tau neutrinos at all! - just for consistency with thermal
      neutrino emission routines */

  if(flavour==1) return ffn_rate_nue(T9,lg_rhoYe,nuclide,0)*ffn_rate_nue(T9,lg_rhoYe,nuclide,1);
  if(flavour==2) return ffn_rate_nuebar(T9,lg_rhoYe,nuclide,0)*ffn_rate_nuebar(T9,lg_rhoYe,nuclide,1);

  if(flavour==3) return 0.0;
  if(flavour==4) return 0.0;

  return 0.0;
}

double 
weak_R(double T9, double lg_rhoYe, int nuclide, int flavour)
{

  /* no muon/tau neutrinos at all! - just for consistency with thermal
      neutrino emission routines */

  if(flavour==1) return ffn_rate_nue(T9,lg_rhoYe,nuclide,0);
  if(flavour==2) return ffn_rate_nuebar(T9,lg_rhoYe,nuclide,0);

  if(flavour==3) return 0.0;
  if(flavour==4) return 0.0;

  return 0.0;
}


//Normalized Fermi-Dirac neutrino spectrum
double 
fd_norm(double Enu, double Eavg)
{
  const double a=17.3574, b = 3.15137;
       
  return a*Enu*Enu/Eavg/Eavg/Eavg/(1.0 + exp(b*Enu/Eavg) );

}


/* for electron capture, change sign of mu for positron capture */
inline double 
capture(double Enu, double Q_eff, double kT, double mu)
{
  const double me=0.511;
  
  if(Enu<=0.0) return 0.0;
  if(Enu<=Q_eff+me) return 0.0;
  
  return log(2.0)*Enu*Enu*(Enu-Q_eff)*sqrt( (Enu-Q_eff-me)*(Enu-Q_eff+me) )/
                      (1.0+exp( (Enu-Q_eff-mu)/kT ) )/me/me/me/me/me;
}

/*NOTE: mu sign is for electrons, i.e. always positive */
double proton(double Enu, double kT, double mu, int flavour)
{

  const double me=0.511, Q=-1.29, half_life=1083.0; // ustalic ile obecnie wynosi ten wspolczynnik !!!!!
  double coeff;
  
  if(flavour>1) return 0.0;

  if(Enu<=0.0) return 0.0;
  
  coeff = log(2.0)/half_life/me/me/me/me/me;

  return coeff*Enu*Enu*(Enu-Q)*sqrt( (Enu-Q)*(Enu-Q) - me*me )/( 1.0 + exp( (Enu-Q-mu)/kT ) );

}

/* for beta minus decay, change sign of mu for beta plus decay */
inline double 
decay(double Enu, double Q_eff, double kT, double mu)
{
  const double me=0.511;

  if(Enu<=0.0) return 0.0;
  if(Enu>=Q_eff-me) return 0.0;

  return log(2.0)*Enu*Enu*(Q_eff-Enu)*sqrt( (Enu-Q_eff-me)*(Enu-Q_eff+me) )/
                      (1.0+exp( (Enu-Q_eff+mu)/kT ) )/me/me/me/me/me;
}

/*NOTE: mu sign is for electrons, i.e. always positive */
double 
neutron(double Enu, double kT, double mu, int flavour)
{

  const double me=0.511, Q=1.29, half_life=1083.0; 
  double coeff, spectrum=0.0;

  if(flavour!=2) return 0.0;

  coeff = log(2.0)/half_life/me/me/me/me/me;
  
  if( (Enu<=Q+me) && (Enu>=Q-me) ) return  0.0;
  if(Enu>Q+me)  spectrum = coeff*Enu*Enu*(Enu-Q)*sqrt( (Enu-Q-me)*(Enu-Q+me) )/( 1.0 + exp( (Enu-Q+mu)/kT ) );
  if(Enu<Q-me)  spectrum = coeff*Enu*Enu*(Q-Enu)*sqrt( (Enu-Q-me)*(Enu-Q+me) )/( 1.0 + exp( (Enu-Q+mu)/kT ) );
  

  return spectrum;

}


double 
neutron_R(double kT, double mu)
{
  
  
  double alpha,beta,eta,theta;
  double fd,fdd1,fdd2;
  // Fermi-Dirac integrals degree required
  double k_12=0.5,k_32=1.5,k_52=2.5,k_72=3.5;
  const double me=0.511, Q=1.29, half_life=1083.0;
  double coeff=1.0,rate_decay=0.0,rate_capture=0.0;
  /* 5-point Gaussian quadrature for Enu<Q-me */
  const double w[5]={0.000176799,0.00659261,0.00659261,0.0214928,0.00547821};
  const double x[5]={-1.25346,-1.11023,-0.9005,-0.690766 ,-0.547543};
  int ii;
  
  
  for(ii=0;ii<80;ii++)
    rate_decay = rate_decay + neutron(ii*0.01,kT,mu,2)*0.01;
  
  eta=-(mu+me)/kT;
  theta=kT/me;
  
  
  
  //dfermi_(&k_12,&eta,&theta,&fd,&fdd1,&fdd2);
  fd = Ffermi(k_12,eta,theta);
  
  rate_capture = rate_capture+(me*me*me+2.0*me*me*Q+me*Q*Q)*fd;
  
  //dfermi_(&k_32,&eta,&theta,&fd,&fdd1,&fdd2);
  fd = Ffermi(k_32,eta,theta);
  
  rate_capture = rate_capture+(3.0*kT*me*me+4.0*kT*me*Q+kT*Q*Q)*fd;
  
  //dfermi_(&k_52,&eta,&theta,&fd,&fdd1,&fdd2);
  fd = Ffermi(k_52,eta,theta);
  
  rate_capture = rate_capture+(3.0*kT*kT*me+2.0*kT*kT*Q)*fd;

  //dfermi_(&k_72,&eta,&theta,&fd,&fdd1,&fdd2);
  fd = Ffermi(k_72,eta,theta);
  
  rate_capture = rate_capture+kT*kT*kT*fd;
  
  rate_capture = rate_capture*sqrt(2.0*me)*pow(kT,1.5);
  
  coeff = log(2.0)/half_life/me/me/me/me/me;
  
  
  
  return rate_decay+rate_capture*coeff;
#if 0  
  const double me=0.511, Q=1.29, half_life=1083.0; 
  double coeff, rate1=0.0, rate2=0.0;
  /* 5-point Gaussian quadrature for Enu<Q-me */
  const double w1[5]={0.000176799,0.00659261,0.00659261,0.0214928,0.00547821};
  const double x1[5]={-1.25346,-1.11023,-0.9005,-0.690766 ,-0.547543};
  /* 5-point custom quadrature for Integrate[g[x]/(1+Exp[x/a]),{x,1,Infinity}]
  with (-1+Log[1+Exp[1/a]) scalling
  */
  const double w2[5]={0.4678602960432629,0.43309001130094066,
   0.0943535707356209,0.004665701305738638,
   0.000030420613986481256};
  const double x2[5]={0.29997378585182466,1.5225280900298397,
   3.7119291537368135,7.188390073437786,
   12.738063782307867};
  double y,a;
  int ii;
  
 
  coeff = log(2.0)/half_life/me/me/me/me/me;
  
  a=kT/(me+mu);
  
  for(ii=0;ii<5;ii++)
    rate1 = rate1 + w1[ii]/(1.0 + exp((x1[ii] + mu)/kT));
  
  for(ii=0;ii<5;ii++)
  {
    y=me+x2[ii]*kT;
    rate2 = rate2 + w2[ii]*(y+Q)*(y+Q)*y*sqrt(y*y-me*me);
  }
    rate2 = rate2*(me+mu)*(a*log(1.0+exp(1.0/a))-1.0);
    
  return (rate1+rate2)*coeff;
#endif /* #if 0 ends */  
}


double 
neutron_Q(double kT, double mu)
{
  
  
  double alpha,beta,eta,theta;
  double fd,fdd1,fdd2;
  // Fermi-Dirac integrals degree required
  double k_12=0.5,k_32=1.5,k_52=2.5,k_72=3.5,k_92=4.5;
  const double me=0.511, Q=1.29, half_life=1083.0;
  double coeff=1.0,rate_decay=0.0,rate_capture=0.0;
  int ii;
  
  
  for(ii=0;ii<80;ii++)
    rate_decay = rate_decay + (ii*0.01)*neutron(ii*0.01,kT,mu,2)*0.01;
  
  eta=-(mu+me)/kT;
  theta=kT/me;
  
  //dfermi_(&k_12,&eta,&theta,&fd,&fdd1,&fdd2);
  fd = Ffermi(k_12,eta,theta);
  
  rate_capture = rate_capture+(me*me*me*me+3.0*me*me*me*Q+3.0*me*me*Q*Q+me*Q*Q*Q)*fd;
  
  //dfermi_(&k_32,&eta,&theta,&fd,&fdd1,&fdd2);
  fd = Ffermi(k_32,eta,theta);
  
  rate_capture = rate_capture+(4.0*kT*me*me*me+9.0*kT*me*me*Q+6.0*kT*me*Q*Q+kT*Q*Q*Q)*fd;
  
  //dfermi_(&k_52,&eta,&theta,&fd,&fdd1,&fdd2);
  fd = Ffermi(k_52,eta,theta);
  
  rate_capture = rate_capture+(6.0*kT*kT*me*me+9.0*kT*kT*me*Q+3.0*kT*kT*Q*Q)*fd;
  
  //dfermi_(&k_72,&eta,&theta,&fd,&fdd1,&fdd2);
  fd = Ffermi(k_72,eta,theta);
  
  rate_capture = rate_capture+(4.0*kT*kT*kT*me+3.0*kT*kT*kT*Q)*fd;
  
  //dfermi_(&k_92,&eta,&theta,&fd,&fdd1,&fdd2);
  fd = Ffermi(k_92,eta,theta);
  
  rate_capture = rate_capture+kT*kT*kT*kT*fd;
  
  rate_capture = rate_capture*sqrt(2.0*me)*pow(kT,1.5);
  
  coeff = log(2.0)/half_life/me/me/me/me/me;
  
  
  
  return (rate_decay+rate_capture*coeff);
}  

double 
proton_R(double kT, double mu)
{
  
  
  double alpha,beta,eta,theta;
  double fd,fdd1,fdd2;
  // Fermi-Dirac integrals degree required
  double k_12=0.5,k_32=1.5,k_52=2.5,k_72=3.5;
  const double me=0.511, Q=1.29, half_life=1083.0;
  double coeff=1.0,rate_decay=0.0,rate_capture=0.0;
  /* 5-point Gaussian quadrature for Enu<Q-me */
  const double w[5]={0.000176799,0.00659261,0.00659261,0.0214928,0.00547821};
  const double x[5]={-1.25346,-1.11023,-0.9005,-0.690766 ,-0.547543};
  int ii;
  double Enu;
  
  for(ii=0;ii<779;ii++)
  {  
    Enu=ii*0.001;
//    rate_decay = rate_decay + exp((Enu-Q+mu)/kT)*neutron(Enu,kT,mu,2)*0.01;
    rate_decay = rate_decay + 0.001*Enu*Enu*(Q-Enu)*sqrt( (Enu-Q-me)*(Enu-Q+me) )/( 1.0 + exp( (-Enu+Q-mu)/kT ) );
  }
  
  
  eta=(mu-me)/kT;
  theta=kT/me;
  
  
  
  //dfermi_(&k_12,&eta,&theta,&fd,&fdd1,&fdd2);
  fd = Ffermi(k_12,eta,theta);
  
  rate_capture = rate_capture+(me*me*me-2.0*me*me*Q+me*Q*Q)*fd;
  
  //dfermi_(&k_32,&eta,&theta,&fd,&fdd1,&fdd2);
  fd = Ffermi(k_32,eta,theta);
  
  rate_capture = rate_capture+(3.0*kT*me*me-4.0*kT*me*Q+kT*Q*Q)*fd;
  
  //dfermi_(&k_52,&eta,&theta,&fd,&fdd1,&fdd2);
  fd = Ffermi(k_52,eta,theta);
  
  rate_capture = rate_capture+(3.0*kT*kT*me-2.0*kT*kT*Q)*fd;
  
  //dfermi_(&k_72,&eta,&theta,&fd,&fdd1,&fdd2);
  fd = Ffermi(k_72,eta,theta);
  
  rate_capture = rate_capture+kT*kT*kT*fd;
  
  rate_capture = rate_capture*sqrt(2.0*me)*pow(kT,1.5);
  
  coeff = log(2.0)/half_life/me/me/me/me/me;
  
  
  
  return (-rate_decay+rate_capture)*coeff;
}  

double 
proton_Q(double kT, double mu)
{
  
  
  double alpha,beta,eta,theta;
  double fd,fdd1,fdd2;
  // Fermi-Dirac integrals degree required
  double k_12=0.5,k_32=1.5,k_52=2.5,k_72=3.5,k_92=4.5;
  const double me=0.511, Q=1.29, half_life=1083.0;
  double coeff=1.0,rate_decay=0.0,rate_capture=0.0;
  int ii;
  double Enu;
  
  for(ii=0;ii<779;ii++)
  {  
    Enu=ii*0.001;
    //    rate_decay = rate_decay + exp((Enu-Q+mu)/kT)*neutron(Enu,kT,mu,2)*0.01;
    rate_decay = rate_decay + 0.001*Enu*Enu*Enu*(Q-Enu)*sqrt( (Enu-Q-me)*(Enu-Q+me) )/( 1.0 + exp( (-Enu+Q-mu)/kT ) );
  }
  
  
  eta=(mu-me)/kT;
  theta=kT/me;
  
  //dfermi_(&k_12,&eta,&theta,&fd,&fdd1,&fdd2);
  fd = Ffermi(k_12,eta,theta);
  
  rate_capture = rate_capture+(me*me*me*me - 3.0*me*me*me*Q+3.0*me*me*Q*Q-me*Q*Q*Q)*fd;
  
  //dfermi_(&k_32,&eta,&theta,&fd,&fdd1,&fdd2);
  fd = Ffermi(k_32,eta,theta);
  
  rate_capture = rate_capture+(4.0*kT*me*me*me-9.0*kT*me*me*Q+6.0*kT*me*Q*Q-kT*Q*Q*Q)*fd;
  
  //dfermi_(&k_52,&eta,&theta,&fd,&fdd1,&fdd2);
  fd = Ffermi(k_52,eta,theta);
  
  rate_capture = rate_capture+(6.0*kT*kT*me*me-9.0*kT*kT*me*Q+3.0*kT*kT*Q*Q)*fd;
  
  //dfermi_(&k_72,&eta,&theta,&fd,&fdd1,&fdd2);
  fd = Ffermi(k_72,eta,theta);
  
  rate_capture = rate_capture+(4.0*kT*kT*kT*me-3.0*kT*kT*kT*Q)*fd;
  
  //dfermi_(&k_92,&eta,&theta,&fd,&fdd1,&fdd2);
  fd = Ffermi(k_92,eta,theta);
  
  rate_capture = rate_capture+kT*kT*kT*kT*fd;
  
  rate_capture = rate_capture*sqrt(2.0*me)*pow(kT,1.5);
  
  coeff = log(2.0)/half_life/me/me/me/me/me;
  
  
  
  return (-rate_decay+rate_capture)*coeff;
}  


double 
Enu_avg_capture(double Q_eff, double kT, double mu)
{
  double int_q=0.0,int_r=0.0;
  double Enu_max=100.0; //In principle infinite range, many nuclei has Enu_MIN>20 MeV , mu>25 for high densities
  double Enu_min=0.0,dEnu,Enu=0.01;
  int ii,Num=1000;
  
  if((Q_eff+0.511)>0.0) Enu_min=Q_eff+0.511; else Enu_min=0.0;
  Enu_max=10.0*kT+2.0*fabs(mu)+2.0*fabs(Q_eff);

  dEnu=(Enu_max-Enu_min)/Num;
  if(dEnu<=0.0) return 0.0;
  
  for(ii=1;ii<Num;ii++){
    Enu=Enu_min + ii*dEnu;
    int_r = int_r + dEnu    *capture(Enu,Q_eff,kT,mu ) ;
    int_q = int_q + dEnu*Enu*capture(Enu,Q_eff,kT,mu ) ;
  }

  Enu=Enu_max;
  int_r = int_r + 0.5*dEnu    *capture(Enu,Q_eff,kT,mu );
  int_q = int_q + 0.5*dEnu*Enu*capture(Enu,Q_eff,kT,mu );

  //Analytical formula for tail

  int_r = int_r + log(2.0)/pow(0.511,5.0)*exp((Q_eff+mu-Enu_max)/kT)*kT*
               ( pow(Enu_max,4.0) + 4.0*pow(Enu_max,3.0)*kT + 12.0*Enu_max*Enu_max*kT*kT + 24.0*Enu_max*kT*kT*kT + 
                  24.0*kT*kT*kT*kT);

  int_q = int_q + log(2.0)/pow(0.511,5.0)*exp((Q_eff+mu-Enu_max)/kT)*kT*
               ( pow(Enu_max,5.0) + 5.0*pow(Enu_max,4.0)*kT + 20.0*Enu_max*Enu_max*Enu_max*kT*kT 
                 + 60.0*Enu_max*Enu_max*kT*kT*kT + 120.0*Enu_max*kT*kT*kT*kT+120.0*pow(kT,5.0));

  return int_q/int_r;
}

double 
Enu_avg_decay(double Q_eff, double kT, double mu)
{
  const double me=0.511;
  double int_q=0.0,int_r=0.0;
  double Enu_min, Enu_max, dEnu,Enu;
  int ii,Num=1024;
  
  Enu_min=0.0;
  Enu_max=Q_eff-me;

  dEnu=(Enu_max-Enu_min)/Num;

  if(dEnu<0.0) return 0.0;
  
  for(ii=1;ii<Num;ii++){
    Enu=Enu_min + ii*dEnu;
    int_r = int_r + dEnu*decay(Enu,Q_eff,kT,mu ) ;
    int_q = int_q + dEnu*Enu*decay(Enu,Q_eff,kT,mu ) ;
  }

  return int_q/int_r;


}

double 
Enu_avg_capture_and_decay(double Q_eff, double kT, double mu)
{

  int ii;
  double Enu;
  const double me=0.511;

/* decay */
  double int_q_decay=0.0,int_r_decay=0.0;
  double Enu_min_decay, Enu_max_decay, dEnu_decay;
  int Num_decay=256;
/* capture */
  double int_q_capture=0.0, int_r_capture=0.0;
  double Enu_max_capture=50.0; //In principle infinite range, many nuclei has Enu_MIN>20 MeV , mu>25 for high densities
  double Enu_min_capture=0.0, dEnu_capture;
  int Num_capture=256;


/* capture integral start */
  
  if((Q_eff+me)>0.0) Enu_min_capture=Q_eff+me; else Enu_min_capture=0.0;

  dEnu_capture=(Enu_max_capture-Enu_min_capture)/Num_capture;

  
  for(ii=1;ii<Num_capture;ii++){
    Enu=Enu_min_capture + ii*dEnu_capture;
    int_r_capture = int_r_capture + dEnu_capture*capture(Enu,Q_eff,kT,mu ) ;
    int_q_capture = int_q_capture + dEnu_capture*Enu*capture(Enu,Q_eff,kT,mu ) ;
  }

  Enu=Enu_max_capture;
  int_r_capture = int_r_capture + 0.5*dEnu_capture*capture(Enu,Q_eff,kT,mu );
  int_q_capture = int_q_capture + 0.5*dEnu_capture*Enu*capture(Enu,Q_eff,kT,mu );

/* capture integral done */

  mu = 0.0-mu;  //if nuclei capture electrons, positrons are produced in reverse reaction, AND vice-versa

/* decay integral start */
  Enu_min_decay=0.0;
  Enu_max_decay=Q_eff-me;

  dEnu_decay=(Enu_max_decay-Enu_min_decay)/Num_decay;

  for(ii=1;ii<Num_decay;ii++){
    Enu=Enu_min_decay + ii*dEnu_decay;
    int_r_decay = int_r_decay + dEnu_decay*decay(Enu,Q_eff,kT,mu ) ;
    int_q_decay = int_q_decay + dEnu_decay*Enu*decay(Enu,Q_eff,kT,mu ) ;
  }

  /* decay intergral done */

  
  return (int_q_capture + int_q_decay)/(int_r_capture + int_r_decay);

}


/* Function calculates effective Q-values for electron captures
from given temperature kT, chemical potential mu and average neutrino energy
Enu_avg (all in MeV). This could give good fit for the neutrino
spectra as explained by Langanke (), except for multi peaked cases */
double 
Q_eff_capture(double Enu_avg, double  kT, double mu, double EPSABS)
{
  /* i hope no nuclei with Q-value larger than 30 MeV exist */
  double Q_eff=0.75, Q_eff_left=-30.0, Q_eff_right=30.0;
  
  while((Q_eff_right-Q_eff_left)>EPSABS  )
  { 

    if(Enu_avg_capture(Q_eff,kT,mu)>Enu_avg) 
      Q_eff_right=Q_eff; 
    else Q_eff_left=Q_eff;

    Q_eff=0.5*Q_eff_left+0.5*Q_eff_right;
  }
      
  return Q_eff;
}


double 
Q_eff_decay(double Enu_avg, double  kT, double mu, double EPSABS)
{
  /* i hope no nuclei with Q-value larger than 30 MeV exist */
  double Q_eff=3.0, Q_eff_left=0.511, Q_eff_right=30.0;
  
  while((Q_eff_right-Q_eff_left)>EPSABS  )
  { 
    if(Enu_avg_decay(Q_eff,kT,mu)>Enu_avg) Q_eff_right=Q_eff; 
        else Q_eff_left=Q_eff;
    Q_eff=0.5*Q_eff_left+0.5*Q_eff_right;
  }
      
  return Q_eff;
}


/* WARNING! this is non-monotonic function, solution is NOT UNIQUE! EXPERIMENTAL !*/
double 
Q_eff_capture_and_decay(double Enu_avg, double  kT, double mu,  double Q_guess, double EPSABS)
{
  /* i hope no nuclei with Q-value larger than 25 MeV exist */
  double Q_eff=3.0*kT, Q_eff_left=0.511, Q_eff_right=30.0;

  Q_eff_left= Q_guess/2.0;
  Q_eff_right =   Q_guess*2.0;
  Q_eff = Q_guess;
  
  while((Q_eff_right-Q_eff_left)>EPSABS  )
  { 
    if(Enu_avg_capture_and_decay(Q_eff,kT,mu)>Enu_avg) Q_eff_right=Q_eff; 
        else Q_eff_left=Q_eff;
    Q_eff=0.5*Q_eff_left+0.5*Q_eff_right;
  }
      
  return Q_eff;


}

double 
effective_capture_rate(double Q_eff, double kT, double mu)
{
  double int_r=0.0;
  double Enu_min=0.0,Enu_max=100.0, dEnu=0.01, Enu;
  int ii,Num=1000;
  
/* Trapezoidal integration rule up to Enu_max*/
  Enu_min=0.0;
  Enu_max=10.0*kT+2.0*mu+2.0*fabs(Q_eff);
  if((Q_eff+0.511)>0.0) Enu_min=Q_eff+0.511;
  
  dEnu=(Enu_max-Enu_min)/Num;
  
  for(ii=1;ii<Num;ii++){
    Enu=Enu_min + ii*dEnu;
    int_r = int_r + dEnu*capture(Enu,Q_eff,kT,mu );
  }

  Enu=Enu_max;
  int_r = int_r + 0.5*dEnu*capture(Enu,Q_eff,kT,mu );

/* Analytical formula for integral of the tail */

  int_r = int_r + log(2.0)/pow(0.511,5.0)*exp((Q_eff+mu-Enu_max)/kT)*kT*
               ( pow(Enu_max,4.0) + 4.0*pow(Enu_max,3.0)*kT + 12.0*Enu_max*Enu_max*kT*kT + 24.0*Enu_max*kT*kT*kT + 
                  24.0*kT*kT*kT*kT);


  return int_r;

}

double 
effective_decay_rate(double Q_eff, double kT, double mu)
{
  double int_r=0.0;
  double Enu_min, Enu_max, dEnu,Enu;
  int ii,Num=1024;
  
  Enu_min =  0.0;
  Enu_max =  Q_eff-0.511;
  
  dEnu=(Enu_max-Enu_min)/Num;
  
  for(ii=1;ii<Num;ii++){
    Enu=Enu_min + ii*dEnu;
    int_r = int_r + dEnu*decay(Enu,Q_eff,kT,mu );
  }

  return int_r;

}


/*
  NOTE: FUNCTION is not intended to be used in large-scale processing! Use 
  NSE_neutrino_spectrum_table() below. Useful only for testing/debugging 
  and producing single (rho, kT, Ye) triad graphs 
*/ 
double 
NSE_neutrino_spectrum(double Enu, double kT, double rho, double Ye, int flavour)
{
  double   spectrum=0.0;
  int inuc=0; 
  double T9, lg_rhoYe, rhoYe, lg_rho, lg_T, mu=0.0;
  double Q_weak_nue, R_weak_nue,    weak_rate_nue, avg_neutrino_energy_weak_nue , X_NSE_nue;
  double lambda_nue, lambda_nuebar;
  double Q_weak_nuebar, R_weak_nuebar, weak_rate_nuebar, avg_neutrino_energy_weak_nuebar, X_NSE_nuebar;

  const  double proton_mass = 1.67262158E-24, qe=1.6021892E-6;

  double eff_rate, Q_value_eff, capture_rate, decay_rate; 

  if(flavour>2) return 0.0;
  
  inuc = 0; //protons and neutrons are handled without use of FFN tables here

  lg_rho=log10(rho);
  lg_rhoYe = log10(rho*Ye);
  rhoYe=rho*Ye;
  T9=kT*11.6045;
  lg_T = log10(T9*1E9);
  mu = chemical_potential(kT, rhoYe, 1E-12);  

  if(flavour==1){ 
    X_NSE_nue    =  NSE(T9,rho,Ye,1,0);
    spectrum = spectrum +   proton(Enu,kT, mu, flavour)*X_NSE_nue*rho/proton_mass;
  }


  if(flavour==2){ 
    X_NSE_nuebar    =  NSE(T9,rho,Ye,0,1);
    spectrum= spectrum +   neutron(Enu,kT, mu, flavour)*X_NSE_nuebar*rho/proton_mass;
  }

  for(inuc=1;inuc<189;inuc++){
    
    if(flavour==1){ 
      X_NSE_nue    =  NSE(T9,rho,Ye,Z_nue[inuc],A_nue[inuc]-Z_nue[inuc]);
      avg_neutrino_energy_weak_nue = ffn_rate_nue(T9,lg_rhoYe,inuc,1);
      weak_rate_nue                = ffn_rate_nue(T9,lg_rhoYe,inuc,0);
      
      lambda_nue=weak_rate_nue*X_NSE_nue/A_nue[inuc];
      R_weak_nue=weak_rate_nue*X_NSE_nue*rho/proton_mass/A_nue[inuc];
      Q_weak_nue=R_weak_nue*avg_neutrino_energy_weak_nue*qe;
      capture_rate             = ffn_rate_nue(T9,lg_rhoYe,inuc,3);
      decay_rate               = ffn_rate_nue(T9,lg_rhoYe,inuc,2);

      if( (capture_rate>=decay_rate) && (avg_neutrino_energy_weak_nue>3.0*kT) ) 
        Q_value_eff = Q_eff_capture(avg_neutrino_energy_weak_nue,kT,mu,1E-12);
      else  
        Q_value_eff = Q_eff_decay(avg_neutrino_energy_weak_nue,kT,-mu,1E-12);
    
      if( (capture_rate>=decay_rate) && (avg_neutrino_energy_weak_nue>3.0*kT)) 
        eff_rate = effective_capture_rate(Q_value_eff,kT,mu)/weak_rate_nue;
      else
        eff_rate =   effective_decay_rate(Q_value_eff,kT,-mu)/weak_rate_nue;

      if((capture_rate>=decay_rate) && (avg_neutrino_energy_weak_nue>3.0*kT) )  
          spectrum = spectrum +   capture(Enu,Q_value_eff, kT,  mu)/eff_rate*X_NSE_nue*rho/proton_mass/A_nue[inuc];
      else
          spectrum = spectrum +     decay(Enu,Q_value_eff, kT, -mu)/eff_rate*X_NSE_nue*rho/proton_mass/A_nue[inuc];

    }

    
    if(flavour==2){ /* antineutrinos */
      X_NSE_nuebar    =  NSE(T9,rho,Ye,Z_nuebar[inuc],A_nuebar[inuc]-Z_nuebar[inuc]);
      avg_neutrino_energy_weak_nuebar = ffn_rate_nuebar(T9,lg_rhoYe,inuc,1);
      weak_rate_nuebar                = ffn_rate_nuebar(T9,lg_rhoYe,inuc,0);
      lambda_nuebar=weak_rate_nuebar*X_NSE_nuebar/A_nuebar[inuc];
      R_weak_nuebar=weak_rate_nuebar*X_NSE_nuebar*rho/proton_mass/A_nuebar[inuc];
      Q_weak_nuebar=R_weak_nuebar*avg_neutrino_energy_weak_nuebar*qe;
      capture_rate             = ffn_rate_nuebar(T9,lg_rhoYe,inuc,3);
      decay_rate               = ffn_rate_nuebar(T9,lg_rhoYe,inuc,2);

      
      if( (capture_rate>=decay_rate) && (avg_neutrino_energy_weak_nuebar>3.0*kT) ) 
        Q_value_eff = Q_eff_capture(avg_neutrino_energy_weak_nuebar,kT,-mu,1E-12);
      else
        Q_value_eff = Q_eff_decay(avg_neutrino_energy_weak_nuebar,kT,mu,1E-12);

     if( (capture_rate>=decay_rate) && (avg_neutrino_energy_weak_nuebar>3.0*kT) ) 
        eff_rate = effective_capture_rate(Q_value_eff,kT,-mu)/weak_rate_nuebar;
      else   
        eff_rate = effective_decay_rate(Q_value_eff,kT,mu)/weak_rate_nuebar;

      if( (capture_rate>=decay_rate) && (avg_neutrino_energy_weak_nuebar>3.0*kT) )  
          spectrum = spectrum + capture(Enu,Q_value_eff, kT, -mu)/eff_rate*X_NSE_nuebar*rho/proton_mass/A_nuebar[inuc];
      else
          spectrum = spectrum +   decay(Enu,Q_value_eff, kT,  mu)/eff_rate*X_NSE_nuebar*rho/proton_mass/A_nuebar[inuc];
    }
  }

  return spectrum;
}


/* 
  NOTE: FUNCTION is not intended to be used in large-scale processing! Use 
  NSE_neutrino_spectrum_table() below. Useful only for testing/debugging 
  and producing single (rho, kT, Ye) triad graphs 
*/ 
double 
NSE_neutrino_spectrum_tbl(double Enu, double kT, double rho, double Ye, int flavour)
{
  double   spectrum=0.0;
  int inuc=0; 
  double T9, lg_rhoYe, rhoYe, lg_rho, lg_T, mu=0.0;
  double Q_weak_nue, R_weak_nue,    weak_rate_nue, avg_neutrino_energy_weak_nue , X_NSE_nue;
  double lambda_nue, lambda_nuebar;
  double Q_weak_nuebar, R_weak_nuebar, weak_rate_nuebar, avg_neutrino_energy_weak_nuebar, X_NSE_nuebar;

  const  double proton_mass = 1.67262158E-24, qe=1.6021892E-6;

  double eff_rate, Q_value_eff, capture_rate, decay_rate; 

  if(flavour>2) return 0.0;

  inuc = 0; //protons and neutrons are handled without use of FFN tables here

  lg_rho=log10(rho);
  lg_rhoYe = log10(rho*Ye);
  rhoYe=rho*Ye;
  T9=kT*11.6045;
  lg_T = log10(T9*1E9);
  mu = chemical_potential(kT, rhoYe, 1E-12);  

  if(flavour==1){ 
    X_NSE_nue    =  NSE(T9,rho,Ye,1,0);
    spectrum = spectrum +   proton(Enu,kT, mu, flavour)*X_NSE_nue*rho/proton_mass;
  }
 

  if(flavour==2){ 
    X_NSE_nuebar    =  NSE(T9,rho,Ye,0,1);
    spectrum= spectrum +   neutron(Enu,kT, mu, flavour)*X_NSE_nuebar*rho/proton_mass;
  }


  for(inuc=1;inuc<189;inuc++){
    
    if(flavour==1){ 
      X_NSE_nue    =  NSE(T9,rho,Ye,Z_nue[inuc],A_nue[inuc]-Z_nue[inuc]);
      avg_neutrino_energy_weak_nue = ffn_rate_nue(T9,lg_rhoYe,inuc,1);
      weak_rate_nue                = ffn_rate_nue(T9,lg_rhoYe,inuc,0);
      lambda_nue=weak_rate_nue*X_NSE_nue/A_nue[inuc];
      R_weak_nue=weak_rate_nue*X_NSE_nue*rho/proton_mass/A_nue[inuc];
      Q_weak_nue=R_weak_nue*avg_neutrino_energy_weak_nue*qe;
      capture_rate             = ffn_rate_nue(T9,lg_rhoYe,inuc,3);
      decay_rate               = ffn_rate_nue(T9,lg_rhoYe,inuc,2);

      Q_value_eff = ffn_rate_nue(T9,lg_rhoYe,inuc,5);
      eff_rate    = ffn_rate_nue(T9,lg_rhoYe,inuc,4);

      if( (capture_rate>=decay_rate)  && (avg_neutrino_energy_weak_nue>3.0*kT) )  
          spectrum = spectrum + capture(Enu,Q_value_eff, kT,  mu)/eff_rate*X_NSE_nue*rho/proton_mass/A_nue[inuc];
      else
          spectrum = spectrum +   decay(Enu,Q_value_eff, kT, -mu)/eff_rate*X_NSE_nue*rho/proton_mass/A_nue[inuc];

    }

    if(flavour==2){ /* antineutrinos */
      X_NSE_nuebar    =  NSE(T9,rho,Ye,Z_nuebar[inuc],A_nuebar[inuc]-Z_nuebar[inuc]);
      avg_neutrino_energy_weak_nuebar = ffn_rate_nuebar(T9,lg_rhoYe,inuc,1);
      weak_rate_nuebar                = ffn_rate_nuebar(T9,lg_rhoYe,inuc,0);
      lambda_nuebar=weak_rate_nuebar*X_NSE_nuebar/A_nuebar[inuc];
      R_weak_nuebar=weak_rate_nuebar*X_NSE_nuebar*rho/proton_mass/A_nuebar[inuc];
      Q_weak_nuebar=R_weak_nuebar*avg_neutrino_energy_weak_nuebar*qe;
      capture_rate             = ffn_rate_nuebar(T9,lg_rhoYe,inuc,3);
      decay_rate               = ffn_rate_nuebar(T9,lg_rhoYe,inuc,2);

      Q_value_eff = ffn_rate_nuebar(T9,lg_rhoYe,inuc,5);     
      eff_rate    = ffn_rate_nuebar(T9,lg_rhoYe,inuc,4);

      if((capture_rate>=decay_rate) && (avg_neutrino_energy_weak_nuebar>3.0*kT))  
          spectrum = spectrum + capture(Enu,Q_value_eff, kT, -mu)/eff_rate*X_NSE_nuebar*rho/proton_mass/A_nuebar[inuc];
      else
          spectrum = spectrum +   decay(Enu,Q_value_eff, kT,  mu)/eff_rate*X_NSE_nuebar*rho/proton_mass/A_nuebar[inuc];
    }
  }

  return spectrum;
}



/*
Function fill table given as a first argument with neutrino spectrum
determined by triad kT, rho, Ye (2nd,3th,4th argument) with flavour
as 5th arg, and parametrs of the spectrum as remaining args:
Enu_min, Enu_max, Npoints, and logscale: 0 for linear, 1 for logscale

TODO: pass as an argument pointer to table filled with
wanted neutrino energies, rather than Emin, Emax and 
spacing logscale (?, maybe linear and logscale is enough ?) 

*/
void 
NSE_neutrino_spectrum_table(double *spectrum, double kT, double mu, double rho, double Ye, int flavour, double Enu_min, double Enu_max, int num, int logscale)
{

  int inuc=0, ii; 
  double T9, lg_rhoYe, rhoYe, lg_rho, lg_T;
  double Q_weak_nue, R_weak_nue,    weak_rate_nue, avg_neutrino_energy_weak_nue , X_NSE_nue;
  double lambda_nue, lambda_nuebar;
  double Q_weak_nuebar, R_weak_nuebar, weak_rate_nuebar, avg_neutrino_energy_weak_nuebar, X_NSE_nuebar;

  const  double proton_mass = 1.67262158E-24, qe=1.6021892E-6;

  double eff_rate, Q_value_eff, capture_rate, decay_rate; 
  double Enu, dEnu=(Enu_max-Enu_min)/num, lg_Enu_min, lg_Enu_max;
  
  lg_Enu_min = log10(Enu_min); 
  lg_Enu_max = log10(Enu_max);

  for(ii=0;ii<=num;ii++) spectrum[ii]=0.0;

  NSE_nuclei_neutrino_spectrum_table  (spectrum, kT, mu, rho, Ye, flavour, Enu_min, Enu_max, num, logscale); 
  NSE_nucleons_neutrino_spectrum_table(spectrum, kT, mu, rho, Ye, flavour, Enu_min, Enu_max, num, logscale); 
  

}

void 
NSE_nucleons_neutrino_spectrum_table(double *spectrum, double kT, double mu, double rho, double Ye, int flavour, double Enu_min, double Enu_max, int num, int logscale)
{
  
  int ii; 
  double T9, lg_rhoYe, rhoYe, lg_rho, lg_T;
  double Q_weak_nue, R_weak_nue,    weak_rate_nue, avg_neutrino_energy_weak_nue , X_NSE_nue;
  double lambda_nue, lambda_nuebar;
  double Q_weak_nuebar, R_weak_nuebar, weak_rate_nuebar, avg_neutrino_energy_weak_nuebar, X_NSE_nuebar;
  
  const  double proton_mass = 1.67262158E-24, qe=1.6021892E-6;
  
  double eff_rate, Q_value_eff, capture_rate, decay_rate; 
  double Enu, dEnu=(Enu_max-Enu_min)/num, lg_Enu_min, lg_Enu_max;
  
  lg_Enu_min = log10(Enu_min); 
  lg_Enu_max = log10(Enu_max);
  
/*
NSE_nucleons_* and NSE_nuclei_* are NOT intended to be called
directly - user must initialize table spectrum[] !!!
*/
//  for(ii=0;ii<=num;ii++) spectrum[ii]=0.0; 
  
  

  
  lg_rho=log10(rho);
  lg_rhoYe = log10(rho*Ye);
  rhoYe=rho*Ye;
  T9=kT*11.6045;
  lg_T = log10(T9*1E9);
  //  mu = chemical_potential(kT, rhoYe, 1E-12);  
  //mu = chemical_potential_tbl(lg_T, lg_rhoYe);  
  
  if(flavour==1){ 
    X_NSE_nue    =  NSE(T9,rho,Ye,1,0); //Z=1, N=0 for protons
    //    #pragma omp parallel for
    for(ii=0;ii<=num;ii++){ 
      if(!logscale) Enu=Enu_min+ii*dEnu; else Enu=pow(10.0, lg_Enu_min + (lg_Enu_max-lg_Enu_min)*ii/num );
      spectrum[ii] = spectrum[ii] +  proton(Enu,kT, mu, flavour)*X_NSE_nue*rho/proton_mass;
    }
  }
  
  if(flavour==2){ 
    X_NSE_nuebar    =  NSE(T9,rho,Ye,0,1); //Z=0,N=1 for neutrons
    //    #pragma omp parallel for    
    for(ii=0;ii<=num;ii++){ 
      if(!logscale) Enu=Enu_min+ii*dEnu; else Enu=pow(10.0, lg_Enu_min + (lg_Enu_max-lg_Enu_min)*ii/num );
      spectrum[ii] = spectrum[ii] +  neutron(Enu,kT, mu, flavour)*X_NSE_nuebar*rho/proton_mass;
    }
  }
  
  
}

void 
NSE_nuclei_neutrino_spectrum_table(double *spectrum, double kT, double mu, double rho, double Ye, int flavour, double Enu_min, double Enu_max, int num, int logscale)
{
  
  int inuc=0, ii; 
  double T9, lg_rhoYe, rhoYe, lg_rho, lg_T;
  double Q_weak_nue, R_weak_nue,    weak_rate_nue, avg_neutrino_energy_weak_nue , X_NSE_nue;
  double lambda_nue, lambda_nuebar;
  double Q_weak_nuebar, R_weak_nuebar, weak_rate_nuebar, avg_neutrino_energy_weak_nuebar, X_NSE_nuebar;
  
  const  double proton_mass = 1.67262158E-24, qe=1.6021892E-6;
  
  double eff_rate, Q_value_eff, capture_rate, decay_rate; 
  double Enu, dEnu=(Enu_max-Enu_min)/num, lg_Enu_min, lg_Enu_max;
  
  lg_Enu_min = log10(Enu_min); 
  lg_Enu_max = log10(Enu_max);
  
  lg_rho=log10(rho);
  lg_rhoYe = log10(rho*Ye);
  rhoYe=rho*Ye;
  T9=kT*11.6045;
  lg_T = log10(T9*1E9);
  //  mu = chemical_potential(kT, rhoYe, 1E-12);  
  //mu = chemical_potential_tbl(lg_T, lg_rhoYe);  
  
  for(inuc=1;inuc<189;inuc++){
    
    if(flavour==1){ 
      X_NSE_nue    =  NSE(T9,rho,Ye,Z_nue[inuc],A_nue[inuc]-Z_nue[inuc]);
      avg_neutrino_energy_weak_nue = ffn_rate_nue(T9,lg_rhoYe,inuc,1);
      weak_rate_nue                = ffn_rate_nue(T9,lg_rhoYe,inuc,0);
      lambda_nue=weak_rate_nue*X_NSE_nue/A_nue[inuc];
      R_weak_nue=weak_rate_nue*X_NSE_nue*rho/proton_mass/A_nue[inuc];
      Q_weak_nue=R_weak_nue*avg_neutrino_energy_weak_nue*qe;
      capture_rate             = ffn_rate_nue(T9,lg_rhoYe,inuc,3);
      decay_rate               = ffn_rate_nue(T9,lg_rhoYe,inuc,2);
      
      Q_value_eff = ffn_rate_nue(T9,lg_rhoYe,inuc,5);
      eff_rate    = ffn_rate_nue(T9,lg_rhoYe,inuc,4);
      
      if( (capture_rate>=decay_rate)  && (avg_neutrino_energy_weak_nue>3.0*kT) )  
        for(ii=0;ii<=num;ii++){ 
          if(!logscale) Enu=Enu_min+ii*dEnu; else Enu=pow(10.0, lg_Enu_min + (lg_Enu_max-lg_Enu_min)*ii/num );
          spectrum[ii] = spectrum[ii] +  capture(Enu,Q_value_eff, kT,  mu)/eff_rate*X_NSE_nue*rho/proton_mass/A_nue[inuc];
        }
        else
        for(ii=0;ii<=num;ii++){ 
          if(!logscale) Enu=Enu_min+ii*dEnu; else Enu=pow(10.0, lg_Enu_min + (lg_Enu_max-lg_Enu_min)*ii/num );
          spectrum[ii] = spectrum[ii] +  decay(Enu,Q_value_eff, kT,  -mu)/eff_rate*X_NSE_nue*rho/proton_mass/A_nue[inuc];
        }  
          
    }
    
    if(flavour==2){ /* antineutrinos */
      X_NSE_nuebar    =  NSE(T9,rho,Ye,Z_nuebar[inuc],A_nuebar[inuc]-Z_nuebar[inuc]);
      avg_neutrino_energy_weak_nuebar = ffn_rate_nuebar(T9,lg_rhoYe,inuc,1);
      weak_rate_nuebar                = ffn_rate_nuebar(T9,lg_rhoYe,inuc,0);
      lambda_nuebar=weak_rate_nuebar*X_NSE_nuebar/A_nuebar[inuc];
      R_weak_nuebar=weak_rate_nuebar*X_NSE_nuebar*rho/proton_mass/A_nuebar[inuc];
      Q_weak_nuebar=R_weak_nuebar*avg_neutrino_energy_weak_nuebar*qe;
      capture_rate             = ffn_rate_nuebar(T9,lg_rhoYe,inuc,3);
      decay_rate               = ffn_rate_nuebar(T9,lg_rhoYe,inuc,2);
      
      eff_rate = ffn_rate_nuebar(T9,lg_rhoYe,inuc,4);
      Q_value_eff = ffn_rate_nuebar(T9,lg_rhoYe,inuc,5);     
      
      if( (capture_rate>=decay_rate)  && (avg_neutrino_energy_weak_nuebar>3.0*kT) )  
        for(ii=0;ii<=num;ii++){ 
          if(!logscale) Enu=Enu_min+ii*dEnu; else Enu=pow(10.0, lg_Enu_min + (lg_Enu_max-lg_Enu_min)*ii/num );
          spectrum[ii] = spectrum[ii] +  capture(Enu,Q_value_eff, kT,  -mu)/eff_rate*X_NSE_nuebar*rho/proton_mass/A_nuebar[inuc];
        }
        else
          for(ii=0;ii<=num;ii++){ 
            if(!logscale) Enu=Enu_min+ii*dEnu; else Enu=pow(10.0, lg_Enu_min + (lg_Enu_max-lg_Enu_min)*ii/num );
            spectrum[ii] = spectrum[ii] +    decay(Enu,Q_value_eff, kT,  mu)/eff_rate*X_NSE_nuebar*rho/proton_mass/A_nuebar[inuc];
          }  
    }
  }
  
  
}



double 
NSE_neutrino_Q_old(double kT, double rho, double Ye, int flavour, int iz, int in)
{

  int inuc,ii; 
  double T9, lg_rhoYe, rhoYe, lg_rho, lg_T, mu=0.0;
  double Q_weak_nue, R_weak_nue,    weak_rate_nue, avg_neutrino_energy_weak_nue , X_NSE_nue;
  double lambda_nue, lambda_nuebar;
  double Q_weak_nuebar, R_weak_nuebar, weak_rate_nuebar, avg_neutrino_energy_weak_nuebar, X_NSE_nuebar;

  const  double proton_mass = 1.67262158E-24, qe=1.6021892E-6;
  double eff_rate, Q_value_eff, capture_rate, decay_rate; 
  double Q_total=0.0;

  lg_rho=log10(rho);
  lg_rhoYe = log10(rho*Ye);
  rhoYe=rho*Ye;
  T9=kT*11.6045;
  lg_T = log10(T9*1E9);



  
  if(flavour==1){ 
    inuc = ffn_inuc_nue_ZN[iz][in];
    if(inuc==-1) return 0.0;  // nuclei not in FFN tables
    X_NSE_nue    =  NSE(T9,rho,Ye,iz,in);
    if(X_NSE_nue==0.0) return 0.0; // nuclei not in NSE
    avg_neutrino_energy_weak_nue = ffn_rate_nue(T9,lg_rhoYe,inuc,1);
    weak_rate_nue                = ffn_rate_nue(T9,lg_rhoYe,inuc,0);
    lambda_nue=weak_rate_nue*X_NSE_nue/A_nue[inuc];
    R_weak_nue=weak_rate_nue*X_NSE_nue*rho/proton_mass/A_nue[inuc];
    Q_weak_nue=R_weak_nue*avg_neutrino_energy_weak_nue*qe;
    capture_rate             = ffn_rate_nue(T9,lg_rhoYe,inuc,3);
    decay_rate               = ffn_rate_nue(T9,lg_rhoYe,inuc,2);
    Q_total = Q_total+ Q_weak_nue;

  }

  if(flavour==2){ /* antineutrinos */
    inuc = ffn_inuc_nuebar_ZN[iz][in];
    if(inuc==-1) return 0.0; // nuclei not in FFN tables
    X_NSE_nuebar    =  NSE(T9,rho,Ye,iz,in);
    if(X_NSE_nuebar==0.0) return 0.0; // nuclei not in NSE
    avg_neutrino_energy_weak_nuebar = ffn_rate_nuebar(T9,lg_rhoYe,inuc,1);
    weak_rate_nuebar                = ffn_rate_nuebar(T9,lg_rhoYe,inuc,0);
    lambda_nuebar=weak_rate_nuebar*X_NSE_nuebar/A_nuebar[inuc];
    R_weak_nuebar=weak_rate_nuebar*X_NSE_nuebar*rho/proton_mass/A_nuebar[inuc];
    Q_weak_nuebar=R_weak_nuebar*avg_neutrino_energy_weak_nuebar*qe;
    capture_rate             = ffn_rate_nuebar(T9,lg_rhoYe,inuc,3);
    decay_rate               = ffn_rate_nuebar(T9,lg_rhoYe,inuc,2);
    Q_total = Q_total+ Q_weak_nuebar;

  }

  return Q_total;


}


double 
NSE_neutrino_R_old(double kT, double rho, double Ye, int flavour, int iz, int in)
{

  int inuc,ii; 
  double T9, lg_rhoYe, rhoYe, lg_rho, lg_T, mu=0.0;
  double Q_weak_nue, R_weak_nue,    weak_rate_nue, avg_neutrino_energy_weak_nue , X_NSE_nue;
  double lambda_nue, lambda_nuebar;
  double Q_weak_nuebar, R_weak_nuebar, weak_rate_nuebar, avg_neutrino_energy_weak_nuebar, X_NSE_nuebar;

  const  double proton_mass = 1.67262158E-24, qe=1.6021892E-6;
  double eff_rate, Q_value_eff, capture_rate, decay_rate; 
  double R_total=0.0;

  lg_rho=log10(rho);
  lg_rhoYe = log10(rho*Ye);
  rhoYe=rho*Ye;
  T9=kT*11.6045;
  lg_T = log10(T9*1E9);

  if(flavour==1){ 
    inuc = ffn_inuc_nue_ZN[iz][in];
    if(inuc==-1) return 0.0;
    X_NSE_nue    =  NSE(T9,rho,Ye,iz,in);
    if(X_NSE_nue==0.0) return 0.0;
    
    avg_neutrino_energy_weak_nue = ffn_rate_nue(T9,lg_rhoYe,inuc,1);
    weak_rate_nue                = ffn_rate_nue(T9,lg_rhoYe,inuc,0);
    lambda_nue=weak_rate_nue*X_NSE_nue/A_nue[inuc];
    R_weak_nue=weak_rate_nue*X_NSE_nue*rho/proton_mass/A_nue[inuc];
    Q_weak_nue=R_weak_nue*avg_neutrino_energy_weak_nue*qe;
    capture_rate             = ffn_rate_nue(T9,lg_rhoYe,inuc,3);
    decay_rate               = ffn_rate_nue(T9,lg_rhoYe,inuc,2);
    R_total = R_total+ R_weak_nue;

  }

  if(flavour==2){ /* antineutrinos */
    inuc = ffn_inuc_nuebar_ZN[iz][in];
    if(inuc==-1) return 0.0;
    X_NSE_nuebar    =  NSE(T9,rho,Ye,iz,in);
    if(X_NSE_nuebar==0.0) return 0.0;
    avg_neutrino_energy_weak_nuebar = ffn_rate_nuebar(T9,lg_rhoYe,inuc,1);
    weak_rate_nuebar                = ffn_rate_nuebar(T9,lg_rhoYe,inuc,0);
    lambda_nuebar=weak_rate_nuebar*X_NSE_nuebar/A_nuebar[inuc];
    R_weak_nuebar=weak_rate_nuebar*X_NSE_nuebar*rho/proton_mass/A_nuebar[inuc];
    Q_weak_nuebar=R_weak_nuebar*avg_neutrino_energy_weak_nuebar*qe;
    capture_rate             = ffn_rate_nuebar(T9,lg_rhoYe,inuc,3);
    decay_rate               = ffn_rate_nuebar(T9,lg_rhoYe,inuc,2);
    R_total = R_total+ R_weak_nuebar;

  }

  return R_total;


}





double
NSE_neutrino_Q(double kT, double rho, double Ye, double mu, int flavour, int iz, int in)
{
  
  int inuc,ii; 
  double T9, lg_rhoYe, rhoYe, lg_rho, lg_T;
  double Q_weak_nue, R_weak_nue,    weak_rate_nue, avg_neutrino_energy_weak_nue , X_NSE_nue;
  double lambda_nue, lambda_nuebar;
  double Q_weak_nuebar, R_weak_nuebar, weak_rate_nuebar, avg_neutrino_energy_weak_nuebar, X_NSE_nuebar;
  
  const  double proton_mass = 1.67262158E-24, qe=1.6021892E-6;
  
  double eff_rate, Q_value_eff, capture_rate, decay_rate; 
  double Enu, dEnu, Enu_max=100.0, Enu_min=0.0;
  int num=1000;
  double spectrum=0.0;
  
  dEnu=(Enu_max-Enu_min)/((double) num);
  
  lg_rho=log10(rho);
  lg_rhoYe = log10(rho*Ye);
  rhoYe=rho*Ye;
  T9=kT*11.6045;
  lg_T = log10(T9*1E9);
  //  mu = chemical_potential(kT, rhoYe, 1E-12);  
  //mu = chemical_potential_tbl(lg_T, lg_rhoYe);  
  
  if (flavour==1) /* electron neutrinos */
  {
    inuc = ffn_inuc_nue_ZN[iz][in];
    if(inuc==-1) return 0.0;  // nuclei not in FFN tables
    X_NSE_nue    =  NSE(T9,rho,Ye,iz,in);
    if(X_NSE_nue==0.0) return 0.0; // nuclei not in NSE
    avg_neutrino_energy_weak_nue = ffn_rate_nue(T9,lg_rhoYe,inuc,1);
    weak_rate_nue                = ffn_rate_nue(T9,lg_rhoYe,inuc,0);
    lambda_nue=weak_rate_nue*X_NSE_nue/A_nue[inuc];
    R_weak_nue=weak_rate_nue*X_NSE_nue*rho/proton_mass/A_nue[inuc];
    Q_weak_nue=R_weak_nue*avg_neutrino_energy_weak_nue*qe;
    capture_rate             = ffn_rate_nue(T9,lg_rhoYe,inuc,3);
    decay_rate               = ffn_rate_nue(T9,lg_rhoYe,inuc,2);
    
    Q_value_eff = ffn_rate_nue(T9,lg_rhoYe,inuc,5);
    eff_rate    = ffn_rate_nue(T9,lg_rhoYe,inuc,4);
    
    if( (capture_rate>=decay_rate)  && (avg_neutrino_energy_weak_nue>3.0*kT) )  
      for(ii=0;ii<=num;ii++)
      {
        Enu=Enu_min+ii*dEnu;
        spectrum = spectrum +  Enu*dEnu*capture(Enu,Q_value_eff, kT,  mu)/eff_rate*X_NSE_nue*rho/proton_mass/A_nue[inuc];
      }
    else
      for(ii=0;ii<=num;ii++)
      {
        Enu=Enu_min+ii*dEnu;
        spectrum = spectrum +  Enu*dEnu*decay(Enu,Q_value_eff, kT,  -mu)/eff_rate*X_NSE_nue*rho/proton_mass/A_nue[inuc];
      }
    return spectrum*qe; /* electron neutrinos */      
  }
  
  if (flavour==2) /* electron antineutrinos */
  { 
    inuc = ffn_inuc_nuebar_ZN[iz][in];
    if(inuc==-1) return 0.0;  // nuclei not in FFN tables
      X_NSE_nuebar    =  NSE(T9,rho,Ye,iz,in);
    if(X_NSE_nuebar==0.0) return 0.0; // nuclei not in NSE
    avg_neutrino_energy_weak_nuebar = ffn_rate_nuebar(T9,lg_rhoYe,inuc,1);
    weak_rate_nuebar                = ffn_rate_nuebar(T9,lg_rhoYe,inuc,0);
    lambda_nuebar=weak_rate_nuebar*X_NSE_nuebar/A_nuebar[inuc];
    R_weak_nuebar=weak_rate_nuebar*X_NSE_nuebar*rho/proton_mass/A_nuebar[inuc];
    Q_weak_nuebar=R_weak_nuebar*avg_neutrino_energy_weak_nuebar*qe;
    capture_rate             = ffn_rate_nuebar(T9,lg_rhoYe,inuc,3);
    decay_rate               = ffn_rate_nuebar(T9,lg_rhoYe,inuc,2);
    
    eff_rate = ffn_rate_nuebar(T9,lg_rhoYe,inuc,4);
    Q_value_eff = ffn_rate_nuebar(T9,lg_rhoYe,inuc,5);     
    
    if( (capture_rate>=decay_rate)  && (avg_neutrino_energy_weak_nuebar>3.0*kT) )  
      for(ii=0;ii<=num;ii++)
      {
        Enu=Enu_min+ii*dEnu;
        spectrum = spectrum +  Enu*dEnu*capture(Enu,Q_value_eff, kT,  -mu)/eff_rate*X_NSE_nuebar*rho/proton_mass/A_nuebar[inuc];
      }
    else
      for(ii=0;ii<=num;ii++)
      {
        Enu=Enu_min+ii*dEnu;
        spectrum = spectrum +  Enu*dEnu*decay(Enu,Q_value_eff, kT,  mu)/eff_rate*X_NSE_nuebar*rho/proton_mass/A_nuebar[inuc];
      }

    return spectrum*qe; /* electron antineutrinos */    
  }

  return 0.0;

}


double 
NSE_neutrino_R(double kT, double rho, double Ye, double mu, int flavour, int iz, int in)
{
  
  int inuc,ii; 
  double T9, lg_rhoYe, rhoYe, lg_rho, lg_T;
  double Q_weak_nue, R_weak_nue,    weak_rate_nue, avg_neutrino_energy_weak_nue , X_NSE_nue;
  double lambda_nue, lambda_nuebar;
  double Q_weak_nuebar, R_weak_nuebar, weak_rate_nuebar, avg_neutrino_energy_weak_nuebar, X_NSE_nuebar;
  
  const  double proton_mass = 1.67262158E-24, qe=1.6021892E-6;
  
  double eff_rate, Q_value_eff, capture_rate, decay_rate; 
  double Enu, dEnu, Enu_max=100.0, Enu_min=0.0;
  int num=1000;
  double spectrum=0.0;
  
  dEnu=(Enu_max-Enu_min)/((double) num);
  
  lg_rho=log10(rho);
  lg_rhoYe = log10(rho*Ye);
  rhoYe=rho*Ye;
  T9=kT*11.6045;
  lg_T = log10(T9*1E9);
  //  mu = chemical_potential(kT, rhoYe, 1E-12);  
  //mu = chemical_potential_tbl(lg_T, lg_rhoYe);  
  
  if (flavour==1) /* electron neutrinos */
  {
    inuc = ffn_inuc_nue_ZN[iz][in];
    if(inuc==-1) return 0.0;  // nuclei not in FFN tables
      X_NSE_nue    =  NSE(T9,rho,Ye,iz,in);
    if(X_NSE_nue==0.0) return 0.0; // nuclei not in NSE
    avg_neutrino_energy_weak_nue = ffn_rate_nue(T9,lg_rhoYe,inuc,1);
    weak_rate_nue                = ffn_rate_nue(T9,lg_rhoYe,inuc,0);
    lambda_nue=weak_rate_nue*X_NSE_nue/A_nue[inuc];
    R_weak_nue=weak_rate_nue*X_NSE_nue*rho/proton_mass/A_nue[inuc];
    Q_weak_nue=R_weak_nue*avg_neutrino_energy_weak_nue*qe;
    capture_rate             = ffn_rate_nue(T9,lg_rhoYe,inuc,3);
    decay_rate               = ffn_rate_nue(T9,lg_rhoYe,inuc,2);
    
    Q_value_eff = ffn_rate_nue(T9,lg_rhoYe,inuc,5);
    eff_rate    = ffn_rate_nue(T9,lg_rhoYe,inuc,4);
    
    if( (capture_rate>=decay_rate)  && (avg_neutrino_energy_weak_nue>3.0*kT) )  
      for(ii=0;ii<=num;ii++)
      {
        Enu=Enu_min+ii*dEnu;
        spectrum = spectrum +  dEnu*capture(Enu,Q_value_eff, kT,  mu)/eff_rate*X_NSE_nue*rho/proton_mass/A_nue[inuc];
      }
    else
      for(ii=0;ii<=num;ii++)
      {
        Enu=Enu_min+ii*dEnu;
        spectrum = spectrum +  dEnu*decay(Enu,Q_value_eff, kT,  -mu)/eff_rate*X_NSE_nue*rho/proton_mass/A_nue[inuc];
      }  
    return spectrum; /* electron neutrinos */
  }
    
  if (flavour==2) /* electron antineutrinos */
  {
    inuc = ffn_inuc_nuebar_ZN[iz][in];
    if(inuc==-1) return 0.0;  // nuclei not in FFN tables
      X_NSE_nuebar    =  NSE(T9,rho,Ye,iz,in);
    if(X_NSE_nuebar==0.0) return 0.0; // nuclei not in NSE
    avg_neutrino_energy_weak_nuebar = ffn_rate_nuebar(T9,lg_rhoYe,inuc,1);
    weak_rate_nuebar                = ffn_rate_nuebar(T9,lg_rhoYe,inuc,0);
    lambda_nuebar=weak_rate_nuebar*X_NSE_nuebar/A_nuebar[inuc];
    R_weak_nuebar=weak_rate_nuebar*X_NSE_nuebar*rho/proton_mass/A_nuebar[inuc];
    Q_weak_nuebar=R_weak_nuebar*avg_neutrino_energy_weak_nuebar*qe;
    capture_rate             = ffn_rate_nuebar(T9,lg_rhoYe,inuc,3);
    decay_rate               = ffn_rate_nuebar(T9,lg_rhoYe,inuc,2);
    
    eff_rate = ffn_rate_nuebar(T9,lg_rhoYe,inuc,4);
    Q_value_eff = ffn_rate_nuebar(T9,lg_rhoYe,inuc,5);     
    
    if( (capture_rate>=decay_rate)  && (avg_neutrino_energy_weak_nuebar>3.0*kT) )  
      for(ii=0;ii<=num;ii++)
       {
        Enu=Enu_min+ii*dEnu;
        spectrum = spectrum +  dEnu*capture(Enu,Q_value_eff, kT,  -mu)/eff_rate*X_NSE_nuebar*rho/proton_mass/A_nuebar[inuc];
       }
    else
      for(ii=0;ii<=num;ii++)
       {
        Enu=Enu_min+ii*dEnu;
        spectrum = spectrum +  dEnu*decay(Enu,Q_value_eff, kT,  mu)/eff_rate*X_NSE_nuebar*rho/proton_mass/A_nuebar[inuc];
       }  
    /* end if-else */
    return spectrum; /* electron antineutrinos */
  }


  return 0.0;

}

double 
capture_integral(double kT, double mu, double Q)
{
  double eta,theta;
  double fd,fdd1,fdd2;
  // Fermi-Dirac integrals degree required
  double k_12=0.5,k_32=1.5,k_52=2.5,k_72=3.5;
  double rate_capture=0.0;
  const double me=0.511;
  double Enu,dEnu,Enu_min,Enu_max=100.0;;
  int ii,num=100;
  
  
  eta=(mu-me)/kT;
  theta=kT/me;
  
  //dfermi_(&k_12,&eta,&theta,&fd,&fdd1,&fdd2);
  fd = Ffermi(k_12,eta,theta);
  
  rate_capture = rate_capture+(me*me*me+2.0*me*me*Q+me*Q*Q)*fd;
  
  //dfermi_(&k_32,&eta,&theta,&fd,&fdd1,&fdd2);
  fd = Ffermi(k_32,eta,theta);
  
  rate_capture = rate_capture+(3.0*kT*me*me+4.0*kT*me*Q+kT*Q*Q)*fd;
  
  //dfermi_(&k_52,&eta,&theta,&fd,&fdd1,&fdd2);
  fd = Ffermi(k_52,eta,theta);
  
  rate_capture = rate_capture+(3.0*kT*kT*me+2.0*kT*kT*Q)*fd;
  
  //dfermi_(&k_72,&eta,&theta,&fd,&fdd1,&fdd2);
  fd = Ffermi(k_72,eta,theta);
  
  rate_capture = rate_capture+kT*kT*kT*fd;
  
  rate_capture = rate_capture*sqrt(2.0*me)*pow(kT,1.5);
  
  if(Q+me>=0.0)
   {
    return rate_capture;
   }
  else
  {
    Enu_min=Q+me;
    dEnu=-Enu_min/((double) num);
    for(ii=1;ii<num;ii++)
    {
      Enu=Enu_min+ii*dEnu;
      rate_capture = rate_capture -  dEnu*capture(Enu, Q, kT, mu);
    }
  }
  
  return rate_capture;
}