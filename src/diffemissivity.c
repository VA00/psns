#include <cuba.h>
#include "thermal.h"
#include "constants.h"
#include <stdio.h>
#include <math.h>

#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.0000001
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0




double Emax; 
double m;
double CV;
double CA;
double G;
double kT;
double mu;
double k0;




enum Integrands { PAIR_1, PAIR_2 , PHOTO_1, PHOTO_2};

static int Integrand_pair_1(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
{

  ff[0] = integr_pair1( (Emax - m) * xx[0] + m, (Emax - m) * xx[1] + m, k0, 
			 2.0 * xx[2] - 1.0, m, CV, CA, kT, mu)
		* (Emax - m) * (Emax - m) * 2.0;  
}


//--------------------------------------------------------------------------------
static int Integrand_pair_2(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
{

  ff[0] = integr_pair2( (Emax - m) * xx[0] + m, (Emax - m) * xx[1] + m, k0, 
			 2.0 * xx[2] - 1.0, m, CV, CA, kT, mu)
		* (Emax - m) * (Emax - m) * 2.0;  
}


static int Integrand_photo_1(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
{

ff[0] = integr_photo(2.0 * xx[0] - 1.0 , k0, Emax*xx[1], 2.0 * xx[2] - 1.0, 2.0*M_PI*xx[3],
    2.0*M_PI*xx[4], Emax * xx[5], Emax*xx[6], kT, mu, G, ee, CV, CA, m)* Emax*(Emax*Emax) * 2.0*2.0*2.0*M_PI*2.0*M_PI;  ;

}

static int Integrand_photo_2(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
{

ff[0] = integr_photo(2.0 * xx[0] - 1.0 ,  Emax*xx[1], k0, 2.0 * xx[2] - 1.0, 2.0*M_PI*xx[3],
    2.0*M_PI*xx[4], Emax * xx[5], Emax*xx[6], kT, mu, G, ee, CV, CA, m)* Emax * (Emax*Emax) * 2.0*2.0*2.0*M_PI*2.0*M_PI;  ;

}







//-------------------------------------------------------------------------------------------------------------------

static int (*GetIntegrand(int num, int *dim, double *epsrel, 
	        double p1, double p2, double p3, double p4, double p5, double p6,double p7))
		(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
{
    
switch (num) {
	
  case PAIR_1:
    Emax = p1;
    m    = p2;
    CV   = p3;
    CA   = p4;
    kT   = p5;
    mu   = p6;
    k0   = p7;
    *dim = 3;
    *epsrel = 0.01;
  return &Integrand_pair_1;
    	
  case PAIR_2 :
    Emax = p1;
    m    = p2;
    CV   = p3;
    CA   = p4;
    kT   = p5;
    mu   = p6;
    k0   = p7;
    *dim = 3;
    *epsrel = 0.01;
  return &Integrand_pair_2; 
	
  case PHOTO_1:
    Emax = p1;
    m    = p2;
    CV   = p3;
    CA   = p4;
    kT   = p5;
    mu   = p6;
    k0   = p7;
    *dim = 7;
    *epsrel = 0.01;
  return &Integrand_photo_1;  				    
	
  case PHOTO_2:
    Emax = p1;
    m    = p2;
    CV   = p3;
    CA   = p4;
    kT   = p5;
    mu   = p6;
    k0   = p7;
    *dim = 7;
    *epsrel = 0.01;
  return &Integrand_photo_2;
}

return 0;
}


double pair(double Enu, double kT, double mu, int flavour, int method)
{

int KEY=0;
int NSTART=1000;
int NINCREASE=500;
int CUBA_VERBOSE=0;
int LAST=4;
int MINEVAL=2000;
int MAXEVAL=500000;

int NNEW=1000;
double FLATNESS=25.0;

double sin2_theta_W = SQUARED_SIN_WEINBERG_ANGLE;
double E_cutoff = 10.0;

double CV = 0.5 + 2.0 * sin2_theta_W;
double CA = 0.5;
double m = 0.511;
                
integrand_t Integrand_pair;
int nregions, neval, fail;
int NDIM, NCOMP=1;
double EPSREL=0.1, EPSABS=1e-100;
int INT_METHOD=4;
int PRINT_INT=0;
double integral[NCOMP], error[NCOMP], prob[NCOMP];
int C=PAIR_2;
double emmi=-1.0;
	
        
//1 - nu_e, 2-anty-nu_e, 3 - nu_mutau, 4- anty-nu_mutau
switch(flavour){
  case 1: CV =  0.5 + 2.0 * sin2_theta_W; C = PAIR_1;
  break;

  case 2: CV =  0.5 + 2.0 * sin2_theta_W; C = PAIR_2;
  break;
	
  case 3: CV = 0.5 - 2.0 * sin2_theta_W; CA=0.5; C = PAIR_1;
  break;

  case 4: CV = 0.5 - 2.0 * sin2_theta_W; CA=0.5; C = PAIR_2;
  break;
};

switch(method){
  //1 - FAST (DIVONNE, small accuracy)
  case 1: INT_METHOD=4;EPSREL=0.1;
  break;
  //  - NORMAL (VEGAS)
  case 2: INT_METHOD=2;EPSREL=0.01;MINEVAL=2000;MAXEVAL=50000;
  break;
  // - ACCURATE (CUHRE)
  case 3: INT_METHOD=1;EPSREL=0.03;MINEVAL=200000;MAXEVAL=1000000;
  break;
};

//Simple but working cutoff estimate for incoming electron and positron
E_cutoff = 10.0*kT + 2.0*mu+2.0*Enu;

Integrand_pair = GetIntegrand(C, &NDIM, &EPSREL, E_cutoff, m, CV, CA, kT, mu, Enu);

switch(INT_METHOD){

  case 1:  
    Cuhre(NDIM, NCOMP, Integrand_pair, NULL, 1, EPSREL, EPSABS, CUBA_VERBOSE | LAST, 
    MINEVAL, MAXEVAL, KEY, NULL, NULL, &nregions, &neval, &fail, integral, error, prob);
	
    if(PRINT_INT){	
      printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n", nregions, neval, fail);
      printf("For k0=%lf:\t%e\t+-\t%e\n", Enu , integral[0], error[0]);
    };
  break;

  case 2: 
    Vegas(NDIM,NCOMP,Integrand_pair, NULL, 1,
		EPSREL, EPSABS, 
		CUBA_VERBOSE, 0, MINEVAL, MAXEVAL,	
		NSTART, NINCREASE, 1000, 0, NULL, NULL,
		&neval, &fail, integral, error, prob
    );
/*    
    Vegas(NDIM, NCOMP, Integrand, USERDATA,
		EPSREL, EPSABS, 
		verbose, SEED, MINEVAL, MAXEVAL, 
		NSTART, NINCREASE, NBATCH, GRIDNO, STATEFILE,
		&neval, &fail, integral, error, prob);
*/		
    if(PRINT_INT){	
      printf("VEGAS RESULT:\tnregions %d\tneval %d\tfail %d\n", nregions, neval, fail);
      printf("For k0=%.2f:\t%e\t+-\t%e\n", Enu , integral[0], error[0]);
    };
  break;

  case 3:
    Suave(NDIM, NCOMP, Integrand_pair, NULL, 1,
    EPSREL, EPSABS, CUBA_VERBOSE | LAST, 0, MINEVAL, MAXEVAL,
    NNEW, 2, FLATNESS, NULL, NULL,
    &nregions, &neval, &fail, integral, error, prob);
    
    /*
Cuba version 3.0 changes 
  Suave(NDIM, NCOMP, Integrand, USERDATA,  
		EPSREL, EPSABS, verbose | LAST, SEED, MINEVAL, MAXEVAL, 
		NNEW, FLATNESS,
		&nregions, &neval, &fail, integral, error, prob);

*/

    if(PRINT_INT){	
      printf("SUAVE RESULT:\tnregions %d\tneval %d\tfail %d\n", nregions, neval, fail);
      printf("For k0=%.2f:\t%e\t+-\t%e\n", Enu , integral[0], error[0]);
    };
  break;

  case 4:
    Divonne(NDIM, NCOMP, 
	    Integrand_pair, NULL, 1,
    EPSREL, EPSABS, 
    CUBA_VERBOSE, 0, 
    MINEVAL, MAXEVAL,
    KEY1, KEY2, KEY3, 
    MAXPASS, BORDER, 
    MAXCHISQ, MINDEVIATION,
    NGIVEN, LDXGIVEN, NULL, 
    NEXTRA, NULL, 
    NULL, NULL,
    &nregions, &neval, &fail, integral, error, prob);
   /* Cuba version 3.0 changes 
		   Divonne(NDIM, NCOMP, Integrand, USERDATA,
    EPSREL, EPSABS, verbose, SEED,    MINEVAL, MAXEVAL, 
    KEY1, KEY2, KEY3, MAXPASS,  BORDER, MAXCHISQ, MINDEVIATION,
    NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
    &nregions, &neval, &fail, integral, error, prob);
    
*/ 
    if(PRINT_INT){	
      printf("DIVONNE RESULT:\tnregions %d\tneval %d\tfail %d\n", nregions, neval, fail);
      printf("For k0=%.2f:\t%e\t+-\t%e\n", Enu , integral[0], error[0]);
    };
  break;
}

emmi = integral[0];	



//GF^2 me^9 c^7 / hbar^9 / qe^9 /10^54
return 1.74616E27*emmi;
}

double photo(double Enu, double kT, double mu, int flavour, int method)
{


  const int KEY=0,NSTART=1000,NINCREASE=500;
  const int CUBA_VERBOSE=0,LAST=4;
  int MINEVAL=20000,MAXEVAL=500000;
  const int PRINT_INT=0;   
  int INT_METHOD=2;

  const int NNEW=1000;
  const double FLATNESS=25.0;

  integrand_t Integrand_photo;
  int nregions, neval, fail;
  int NDIM=7, NCOMP=1;
  double EPSREL=0.1, EPSABS=1e-100;
  double integral[NCOMP], error[NCOMP], prob[NCOMP];

  double sin2_theta_W = 0.22280;
  double CV = 0.5 + 2.0 * sin2_theta_W, CA = 0.5;
  double m = 0.511;
  double E_cutoff;
  int C=PHOTO_2;
  double emmi=0.0;
	
        
//1 - nu_e, 2-anty-nu_e, 3 - nu_mutau, 4- anty-nu_mutau
switch(flavour){

 case 1: {CV =  0.5 + 2.0 * sin2_theta_W; C = PHOTO_1;};
 break;
 case 2: {CV =  0.5 + 2.0 * sin2_theta_W; C = PHOTO_2;};
 break;
 case 3: {CV = -0.5 + 2.0 * sin2_theta_W; C = PHOTO_1;};
 break;
 case 4: {CV = -0.5 + 2.0 * sin2_theta_W; C = PHOTO_2;};
 break;

};


switch(method){

  //1 - FAST (DIVONNE, small accuracy)
	case 1: INT_METHOD=4;EPSREL=0.1;
	break;
  //  - NORMAL (VEGAS)
	case 2: INT_METHOD=2;EPSREL=0.005;MINEVAL=2000;MAXEVAL=50000;
	break;
  // - ACCURATE (CUHRE)
	case 3: INT_METHOD=1;EPSREL=0.05;MINEVAL=20000;MAXEVAL=1000000;
	break;
	
};



//Simple but working cutoff estimate for incoming electron and positron
E_cutoff = 4.0 * kT + 1.0*mu + Enu;



	

Integrand_photo = 
GetIntegrand( C , &NDIM, &EPSREL, E_cutoff, m, CV, CA, kT, mu, Enu);



switch(INT_METHOD){

  case 1:  Cuhre( NDIM, NCOMP, Integrand_photo, NULL, 1,
                  EPSREL, EPSABS, CUBA_VERBOSE | LAST, 
         	  MINEVAL, MAXEVAL, KEY, NULL, NULL,
         	  &nregions, &neval, &fail, integral, error, prob
           );
/*	
	     Cuhre(NDIM, NCOMP, Integrand, USERDATA,
    EPSREL, EPSABS, verbose | LAST,
    MINEVAL, MAXEVAL, KEY,
    &nregions, &neval, &fail, integral, error, prob);
*/

  if(PRINT_INT){	
        printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n", 
        nregions, neval, fail);
        printf("For k0=%.2f:\t%e\t+-\t%e\t%lf\n", 
        Enu , integral[0], error[0], error[0]/integral[0]);
  };
  
  break;

  case 2:  EPSREL=0.005;
  	   Vegas(NDIM, NCOMP, Integrand_photo, NULL,1,
  	         EPSREL, EPSABS, 
  	         CUBA_VERBOSE, 0, MINEVAL, MAXEVAL,
    	         NSTART, NINCREASE, 1000, 0, NULL, NULL,
    	         &neval, &fail, integral, error, prob
    	   );
/*	   Cuba 3.0 changes: 

	   Vegas(NDIM, NCOMP, Integrand, USERDATA,
		EPSREL, EPSABS, 
		verbose, SEED, MINEVAL, MAXEVAL, 
		NSTART, NINCREASE, NBATCH, GRIDNO, STATEFILE,
		&neval, &fail, integral, error, prob);
*/

  if(PRINT_INT){ printf("VEGAS RESULT:\tneval %d\tfail %d\n", neval, fail);
                 printf("For k0=%.2f:\t%e\t+-\t%e\t%lf\n", 
                         Enu , integral[0], error[0], error[0]/integral[0]);
  };

  break;


  case 3:  Suave(NDIM, NCOMP, Integrand_photo, NULL, 1,
                 EPSREL, EPSABS, 0, 0, MINEVAL, MAXEVAL,
                 NNEW, 2, FLATNESS, NULL , NULL,
                 &nregions, &neval, &fail, integral, error, prob);
/*
Cuba version 3.0 changes 
	  Suave(NDIM, NCOMP, Integrand, USERDATA,  
		EPSREL, EPSABS, verbose | LAST, SEED, MINEVAL, MAXEVAL, 
		NNEW, FLATNESS,
		&nregions, &neval, &fail, integral, error, prob);

*/

if(PRINT_INT){	
printf("SUAVE RESULT:\tnregions %d\tneval %d\tfail %d\n", nregions, neval, fail);
printf("For k0=%.2f:\t%e\t+-\t%e\n", Enu , integral[0], error[0]);
};

break;



  case 4:EPSREL=0.01;
         Divonne(NDIM, NCOMP, Integrand_photo, NULL, 1,
    		 EPSREL, EPSABS, CUBA_VERBOSE, 0, MINEVAL, MAXEVAL,
                 KEY1, KEY2, KEY3, MAXPASS, BORDER, MAXCHISQ, MINDEVIATION,
                 NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
		 NULL, NULL,
                 &nregions, &neval, &fail, integral, error, prob);
/* Cuba version 3.0 changes 
		   Divonne(NDIM, NCOMP, Integrand, USERDATA,
    EPSREL, EPSABS, verbose, SEED,
    MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
    BORDER, MAXCHISQ, MINDEVIATION,
    NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
    &nregions, &neval, &fail, integral, error, prob);
    
*/

  if(PRINT_INT){	
                 printf("DIVONNE RESULT:\tnregions %d\tneval %d\tfail %d\n", 
                 nregions, neval, fail);
                 printf("For k0=%.2f:\t%e\t+-\t%e\t%lf\n", 
                 Enu , integral[0], error[0], error[0]/integral[0]);
  };

break;



}




emmi = 1.0E22*integral[0];	


return emmi;


}



/*
method 1: analytical formula using omegaL=omega0 (like Adams, Rudderman, Woo)
method 2: Braaten & Segel approximation for dispersion relations
TODO method 3: full polarization functions (B&S gives underestimate error up to 70%
TODO in the pre-supernova core compared to Itoh!!!)
*/
double plasmaL(double Enu, double kT, double mu, int flavour)
{

 
  double f,x,CV=1.0;
  const int method=2;

  int i;
  double integral,integrand,dE2,E1,E2,ct,k,omegaL;
  double V;
  double Zl,Jl,betaL;

  double omega0, omega1, kmax;

  if(Enu<=0.0) return 0.0;
  
  omega0 = plasma_frequency(kT,mu);
  V = electron_velocity(kT,mu);
  omega1 = V*omega0;

  kmax = max_momentum(kT,mu);
              
  switch(flavour){
  case 1: {CV =  0.5 + 2.0 * SQUARED_SIN_WEINBERG_ANGLE;};
  break;
  
  case 2: {CV =  0.5 + 2.0 * SQUARED_SIN_WEINBERG_ANGLE;};
  break;
  
  case 3: {CV = -0.5 + 2.0 * SQUARED_SIN_WEINBERG_ANGLE;};
  break;
  
  case 4: {CV = -0.5 + 2.0 * SQUARED_SIN_WEINBERG_ANGLE;};
  break;
};
  

switch(method){
  case 1: {

  if(Enu>omega0) return 0.0;
  if(!(omega0>0.0)) return 0.0;

  x=Enu/omega0;  


//Singular point, limit computed instead
  if(x==0.5) f=1.0; else
  f=
    4.0*x*(x-1.0)*(8.0*x*x*x*x - 16.0*x*x*x + 2.0*x*x + 6.0*x - 3.0)
    + 
    3.0*(1.0-2.0*x)*(1.0-2.0*x)*log( fabs(1.0-2.0*x)*fabs(1.0-2.0*x) )
  ;

  
//GF^2 /c^9 /hbar^10 qe^8 10^48


  return 2.68697E31*CV*CV*pow(omega0,7.0)*105.0/32.0/1260.0
         /pow(M_PI,4.0)*ALPHA_INV*f/(exp(omega0/kT)-1.0);
};
  break;

  case 2: {


/*
Extremely simple integration using 200 points
*/
  E1=Enu;
  dE2=(kmax-omega0)/200.0;
  integral=0.0;
  
  for(i=1;i<200;i++){
    E2=omega0-E1 + i*dE2;
    
    if( (kmax>E1+E2) && (E1+E2>omega0) )     {
    k=inv_dispersionL(E1+E2,omega0,V,(kmax-omega0)/200.0*0.001);
    
    omegaL=E1+E2;
   
    ct = 0.5*(k*k-E1*E1-E2*E2)/E1/E2;
    Zl = 2.0*omegaL*omegaL/(omegaL*omegaL - k*k)
        *(omegaL*omegaL - V*V*k*k)/(3.0*omega0*omega0 - omegaL*omegaL + V*V * k*k);
    
    betaL=1.5*omega0*omega0/V/V/V*
     (1.5*omegaL/k/k/k*log((omegaL + V*k)/(omegaL - V*k)) - 
omegaL*omegaL/k/k*V/(omegaL*omegaL - V*V * k*k) - 2.0*V/k/k);
    Jl=fabs(k*k/E1/E2/omegaL*(1.0 - betaL)/betaL);
    integrand = Zl*Jl/(exp(omegaL/kT) - 1.0)*
         (omegaL*omegaL - k*k)*(omegaL*omegaL - k*k) 
           *E1*E1*E1*E2*E2*E2*(1.0 - ct*ct)/k/k/(E1 + E2);
    }
    
    else integrand=0.0;
    
    if(4.0*E1*E2-(E1+E2)*(E1+E2)+k*k<0.0) integrand=0.0; 
    integral = integral + integrand;
  }

// GF^2/16/pi^4/alpha/hbar^10/c^9 * erg_to_MeV^8
      return 2.3652059719158757E30*CV*CV*integral*dE2;
         ;};
  break;
}

return 0.0;

}




/*
Simple integral formula for spectrum,
replacing omega0 with mt helps a bit but
errors in pre-sn core are up to 50%,
but decreasing with degeneracy (where plasma dominates)
down to ~10%
TODO Braaten & Segel (disperions not yet implemented)
TODO full polarization ?? depands on how good B&S will work
*/
double plasmaT(double Enu, double kT, double mu, int flavour)
{

  double omega0,mt,x,ct,dct;
  double P1,P2,P3,nB,integral;
  int i;


  double CV=1.0;
  const double sin2_theta_W = SQUARED_SIN_WEINBERG_ANGLE;

  if(Enu<=0.0) return 0.0;

  omega0 = plasma_frequency(kT,mu);
  mt = photon_mass(kT,mu);
  x=Enu/omega0;
               
  switch(flavour){

	case 1: {CV =  0.5 + 2.0 * sin2_theta_W;};
	break;
	case 2: {CV =  0.5 + 2.0 * sin2_theta_W;};
	break;
	case 3: {CV = -0.5 + 2.0 * sin2_theta_W;};
	break;
	case 4: {CV = -0.5 + 2.0 * sin2_theta_W;};
	break;

  };

/*
Extremely simple integration using 200 points
*/
  dct=0.01;
  integral=0.0;
  
  for(i=0;i<200;i++){
  
    ct=-1.0 + i*dct; 
    P1 = 1.0 + 2.0 * (ct-1.0)*(ct-1.0)*x*x*(2.0*x*x-1.0);
    P2 = x*(ct-1.0)*(ct-1.0);
    P3 = 1.0 - 4.0*(ct-1.0)*ct*x*x+4.0*(ct-1.0)*(ct-1.0)*x*x*x*x;
    nB = 1.0/(   exp(  omega0/kT*(x+1.0/2.0/x/(1.0-ct)) ) - 1.0   );
    integral = integral + P1/P2/P3*nB;
  }

//GF^2 /c^9 /hbar^10 Meverg^8 /64 /pi^4/ alpha

return integral*dct*pow(omega0,7.0)*5.9145687147347446E29;
}


/*
    * $Id: diffemissivity.c,v 1.7 2007/02/12 09:06:26 cvsroot Exp $ 
    * $Log: diffemissivity.c,v $
    * Revision 1.7  2007/02/12 09:06:26  cvsroot
    * omega0^7 in plasmaT is indeed CORRECT. CGS coeff adjusted, now plasmaT works perfect. Some comments and TODO's added.
    *
    * Revision 1.6  2007/02/09 10:27:52  cvsroot
    * Braaten-Segel approximation for longitudal plasmon coded
    *
    * Revision 1.5  2007/02/07 20:47:27  cvsroot
    * omega0^8 replaced with omega0^7: check needed
    *
    * $Date: 2007/02/12 09:06:26 $
    * $Author: cvsroot $ 
    * $Revision: 1.7 $ 
*/

