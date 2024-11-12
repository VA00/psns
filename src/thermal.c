/*
Written by M. Misiaszek and A. Odrzywolek
*/
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <stdio.h>
#include <fermidirac.h>
#include "thermal.h"
#define NUM_GAUSS_POINTS 40

// Structure to hold Gaussian quadrature data
typedef struct {
    double   nodes[NUM_GAUSS_POINTS + 1];
    double weights[NUM_GAUSS_POINTS + 1];
} GaussQuadrature;



// Initialize Gaussian quadrature points and weights
static GaussQuadrature gauss_quad = {
    .nodes = {        
	 0.0,  // 0-index unused
	 1.95113832567939976543512341074545479e-2 , 
	 5.85044371524206686289933218834177944e-2 , 
	 9.74083984415845990632784501049369020e-2 ,
	 1.36164022809143886559241078000717067e-1 ,
	 1.74712291832646812559339048011286195e-1 ,
	 2.12994502857666132572388538666321823e-1 ,
	 2.50952358392272120493158816035004797e-1 ,
	 2.88528054884511853109139301434713898e-1 ,
	 3.25664370747701914619112943627358695e-1 ,
	 3.62304753499487315619043286358963588e-1 ,
	 3.98393405881969227024379642517533757e-1 ,
	 4.33875370831756093062386700363181958e-1 ,
	 4.68696615170544477036078364935808657e-1 ,
	 5.02804111888784987593672750367568003e-1 ,
	 5.36145920897131932019857253125400904e-1 ,
	 5.68671268122709784725485786624827158e-1 ,
	 6.00330622829751743154746299164006848e-1 ,
	 6.31075773046871966247928387289336863e-1 ,
	 6.60859898986119801735967122844317234e-1 ,
	 6.89637644342027600771207612438935266e-1 ,
	 7.17365185362099880254068258293815278e-1 ,
	 7.44000297583597272316540527930913673e-1 ,
	 7.69502420135041373865616068749026083e-1 ,
	 7.93832717504605449948639311738454358e-1 ,
	 8.16954138681463470371124994012295707e-1 ,
	 8.38831473580255275616623043902867064e-1 ,
	 8.59431406663111096977192123491656492e-1 ,
	 8.78722567678213828703773343639124407e-1 ,
	 8.96675579438770683194324071967395986e-1 ,
	 9.13263102571757654164733656150947478e-1 ,
	 9.28459877172445795953045959075453133e-1 ,
	 9.42242761309872674752266004500001735e-1 ,
	 9.54590766343634905493481517021029508e-1 ,
	 9.65485089043799251452273155671454998e-1 ,
	 9.74909140585727793385645230069136276e-1 ,
	 9.82848572738629070418288027709116473e-1 ,
	 9.89291302499755531026503167136631385e-1 ,
	 9.94227540965688277892063503664911698e-1 ,
	 9.97649864398237688899494208183122985e-1 ,
	 9.99553822651630629880080499094567184e-1 },
    .weights = {
        0.0,  // 0-index unused
	 3.90178136563066548112804392527540483e-2 ,  
	 3.89583959627695311986255247722608223e-2 ,
	 3.88396510590519689317741826687871658e-2 ,
	 3.86617597740764633270771102671566912e-2 ,
	 3.84249930069594231852124363294901384e-2 ,
	 3.81297113144776383442067915657362019e-2 ,
	 3.77763643620013974897749764263210547e-2 ,
	 3.73654902387304900267053770578386691e-2 ,
	 3.68977146382760088391509965734052192e-2 ,
	 3.63737499058359780439649910465228136e-2 ,
	 3.57943939534160546028615888161544542e-2 ,
	 3.51605290447475934955265923886968812e-2 ,
	 3.44731204517539287943642267310298320e-2 ,
	 3.37332149846115228166751630642387284e-2 ,
	 3.29419393976454013828361809019595361e-2 ,
	 3.21004986734877731480564902872506960e-2 ,
	 3.12101741881147016424428667206035518e-2 ,
	 3.02723217595579806612200100909011747e-2 ,
	 2.92883695832678476927675860195791396e-2 ,
	 2.82598160572768623967531979650145302e-2 ,
	 2.71882275004863806744187066805442598e-2 ,
	 2.60752357675651179029687436002692871e-2 ,
	 2.49225357641154911051178470032198023e-2 ,
	 2.37318828659301012931925246135684162e-2 ,
	 2.25050902463324619262215896861687390e-2 ,
	 2.12440261157820063887107372506131285e-2 ,
	 1.99506108781419989288919287151135633e-2 ,
	 1.86268142082990314287354141521572090e-2 ,
	 1.72746520562693063585842071312909998e-2 ,
	 1.58961835837256880449029092291785257e-2 ,
	 1.44935080405090761169620745834605500e-2 ,
	 1.30687615924013392937868258970563403e-2 ,
	 1.16241141207978269164667699954326348e-2 ,
	 1.01617660411030645208318503524069436e-2 ,
	 8.68394526926085842640945220403428135e-3 ,
	 7.19290476811731275267557086795650747e-3 ,
	 5.69092245140319864926910711716201847e-3 ,
	 4.18031312469489523673930420168135132e-3 ,
	 2.66353358951268166929353583166845546e-3 ,
	 1.14495000318694153454417194131563611e-3 }
};

double gaussian_integrate(double a, double b, double (*integrand)(double x, double kT, double mu, int flavour), double kT, double mu, int flavour) {
    double result = 0.0;

    for (int i = 1; i <= NUM_GAUSS_POINTS; i++) {
        double x = a + (b-a) * gauss_quad.nodes[i];
        result += integrand(x, kT, mu, flavour) * gauss_quad.weights[i];
    }

    return result * (b-a);
}



//Eq. (20) in PHYSICAL REVIEW D 74, 043006 (2006)
double M(int sign, int n, int m, double alpha, double beta, double CV, double CA)
{

  switch(sign){
    case -1:
    return (7.0*CV*CV - 2.0*CA*CA)*Gm(n/2.0 - 1.0/2.0, alpha, beta)*
    Gp(m/2.0 - 1.0/2.0, alpha, beta) +
  9.0*CV*CV*Gm(n/2.0, alpha, beta)*Gp(m/2.0, alpha, beta) +
  (CV*CV+CA*CA)*(4.0*Gm(n/2.0 + 1.0/2.0, alpha, beta)*Gp(m/2.0 + 1.0/2.0, alpha, beta)
 -      Gm(n/2.0 - 1.0/2.0, alpha, beta)*Gp(m/2.0 + 1.0/2.0, alpha, beta) - 
        Gm(n/2.0 + 1.0/2.0, alpha, beta)*Gp(m/2.0 - 1.0/2.0, alpha, beta)
        );
    break;
  
    case 1:
    return (7.0*CV*CV - 2.0*CA*CA)*Gp(n/2.0 - 1.0/2.0, alpha, beta)*
    Gm(m/2.0 - 1.0/2.0, alpha, beta) +
  9.0*CV*CV*Gp(n/2.0, alpha, beta)*Gm(m/2.0, alpha, beta) +
  (CV*CV+CA*CA)*(4.0*Gp(n/2.0 + 1.0/2.0, alpha, beta)*Gm(m/2.0 + 1.0/2.0, alpha, beta)
 -      Gp(n/2.0 - 1.0/2.0, alpha, beta)*Gm(m/2.0 + 1.0/2.0, alpha, beta) - 
        Gp(n/2.0 + 1.0/2.0, alpha, beta)*Gm(m/2.0 - 1.0/2.0, alpha, beta)
        );
    break;
  } 

  return 0.0;
}


/*
Compute density*electron fraction from kT and chem. pot. [MeV]
TODO rewrite using new G_p G_m functions
*/
double rhoYe(double kT, double mu)
{
  double alpha,beta; 
  /*Conversion factor C [g/cm^3] and electron rest mass [MeV]
  
  C = 8 Pi me^3 c^3/h^3*mu /. 
  {
   me -> 9.109383632044565226841657507`10.*^-31, electron mass
   h -> 132521403/200000000000000000000000000000000000000000, Planck constant
   c -> 299792458, speed of light
   mu -> 1.6605390665999999711685831308728`9.220249097464801*^-27 atomic mass unit
  }
  NOTE: above density is in kg/m^3; divide by 1000 to get g/cc
  
  */
  double C=2.92179643E9, // 8 Pi me^3 c^3/h^3*mu in g/cc
         m=510.9989461E-3; // electron mass in MeV
  		  
  alpha=m/kT; beta=mu/kT;

  return C*( Gm(0.0, alpha, beta)-Gp(0.0, alpha, beta) );
}	



/*
Return chemical potential in MeV
TODO some guess for mu_right required
TODO mu_right=10.0 MeV for psns and s15 is enough, but e.g. for Itoh table test 1000 MeV is required
TODO but this cause psns slowdown!!!
*/
double chemical_potential(double kT, double rho_Ye, double EPSABS)
{
	
// We start with mu=1.0 MeV, chem_pot assumed to be less than 10.0 MeV
  double mu=1.0,mu_left=0.0,mu_right=1000.0;

//Classical bisection
  while(mu_right-mu_left>EPSABS)
  {
    mu=0.5*mu_left + 0.5*mu_right;
    if(rhoYe(kT,mu)>rho_Ye) mu_right=mu; else mu_left=mu;  
  }	  

  return mu;
}


/*
Properties of plasmons, all in MeV
*/
double plasma_frequency(double kT, double mu)
{


  double m=0.511,alphaSS=1.0/137.0;
  double alpha=m/kT,beta=mu/kT;
  


  return sqrt(4.0*alphaSS/3.0/M_PI*m*m*(2.0*Gp(-0.5,alpha,beta)+2.0*Gm(-0.5,alpha,beta)
               +Gp(-1.5,alpha,beta)+Gm(-1.5,alpha,beta))
	 );

}


double photon_mass(double kT, double mu)
{
	
	double m=0.511,alphaSS=1.0/137.0;
	double alpha=m/kT,beta=mu/kT;
	
	
	
  return sqrt(4.0*alphaSS/M_PI*m*m*(Gp(-0.5,alpha,beta)+Gm(-0.5,alpha,beta) )
              );

}





double electron_velocity(double kT, double mu)
{
	double m=0.511,alphaSS=1.0/137.0;
        double alpha = m/kT, beta = mu/kT, omega0, omega1_squared;

  
  
  omega0 = plasma_frequency(kT,mu);
  
  omega1_squared = omega0*omega0 - 4.0*alphaSS/M_PI*m*m*( Gp(-2.5,alpha,beta) + Gm(-2.5,alpha,beta)  );
  
  return sqrt(omega1_squared/omega0/omega0);
  
}

double integrand_max_momentum(double p, double kT, double mu, int flavour)
{

  const double m=0.511, alpha = 1/137.0;
  double   v=p/sqrt(p*p+m*m);

  return 4.0*alpha/3.1415926*p*p/sqrt(p*p + m*m)*
				(1.0/v*log((1.0+v)/(1.0-v))-1.0)
				*
				(
				1.0/(1.0 + exp((sqrt(p*p+m*m) - mu)/kT))+1.0/(1.0 + exp((sqrt(p*p + m*m) + mu)/kT))
                 );
}

double max_momentum(double kT, double mu)
{
  return sqrt(gaussian_integrate(0.0, 12.0*kT+mu, &integrand_max_momentum, kT, mu, 0));
}



//Braaten Segel approximation for dispersion relations
double dispersionL(double k, double kmax, double omega0, double V, double EPSABS)
{ 

  double kmax_BS,LHS,RHS;
  double omegaL=omega0,omegaL_left=omega0,omegaL_right=kmax;  
  kmax_BS = sqrt(3.0/V/V*(0.5/V*log((1.0+V)/(1.0-V))-1.0))*omega0;
  if(k>kmax_BS) return 0.0;
  if(k<EPSABS)  return omega0;
 
//Classical bisection


  if(k>omega0) omegaL_left=k;
  while(omegaL_right-omegaL_left>EPSABS)
  {
    omegaL=0.5*omegaL_left + 0.5*omegaL_right;
    LHS=V*V*k*k/3.0/omega0/omega0;
    RHS=(0.5*omegaL/V/k*log((omegaL+V*k)/(omegaL-V*k))-1.0);

    if(LHS>RHS) omegaL_right=omegaL; else omegaL_left=omegaL;  
  }	  

  return omegaL;
 

}			


double dispersionT(double k, double mt, double omega0, double V, double EPSABS)
{ 
	
return sqrt(k*k+ mt*mt); 

}


double inv_dispersionL(double omegaL, double omega0, double V, double EPSABS)
{ 


  double kmax_BS,LHS,RHS;
  double k=0.0,k_left=0.0,k_right=1.0;  
  kmax_BS = sqrt(3.0/V/V*(0.5/V*log((1.0+V)/(1.0-V))-1.0))*omega0;
  
  
k_right=omegaL/V; 
//Classical bisection
  while(k_right-k_left>EPSABS)
  {
    k=0.5*k_left + 0.5*k_right;
    LHS=V*V*k*k/3.0/omega0/omega0;
    RHS=(0.5*omegaL/V/k*log((omegaL+V*k)/(omegaL-V*k))-1.0);

    if(LHS<RHS) k_right=k; else k_left=k;  
  }	  

  return k;
 
 

}			


double inv_dispersionT(double omegaT, double omega0, double V, double EPSABS)
{ 

	return sqrt(omegaT*omegaT-omega0*omega0); 

}			



/*
Returns single flavour Q: 
1- electron neutrino emissivity.
2- electron ant-neutrino emissivity.
3- mu=tau neutrino emissivity.
4- mu=tau anti-neutrino emissivity. 
deltaQ is formula (41) in M. Misiaszek et. al. Phys. Rev. D 74, 043006 (2006) 
http://link.aps.org/abstract/PRD/v74/e043006
*/ 
double pair_Q(double kT, double mu, int flavour)
{ 

  const double sin2_theta_W=0.22280,m=0.511;
  double CV=1.0,CA=0.5;
  double alpha,beta;
  double M_10_p,M_10_m;
  double Q_half_avg, deltaQ, cgs_coeff;		  
  
  switch(flavour){
    case 1: CV =  0.5 + 2.0 * sin2_theta_W;
    break;

    case 2: CV =  0.5 + 2.0 * sin2_theta_W;
    break;
	
    case 3: CV = -0.5 + 2.0 * sin2_theta_W; CA=-0.5;
    break;

    case 4: CV = -0.5 + 2.0 * sin2_theta_W; CA=-0.5;
    break;
  };

  alpha=m/kT; beta=mu/kT;
  

  M_10_m =M(-1,1,0,alpha,beta,CV,CA);
  M_10_p =M( 1,1,0,alpha,beta,CV,CA);
  //GF^2 m^9 c^9/36/\[Pi]^5/hbar^10 - CGS coeff returns Q [erg/cm^3/sec] (like Itoh)
  cgs_coeff = 9.279177945324861E18;
  Q_half_avg = (M_10_m+M_10_p)*cgs_coeff;
  
  deltaQ = 4.0*Gm( 0.0,alpha,beta)*Gp(-0.5,alpha,beta) 
         -     Gm( 1.0,alpha,beta)*Gp(-0.5,alpha,beta) 
         - 4.0*Gm(-0.5,alpha,beta)*Gp( 0.0,alpha,beta) 
         + 7.0*Gm( 0.5,alpha,beta)*Gp( 0.0,alpha,beta) 
         - 7.0*Gm( 0.0,alpha,beta)*Gp( 0.5,alpha,beta) 
         + 4.0*Gm( 1.0,alpha,beta)*Gp( 0.5,alpha,beta) 
         +     Gm(-0.5,alpha,beta)*Gp( 1.0,alpha,beta) 
         - 4.0*Gm( 0.5,alpha,beta)*Gp( 1.0,alpha,beta)
  ;
  
  deltaQ = CV*CA*deltaQ*cgs_coeff;
  
  switch(flavour){
    case 1:
    return Q_half_avg-deltaQ;
    break;
    
    case 2:
    return Q_half_avg+deltaQ;
    break;
    
    case 3:
    return Q_half_avg-deltaQ;
    break;
    
    case 4:
    return Q_half_avg+deltaQ;
    break;
  }
  
  return 0.0;
}


/*
Reaction rate for pair-annihilation process [1/cm^3/sec]
FINAL VERSION. ONLY PHYSICAL CONSTANTS MAY CHANGE
IN FUTURE VERSIONS
*/
double pair_R(double kT, double mu, int flavour)
{ 

  const double sin2_theta_W=0.22280,m=0.510998903;
  double CV=1.0,CA=0.5;
  double alpha,beta;
  double M_00,cgscoeff;		  
  
  switch(flavour){
    case 1: CV =  0.5 + 2.0 * sin2_theta_W;
    break;

    case 2: CV =  0.5 + 2.0 * sin2_theta_W;
    break;
	
    case 3: CV = -0.5 + 2.0 * sin2_theta_W; CA=-0.5;
    break;

    case 4: CV = -0.5 + 2.0 * sin2_theta_W; CA=-0.5;
    break;
  };

  alpha=m/kT; beta=mu/kT;
  

  M_00 = M(1,0,0,alpha,beta,CV,CA);
  
  //GF^2 m^8 c^7/18/\[Pi]^5/hbar^10 - CGS coeff returns R [1/cm^3/sec] (like Schinder)
  cgscoeff = 2.2668740500597664E25;
  return cgscoeff*M_00;

}

double pair_J(double kT, double mu, int flavour)
{
  const double sin2_theta_W=0.22280,m=0.511;
  double alpha=m/kT,beta=mu/kT;  
  double CV=1.0,CA=0.5;
  double G_m_32,G_m_1,G_m_12,G_m_0,G_m_ng;
  double G_p_32,G_p_1,G_p_12,G_p_0,G_p_ng;
  double J2_avg, deltaJ2;
  
   switch(flavour){
    case 1: CV =  0.5 + 2.0 * sin2_theta_W;
    break;

    case 2: CV =  0.5 + 2.0 * sin2_theta_W;
    break;
	
    case 3: CV = -0.5 + 2.0 * sin2_theta_W; CA=-0.5;
    break;

    case 4: CV = -0.5 + 2.0 * sin2_theta_W; CA=-0.5;
    break;
  };

  
  G_m_32 = Gm( 1.5,alpha,beta);
  G_m_1  = Gm( 1.0,alpha,beta);
  G_m_12 = Gm( 0.5,alpha,beta);
  G_m_0  = Gm( 0.0,alpha,beta);
  G_m_ng = Gm(-0.5,alpha,beta);
  
  G_p_32 = Gp( 1.5,alpha,beta);
  G_p_1  = Gp( 1.0,alpha,beta);
  G_p_12 = Gp( 0.5,alpha,beta);
  G_p_0  = Gp( 0.0,alpha,beta);
  G_p_ng = Gp(-0.5,alpha,beta);
  
  
  J2_avg = (CV*CV+CA*CA)*(
                          28.0*(G_m_32*G_p_12+G_p_32*G_m_12)
                          +
                          30.0*G_p_1*G_m_1
                          -
                          7.0*(G_m_32*G_p_ng+G_p_32*G_m_ng)
                          -
                          8.0*(G_m_12*G_p_ng+G_p_12*G_m_ng)
                          )
                   +CV*CV*(
                          70.0*(G_m_12*G_p_ng+G_p_12*G_m_ng)
                          +
                          60.0*(G_m_1*G_p_0+G_p_1*G_m_0)  
                          )
       +30.0*(CV*CV-CA*CA)*G_m_0*G_p_0
       +16.0*(3.0*CV*CV-2.0*CA*CA)*G_m_12*G_p_12
       -(34.0*CV*CV-6.0*CA*CA)*G_m_ng*G_p_ng;
  
  deltaJ2  = CV*CA*(
                    4.0*(G_m_32*G_p_12-G_p_32*G_m_12)
                    -
                    (G_m_32*G_p_ng-G_p_32*G_m_ng)
                    +
                    6.0*(G_m_1*G_p_0-G_p_1*G_m_0)
                    +
                    4.0*(G_m_12*G_p_ng-G_p_12*G_m_ng)
                   );



//CGS: GF^2 m^10 c^11/hbar^10/\[Pi]^5/36 
  J2_avg   = 0.1*J2_avg*7.611409182236904E12;
  deltaJ2  = deltaJ2*7.611409182236904E12;

 switch(flavour){
    case 1:
    return J2_avg-deltaJ2;
    break;
    
    case 2:
    return J2_avg+deltaJ2;
    break;
    
    case 3:
    return J2_avg-deltaJ2;
    break;
    
    case 4:
    return J2_avg+deltaJ2;
    break;
  }
  
  
  return 0.0;
}


double integrand_plasmaL_Q(double k, double kT, double mu, int flavour)
{

  const double sin2_theta_W=0.22280;
  double CV=1.0,CA=0.5;
  double omega_l,Zl,nB,omega0,V,kmax_BS;

   switch(flavour){
    case 1: CV =  0.5 + 2.0 * sin2_theta_W;
    break;

    case 2: CV =  0.5 + 2.0 * sin2_theta_W;
    break;
	
    case 3: CV = -0.5 + 2.0 * sin2_theta_W; CA=-0.5;
    break;

    case 4: CV = -0.5 + 2.0 * sin2_theta_W; CA=-0.5;
    break;
  };
  
  omega0 = plasma_frequency(kT,mu);
  V = electron_velocity(kT,mu);
  kmax_BS = sqrt(3.0/V/V*(0.5/V*log((1.0+V)/(1.0-V))-1.0))*omega0;


  omega_l = dispersionL(k,kmax_BS,omega0,V,(kmax_BS-omega0)*1E-9);
  nB = 1.0/(exp(omega_l/kT)-1.0);
  Zl = 2.0*(omega_l*omega_l-V*V*k*k)/(3.0*omega0*omega0-omega_l*omega_l+V*V*k*k);


  return CV*CV*k*k*Zl*omega_l*omega_l*(omega_l*omega_l-k*k)*(omega_l*omega_l-k*k)*nB;;
}


double plasmaL_Q(double kT, double mu, int flavour)
{

  double integral, omega0, V, kmax_BS;
  omega0 = plasma_frequency(kT,mu);
  V = electron_velocity(kT,mu);
  kmax_BS = sqrt(3.0/V/V*(0.5/V*log((1.0+V)/(1.0-V))-1.0))*omega0;
  
  integral = gaussian_integrate(0.0, kmax_BS, &integrand_plasmaL_Q, kT, mu, flavour);
   
     
/*
Half of the Braaten&Segel emissivity gives
exactly neutrino/antineutrino emissivity as spectra are identical
CGS coefficient = GF^2 /hbar^10/c^9/Pi^4/alpha/96 ergtoMeV^9
*/  
  return 0.5*integral*6.317505411603926E23;
}


double integrand_plasmaL_R(double k, double kT, double mu, int flavour)
{

  const double sin2_theta_W=0.22280;
  double CV=1.0,CA=0.5;
  double omega_l,Zl,nB,omega0,V,kmax_BS;

   switch(flavour){
    case 1: CV =  0.5 + 2.0 * sin2_theta_W;
    break;

    case 2: CV =  0.5 + 2.0 * sin2_theta_W;
    break;
	
    case 3: CV = -0.5 + 2.0 * sin2_theta_W; CA=-0.5;
    break;

    case 4: CV = -0.5 + 2.0 * sin2_theta_W; CA=-0.5;
    break;
  };
  
  omega0 = plasma_frequency(kT,mu);
  V = electron_velocity(kT,mu);
  kmax_BS = sqrt(3.0/V/V*(0.5/V*log((1.0+V)/(1.0-V))-1.0))*omega0;


  omega_l = dispersionL(k,kmax_BS,omega0,V,(kmax_BS-omega0)*1E-9);
  nB = 1.0/(exp(omega_l/kT)-1.0);
  Zl = 2.0*(omega_l*omega_l-V*V*k*k)/(3.0*omega0*omega0-omega_l*omega_l+V*V*k*k);


  return CV*CV*k*k*Zl*omega_l*(omega_l*omega_l-k*k)*(omega_l*omega_l-k*k)*nB;
               
}

double plasmaL_R(double kT, double mu, int flavour)
{ 
 
  double integral, omega0, V, kmax_BS;
  omega0 = plasma_frequency(kT,mu);
  V = electron_velocity(kT,mu);
  kmax_BS = sqrt(3.0/V/V*(0.5/V*log((1.0+V)/(1.0-V))-1.0))*omega0;
  
  integral = gaussian_integrate(0.0, kmax_BS, &integrand_plasmaL_R, kT, mu, flavour);
   
   
// GF^2 /hbar^10/c^9/\[Pi]^4/\[Alpha]/96 ergtoMeV^8     
  return integral*3.943045809823163E29;
}

double integrand_plasmaT_Q(double k, double kT, double mu, int flavour)
{ 
  double omega_t,Zt,nB,omega0,V,mt_BS,mt;
  const double sin2_theta_W=0.22280;
  double CV=1.0,CA=0.5;

  
  
  switch(flavour){
    case 1: CV =  0.5 + 2.0 * sin2_theta_W;
    break;

    case 2: CV =  0.5 + 2.0 * sin2_theta_W;
    break;
	
    case 3: CV = -0.5 + 2.0 * sin2_theta_W; CA=-0.5;
    break;

    case 4: CV = -0.5 + 2.0 * sin2_theta_W; CA=-0.5;
    break;
  };
  
  omega0 = plasma_frequency(kT,mu);
  V = electron_velocity(kT,mu);
  mt_BS = sqrt(3.0/2.0/V/V*(1.0-(1.0-V*V)/2.0/V*log((1.0+V)/(1.0-V))))*omega0;
  mt = photon_mass(kT,mu);

   omega_t = dispersionT(k,mt,omega0,V,12.0*kT*1E-9);
    nB = 1.0/(exp(omega_t/kT)-1.0);
    Zt = 2.0*omega_t*omega_t*(omega_t*omega_t-V*V*k*k)/(3.0*omega0*omega0*omega_t*omega_t+(omega_t*omega_t+k*k)*(omega_t*omega_t-V*V*k*k)-2.0*omega_t*omega_t*(omega_t*omega_t-k*k));
    
   return CV*CV*k*k*Zt*(omega_t*omega_t-k*k)*
   (omega_t*omega_t-k*k)*(omega_t*omega_t-k*k)*
   nB;



}


double plasmaT_Q(double kT, double mu, int flavour)
{ 
  double integral;

  integral = gaussian_integrate(0.0, 12.0*kT, &integrand_plasmaT_Q, kT, mu, flavour);
    
//    GF^2 /hbar^10/c^9/\[Pi]^4/\[Alpha]/96   ergtoMeV^9 
  return integral*6.317505411603926E23;
}


double integrand_plasmaT_R(double k, double kT, double mu, int flavour)
{


  double omega_t,Zt,nB,omega0,V,mt_BS,mt;
  const double sin2_theta_W=0.22280;
  double CV=1.0,CA=0.5;
  double integrand;
    
  switch(flavour){
    case 1: CV =  0.5 + 2.0 * sin2_theta_W;
    break;

    case 2: CV =  0.5 + 2.0 * sin2_theta_W;
    break;
        
    case 3: CV = -0.5 + 2.0 * sin2_theta_W; CA=-0.5;
    break;

    case 4: CV = -0.5 + 2.0 * sin2_theta_W; CA=-0.5;
    break;
  };
  
  omega0 = plasma_frequency(kT,mu);
  V = electron_velocity(kT,mu);
  mt_BS = sqrt(3.0/2.0/V/V*(1.0-(1.0-V*V)/2.0/V*log((1.0+V)/(1.0-V))))*omega0;
  //mt = photon_mass(kT,mu);
  mt = omega0;

  
    omega_t = dispersionT(k,mt,omega0,V,12.0*kT*1E-9);
    nB = 1.0/(exp(omega_t/kT)-1.0);
    Zt = 2.0*omega_t*omega_t*(omega_t*omega_t-V*V*k*k)/(3.0*omega0*omega0*omega_t*omega_t+(omega_t*omega_t+k*k)*(omega_t*omega_t-V*V*k*k)-2.0*omega_t*omega_t*(omega_t*omega_t-k*k));
    integrand  = k*k*Zt*(omega_t*omega_t-k*k)*
              (omega_t*omega_t-k*k)*(omega_t*omega_t-k*k)*nB/omega_t;
 
    return CV*CV*integrand;  

}


double plasmaT_R(double kT, double mu, int flavour)
{

  double integral;

  integral = gaussian_integrate(0.0   ,  6.0*kT, &integrand_plasmaT_R, kT, mu, flavour)+
             gaussian_integrate(6.0*kT, 12.0*kT, &integrand_plasmaT_R, kT, mu, flavour);


// 2 GF^2 /hbar^10/c^9/\[Pi]^4/\[Alpha]/96 ergtoMeV^8 
return integral*7.886091619646326E29;

}