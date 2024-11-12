#include <math.h>
#include <stdio.h>

#ifndef _PSNS_PHOTO_C
#define _PSNS_PHOTO_C


static inline double Sqr(double x)
    {
    return x*x;
    };

//------------------------------------------------------------------------------------------------------------------


double M2_photo_Salpeter(
		double E1,double p1x,double p1y,double p1z,
		double E2,double p2x,double p2y,double p2z,
		double Enu1,double q1x,double q1y,double q1z,
		double Enu2,double q2x,double q2y,double q2z,
		double omega,double kx,double ky,double kz,
		double m)
{
#if 0
Incoming electron 4-momentum	(E1,p1x,p1y,p1z)
Outcoming electron 4-momentum	(E2,p2x,p2y,p2z)
Neutrino 4-momentum	(Enu1,q1x,q1y,q1z)
Anti-Neutrino 4-mom	(Enu2,q2x,q2y,q2z)
Plasmon  4-momentum	(omega,kx,ky,kz)
G - Fermi coonstant
e - electron charge
CV,CA - standard model couplings
m - electron mass
#endif


double M2,P1_P2,P1_Q1,P2_Q2,P1_Q2,P2_Q1,Q1_K,Q2_K,P1_K,P2_K;



//4-momenta scalar products
P1_P2 = E1*E2 - p1x*p2x - p1y*p2y - p1z*p2z;



P1_Q1 = E1*Enu1 - p1x*q1x - p1y*q1y - p1z*q1z;

P1_Q2 = E1*Enu2 - p1x*q2x - p1y*q2y - p1z*q2z;

P2_Q1 = E2*Enu1 - p2x*q1x - p2y*q1y - p2z*q1z;

P2_Q2 = E2*Enu2 - p2x*q2x - p2y*q2y - p2z*q2z;


P1_K =  E1  *omega - kx*p1x - ky*p1y - kz*p1z;

P2_K  = E2  *omega - kx*p2x - ky*p2y - kz*p2z;

Q1_K =  Enu1*omega - kx*q1x - ky*q1y - kz*q1z;

Q2_K  = Enu2*omega - kx*q2x - ky*q2y - kz*q2z;



//Squared matrix element expressed by four-momenta
//dot products.
M2 =

128.0*(

m*m*(P1_Q1*Q2_K - P1_Q1*P2_Q2)/Sqr(P2_K) 
- 
m*m*(P2_Q2*Q1_K + P1_Q1*P2_Q2)/Sqr(P1_K) 
+ 
(P1_Q1*P1_Q2 + Q1_K*P2_Q2 + P1_Q1*P2_Q2)/P1_K 
+ 
(P1_Q1*Q2_K - P1_Q1*P2_Q2 - P2_Q1*P2_Q2)/P2_K
+ 
P1_P2*( Q1_K*P2_Q2-P1_Q1*Q2_K +  2.0*P1_Q1*P2_Q2 )/(P1_K*P2_K) 

);

 

	return M2;
}

//------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------
#if 0
Differential cross section for photoproduction of the
neutrino - antineutrino pair:
d/dcosthetaQQ/dEnu1/dEnu2
NOTE: this is only integrand, you must integrate
to get true cross-section!

Based on:
S. Dutta, S. Ratkovic, M. Prakash,
,,Photoneutrino process in astrophysical systems''
Physical Review D, 69, 023005, (2004)
http://link.aps.org/abstract/PRD/v69/e023005
doi:10.1103/PhysRevD.69.023005


costhetaQQ - cosine of the angle between outcoming neutrino
and antineutrino

Enu1 - outcoming neutrino energy

Enu2 - outcoming antineutrino energy

Integration variables required to compute differential cross-section:
costhetaE=-1..1, phiE=0..2 Pi, phiK=0..2 Pi, 
omega=0..Infinity (Emax cutoff, incoming photon energy),
p=0..Infinity (Emax cutoff, outcoming electron momentum)

Thermodynamic properties:
kT - temperature required to compute Fermi-Dirac (FD) and Bose-Einstein
(BE) functions for photons
mu - chemical potential required to compute FD function for electrons



G - Fermi coonstant
e - electron charge
CV,CA - standard model couplings
m - electron mass
#endif
double integr_photo(
  double costhetaQQ, double Enu1, double Enu2, double costhetaE, double phiE,
  double phiK,double omega,double p, 
  double kT, double mu,  
  double CV, double CA, double m)
{
	
double Msquared;
//Temporary variables for 4-momenta
double _E1, _p1x, _p1y, _p1z; 
double _E2, _p2x, _p2y, _p2z;
double _Enu1, _q1x, _q1y, _q1z; 
double _Enu2, _q2x, _q2y, _q2z; 
double _omega, _kx, _ky, _kz;

//Invariants and components of the total incoming particles momenta
double PPmodul, MM2, EE;
double Px,Py,Pz;
double Px1,Py1,Pz1;
double Px2,Py2,Pz2;
	
//Angles between total incoming momentum and photon momentum 
double costhetaK,sinthetaK;
	
//Fermi-Dirac and Bose-Eistein  
double FD1, FD2, BE;

double sin_phi_P, cos_phi_P, sin_theta_P, cos_theta_P;
double  _kx1, _ky1, _kz1, _kx2, _ky2, _kz2;
	
double sinthetaE, sinthetaQQ;
	
  sinthetaE  = sqrt(1.0-costhetaE*costhetaE);
  sinthetaQQ = sqrt(1.0-costhetaQQ*costhetaQQ);
		
  _E2  = sqrt(p*p+m*m);
  _p2x = p*sinthetaE*cos(phiE);
  _p2y = p*sinthetaE*sin(phiE);
  _p2z = p*costhetaE;
	
  _Enu1 = Enu1;
  _q1x = 0.0;
  _q1y = 0.0;
  _q1z = Enu1;

 _Enu2 = Enu2;
  _q2x  = Enu2*sinthetaQQ;
  _q2y  = 0.0;
  _q2z  = Enu2*costhetaQQ;

/*
FOUR-MOMENTA of OUTCOMING PARTICLES now are expressed by integration variables
*/

	
  EE = _E2 + _Enu1 + _Enu2;
  Px = _p2x + _q1x + _q2x;
  Py = _p2y + _q1y + _q2y;
  Pz = _p2z + _q1z + _q2z;
	
  MM2 = Sqr(_E2 + _Enu1 + _Enu2) - Sqr(_p2x + _q1x + _q2x) 
       -Sqr(_p2y + _q1y + _q2y ) - Sqr(_p2z + _q1z + _q2z);
	
	
  PPmodul = sqrt(EE*EE - MM2);
	

	
  costhetaK = (m*m - MM2  + 2.0 * EE * omega)/2.0/PPmodul/omega;

  if (fabs(costhetaK)>1.0) return 0.0;
	
  sinthetaK = sqrt(1.0 - costhetaK*costhetaK);


  _omega = omega;
  _kx = omega*sin(phiK)*sinthetaK;
  _ky = omega*cos(phiK)*sinthetaK;
  _kz = omega*costhetaK;



// OBROTY 

  sin_phi_P = Py/sqrt(Px*Px+Py*Py);
  cos_phi_P = Px/sqrt(Px*Px+Py*Py);

  sin_theta_P = sqrt(Px*Px+Py*Py)/sqrt(Px*Px+Py*Py+Pz*Pz);
  cos_theta_P = Pz/sqrt(Px*Px+Py*Py+Pz*Pz);



//TU MA BYC obrot wokol osi x o kat theta_P

	
  _kx1 =  _kx*cos_theta_P + _kz*sin_theta_P;
  _ky1 =  _ky;
  _kz1 =  -_kx*sin_theta_P + _kz*cos_theta_P;
	
  Px1 =  Px*cos_phi_P +  Py*sin_phi_P;
  Py1 = -Px*sin_phi_P +  Py*cos_phi_P;
  Pz1 =  Pz;



//A TU MA BYC obrot wokol osi z o kat phi_P- 3/2 Pi

	
  _kx2 =   _kx1*cos_phi_P - _ky1*sin_phi_P;
  _ky2 =   _kx1*sin_phi_P + _ky1*cos_phi_P;
  _kz2 =   _kz1;
	

  Px2 =   -Px1*sin_theta_P - Pz1*cos_theta_P;
  Py2 =   Py1;
  Pz2 =   Px1*cos_theta_P - Pz1*sin_theta_P;
	

//TAK wyglada wektor po rotacji o kat theta_P wokol osi x i phi_P-3/2 Pi wokol osi z

  _kx = _kx2;
  _ky = _ky2;
  _kz = _kz2;



/*
NOW PLASMON FOUR-MOMENTUM IS DETERMINED
*/

  _E1 = _E2 + _Enu1 + _Enu2 - _omega;
  _p1x = _p2x + _q1x +  _q2x  - _kx;
  _p1y = _p2y + _q1y +  _q2y  - _ky;
  _p1z = _p2z + _q1z +  _q2z  - _kz;

/*
ALL FOUR_MOMENTS are expressed by integration vars
*/

  Msquared = M2_photo_Salpeter( _E1, _p1x, _p1y, _p1z, _E2, _p2x, _p2y, _p2z, _Enu1, _q1x, _q1y, _q1z, _Enu2, _q2x, _q2y, _q2z, _omega, _kx, _ky, _kz, m);

//Fermi-Dirac distribution of the incoming electron
 FD1 = 1.0/(1.0+exp( (Enu1 + Enu2 - mu - omega +  sqrt(m*m + p*p) )/kT ));
	
//Fermi blocking factor for outcoming electron
  FD2 = 1.0 - 1.0 / (1.0 + exp( (sqrt(m*m + p*p)-mu)/kT ));
	
//Bose-Einstein distribution of the incoming photon
  BE = 1.0  / ( exp(omega/kT) -1.0 ) ;
	
  if(Msquared<0.0) {return 0.0;}; 

  return Msquared / PPmodul * p*p / sqrt(p*p+m*m) * Enu1 * Enu2 * FD1 * FD2 * BE;

	

}


#endif

