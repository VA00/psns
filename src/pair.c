#include <math.h>
#include <cuba.h>

#ifndef _PSNS_PAIR_C
#define _PSNS_PAIR_C


static inline double Sqr(double x)
    {
    return x*x;
    };

//------------------------------------------------------------------------------------------------------------------

static inline double UnitStep( double x)
    {
    if ( x >= 0.0 ) return 1.0;
    return 0.0;
    }

//------------------------------------------------------------------------------------------------------------------

double H1( double E1, double E2, double k0, double ct, double m)
    {	
    double h1;

    h1= ((-4*(E1*E2*sqrt((E1*E1) - (m*m)) + (m*m)*(sqrt((E1*E1) - (m*m)) + ct*sqrt((E2*E2) - (m*m))) - 
        E1*(k0*sqrt((E1*E1) - (m*m)) + ct*E1*sqrt((E2*E2) - (m*m))))*
        (-((ct*ct)*(E2*E2)*sqrt((E1*E1) - (m*m))) + 
        E2*(E1 - k0)*(sqrt((E1*E1) - (m*m)) + ct*sqrt((E2*E2) - (m*m))) + 
        (m*m)*(sqrt((E1*E1) - (m*m)) + (ct*ct)*sqrt((E1*E1) - (m*m)) + 
        2*ct*sqrt((E2*E2) - (m*m))) - 
        E1*(k0*sqrt((E1*E1) - (m*m)) + ct*E1*sqrt((E2*E2) - (m*m)) + ct*k0*sqrt((E2*E2) - (m*m)))))/
        ((E1*E1) + (E2*E2) - 2*(m*m) + 2*ct*sqrt((E1*E1) - (m*m))*sqrt((E2*E2) - (m*m))) + 
        2*((1 + (ct*ct))*(m*m*m*m) + (E2*E2)*((E1*E1) + (ct*ct)*((E1*E1) - (m*m))) - 
        2*E1*E2*(E1*k0 - (m*m) + ct*sqrt((E1*E1) - (m*m))*sqrt((E2*E2) - (m*m))) + 
        E1*k0*(E1*k0 + 2*ct*sqrt((E1*E1) - (m*m))*sqrt((E2*E2) - (m*m))) - 
        (m*m)*((ct*ct)*(E1*E1) + 2*E1*k0 + 2*ct*sqrt((E1*E1) - (m*m))*sqrt((E2*E2) - (m*m)))) + 
        (((E1*E1) - (m*m))*(3*Sqr((ct*ct)*(E2*E2)*sqrt((E1*E1) - (m*m)) + 
        E2*(-E1 + k0)*(sqrt((E1*E1) - (m*m)) + ct*sqrt((E2*E2) - (m*m))) - 
        (m*m)*(sqrt((E1*E1) - (m*m)) + (ct*ct)*sqrt((E1*E1) - (m*m)) + 
        2*ct*sqrt((E2*E2) - (m*m))) + 
        E1*(k0*sqrt((E1*E1) - (m*m)) + ct*E1*sqrt((E2*E2) - (m*m)) + ct*k0*sqrt((E2*E2) - (m*m)))
        ) - (-(E1*E1) - (E2*E2) + 2*(m*m) - 
        2*ct*sqrt((E1*E1) - (m*m))*sqrt((E2*E2) - (m*m)))*
        ((k0*k0)*(-(E1*E1) + (-1 + (ct*ct))*(m*m)) + 
        2*E1*k0*((m*m) - ct*sqrt((E1*E1) - (m*m))*sqrt((E2*E2) - (m*m))) - 
        2*E2*(-E1 + k0)*(E1*k0 - (m*m) + ct*sqrt((E1*E1) - (m*m))*sqrt((E2*E2) - (m*m))) - 
        (E2*E2)*(E1*(E1 - 2*k0) + (ct*ct)*((E1*E1) + (k0*k0) - (m*m))) - 
        (m*m)*((1 + (ct*ct))*(m*m) - 
        ct*(ct*(E1*E1) + 2*sqrt((E1*E1) - (m*m))*sqrt((E2*E2) - (m*m)))))))/
        Sqr((E1*E1) + (E2*E2) - 2*(m*m) + 2*ct*sqrt((E1*E1) - (m*m))*sqrt((E2*E2) - (m*m))))/
	sqrt((E1*E1) + (E2*E2) - 2*(m*m) + 2*ct*sqrt((E1*E1) - (m*m))*sqrt((E2*E2) - (m*m)));
    
   return h1;
    
   };

//------------------------------------------------------------------------------------------------------------------

double H2( double E1, double E2, double k0, double ct, double m)
    {
    double h2;
    
    h2= (2*(E1*E1)*(k0*k0) - (4*E1*k0*sqrt((E1*E1) - (m*m))*
	((ct*ct)*(E2*E2)*sqrt((E1*E1) - (m*m)) + 
	E2*(-E1 + k0)*(sqrt((E1*E1) - (m*m)) + ct*sqrt((E2*E2) - (m*m))) - 
        (m*m)*(sqrt((E1*E1) - (m*m)) + (ct*ct)*sqrt((E1*E1) - (m*m)) + 
        2*ct*sqrt((E2*E2) - (m*m))) + 
        E1*(k0*sqrt((E1*E1) - (m*m)) + ct*E1*sqrt((E2*E2) - (m*m)) + ct*k0*sqrt((E2*E2) - (m*m)))))/
        ((E1*E1) + (E2*E2) - 2*(m*m) + 2*ct*sqrt((E1*E1) - (m*m))*sqrt((E2*E2) - (m*m))) + 
        (((E1*E1) - (m*m))*(3*Sqr((ct*ct)*(E2*E2)*sqrt((E1*E1) - (m*m)) + 
        E2*(-E1 + k0)*(sqrt((E1*E1) - (m*m)) + ct*sqrt((E2*E2) - (m*m))) - 
        (m*m)*(sqrt((E1*E1) - (m*m)) + (ct*ct)*sqrt((E1*E1) - (m*m)) + 
        2*ct*sqrt((E2*E2) - (m*m))) + 
        E1*(k0*sqrt((E1*E1) - (m*m)) + ct*E1*sqrt((E2*E2) - (m*m)) + ct*k0*sqrt((E2*E2) - (m*m)))
        ) - (-(E1*E1) - (E2*E2) + 2*(m*m) - 
        2*ct*sqrt((E1*E1) - (m*m))*sqrt((E2*E2) - (m*m)))*
        ((k0*k0)*(-(E1*E1) + (-1 + (ct*ct))*(m*m)) + 
        2*E1*k0*((m*m) - ct*sqrt((E1*E1) - (m*m))*sqrt((E2*E2) - (m*m))) - 
        2*E2*(-E1 + k0)*(E1*k0 - (m*m) + ct*sqrt((E1*E1) - (m*m))*sqrt((E2*E2) - (m*m))) - 
        (E2*E2)*(E1*(E1 - 2*k0) + (ct*ct)*((E1*E1) + (k0*k0) - (m*m))) - 
        (m*m)*((1 + (ct*ct))*(m*m) - 
        ct*(ct*(E1*E1) + 2*sqrt((E1*E1) - (m*m))*sqrt((E2*E2) - (m*m)))))))/
        Sqr((E1*E1) + (E2*E2) - 2*(m*m) + 2*ct*sqrt((E1*E1) - (m*m))*sqrt((E2*E2) - (m*m))))/
	sqrt((E1*E1) + (E2*E2) - 2*(m*m) + 2*ct*sqrt((E1*E1) - (m*m))*sqrt((E2*E2) - (m*m)));
   
   return h2;
   
   };

//------------------------------------------------------------------------------------------------------------------
   
double H3( double E1, double E2, double k0, double ct, double m)
   {
   double h3;
   
    h3= (E1*E2 + m*m - ct*sqrt(E1*E1 - m*m)*sqrt(E2*E2 - m*m) )/
	sqrt(E1*E1 + E2*E2 - 2*m*m 
	+ 2*ct*sqrt(E1*E1 - m*m)*sqrt(E2*E2 - m*m));
   
   return h3;   
   };

//------------------------------------------------------------------------------------------------------------------
   
double ov1(double E1, double E2, double k0, double ct, double m, double CV, double CA){

double h1,h2,h3,OV1;    
    
h1 = H1(E1, E2, k0, ct, m);
h2 = H2(E1, E2, k0, ct, m);    
h3 = H3(E1, E2, k0, ct, m);
        
OV1 = (Sqr(CV - CA)*h1 + Sqr(CV + CA)*h2 + 2.0*(CV*CV-CA*CA)*h3*m*m)
		/(8.0*E1*E2*3.14159265);
return OV1;
};

//------------------------------------------------------------------------------------------------------------------

double ov2( double E1, double E2, double k0, double ct, double m, double CV, double CA){
double h1,h2,h3,OV2;

h1 = H1(E1, E2, k0, ct, m);
h2 = H2(E1, E2, k0, ct, m);    
h3 = H3(E1, E2, k0, ct, m);

OV2 = (Sqr(CV + CA)*h1 + Sqr(CV - CA)*h2 + 2.0*(CV*CV-CA*CA)*h3*m*m)
		/(8.0*E1*E2*3.14159265);
return OV2;
};

//------------------------------------------------------------------------------------------------------------------

double ovDicus( double E1, double E2, double ct, double m, double CV, double CA){

double OVDICUS;
    
OVDICUS = ((3*(CV*CV-CA*CA)*(m*m*m*m + 
        	m*m*(E1*E2 - ct*sqrt(E1*E1-m*m)*sqrt(E2*E2-m*m))) + 
    		(Sqr(CA) + Sqr(CV))*((m*m*m*m) + 3*(m*m)*
        	(E1*E2 - ct*sqrt((E1*E1) - (m*m))*sqrt((E2*E2) - (m*m))) + 
        	2*Sqr(E1*E2 - ct*sqrt((E1*E1) - (m*m))*sqrt((E2*E2) - (m*m))))))/(12.*E1*E2*3.14159265);

return OVDICUS;
};

//------------------------------------------------------------------------------------------------------------------

double integr_pair1(double E1, double E2, double k0, double ct, double m, double CV, double CA, double kT, double mu){

double Theta, integrand1, OV1;
	
Theta = UnitStep(
-	
(k0 - (E1 + E2)/2.0  
- sqrt(E1*E1 + E2*E2 - 2.0*m*m + 2.0*ct*sqrt(E1*E1 - m*m)*sqrt(E2*E2 - m*m))/2.0)
*
(k0 - (E1 + E2)/2.0  
+ sqrt(E1*E1 + E2*E2 - 2.0*m*m + 2.0*ct*sqrt(E1*E1 - m*m)*sqrt(E2*E2 - m*m))/2.0)

);

if (Theta == 0.0) return 0.0;

OV1 = ov1( E1, E2, k0, ct, m, CV, CA);
     
integrand1=8.0*3.1415926*3.1415926*(E1*E2*sqrt(E1*E1 - m*m)*sqrt(E2*E2 - m*m)*OV1)
		/(1.0 +  exp((E1 - mu)/kT)  )/(1.0 + exp((E2 + mu)/kT));

return integrand1*Theta;
};

//------------------------------------------------------------------------------------------------------------------

double integr_pair2(double E1, double E2, double k0, double ct, double m, double CV, double CA, double kT, double mu)
    {
    double  Theta, integrand2, OV2, FD1, FD2;



    
Theta = UnitStep(

-
(k0 - (E1 + E2)/2.0  
- sqrt(E1*E1 + E2*E2 - 2.0*m*m + 2*ct*sqrt(E1*E1 - m*m)*sqrt(E2*E2 - m*m))/2.0)
*
(k0 - (E1 + E2)/2.0  
+ sqrt(E1*E1 + E2*E2 - 2.0*m*m + 2*ct*sqrt(E1*E1 - m*m)*sqrt(E2*E2 - m*m))/2.0)

);
   


//This save some computing time
if (Theta == 0.0) return 0.0;	    


    OV2 = ov2( E1, E2, k0, ct, m, CV, CA);


    integrand2 = 
8.0 * 3.1415926 *3.1415926 *
(E1*E2*sqrt(E1*E1 - m*m)*sqrt(E2*E2 - m*m)*OV2)
/
 (
( 1.0 + exp((E1 - mu)/kT) )
*
( 1.0 + exp((E2 + mu)/kT) )
)
;
//return -1.0;    
     return integrand2;
     
    };

//-------------------------------------------------------------------

double integrDicus(double E1, double E2,  double ct, double m, double CV, double CA, double kT, double mu)
    {
    double integrandDicus, OVDICUS;

    OVDICUS = ovDicus( E1, E2, ct, m, CV, CA);

    integrandDicus = 	(E1*E2*sqrt((E1*E1) - (m*m))*sqrt((E2*E2) - (m*m))*OVDICUS)/((1.0 + exp((E1 - mu)/kT))*(1.0 + exp((E2 + mu)/kT)));
    
    return integrandDicus;
    } 

//-------------------------------------------------------------------
#endif

