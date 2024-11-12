double pair(double, double, double, int, int);

double photo(double, double, double, int, int);

double plasmaT(double, double, double, int);

double plasmaL(double, double, double, int);

double integr_photo( double, double, double, double, double, double, double, double, double, double, double, double, double, double, double);

double integr_pair1( double,  double, double, double, double, double, double, double, double);

double integr_pair2( double,  double, double, double, double, double, double, double, double);

double integrDicus( double, double, double, double, double, double, double, double);

double M(int, int, int, double, double, double, double);

double pair_Q(double, double, int);
double pair_R(double, double, int);
double pair_J(double, double, int);

/* versions using pre-calculated tabulated results */
double pair_tbl(double, double, double, int);
double pair_tbl_nue(double, double, double);
double pair_tbl_nuebar(double, double, double);
double pair_tbl_numu(double, double, double);

double pair_Q_tbl(double, double, int);
double pair_R_tbl(double, double, int);
//double pair_J_tbl(double, double, int);

double rhoYe(double, double);
double chemical_potential(double, double, double);
double chemical_potential_tbl(double, double);



double plasma_frequency(double, double);
double photon_mass(double, double);
//omega1/omega0
double electron_velocity(double, double);
double max_momentum(double, double);


//Transverse dispersion NOT YET IMPLEMENTED !!!!!!!
double dispersionL(double, double, double, double, double);
double dispersionT(double, double, double, double, double );
double inv_dispersionL(double, double, double, double );
double inv_dispersionT(double, double, double, double);


double plasmaL_Q(double, double, int);
double plasmaL_R(double, double, int);
double plasmaT_Q(double, double, int);
double plasmaT_R(double, double, int);

/* versions using pre-calculated tabulated results */
double plasmaL_Q_tbl(double, double, int);
double plasmaL_R_tbl(double, double, int);
double plasmaT_Q_tbl(double, double, int);
double plasmaT_R_tbl(double, double, int);


/*
* $Id: thermal.h,v 1.7 2007/02/17 12:24:43 cvsroot Exp $
* $Log: thermal.h,v $
* Revision 1.7  2007/02/17 12:24:43  cvsroot
* Second pair spectrum moment pair_J and plasmon decay  plasmaL_Q etc. added
*
* Revision 1.6  2007/02/11 11:36:42  cvsroot
* New routines for standard F-D integrals, M function, Q, R pair
*
* Revision 1.5  2007/02/09 10:25:38  cvsroot
* Number of arguments for dispesion relations corrected. Some comments updated
*
* Revision 1.4  2007/02/03 12:17:26  cvsroot
* dispersion relations included but not yet implemented
*
* $Date: 2007/02/17 12:24:43 $
* $Author: cvsroot $
* $Revision: 1.7 $
*/
