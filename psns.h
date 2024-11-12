double pair(double, double, double, int, int);

double photo(double, double, double, int, int);

double plasmaT(double, double, double, int);

double plasmaL(double, double, double, int);


double pair_Q(double, double, int);
double pair_R(double, double, int);
double pair_J(double, double, int);

/* versions using pre-calculated tabulated results */
double pair_tbl(double, double, double, int);
double pair_Q_tbl(double, double, int);
double pair_R_tbl(double, double, int);
//double pair_J_tbl(double, double, int); NOT IMPLEMENTED !

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


double ffn_rate_nue(double, double, int, int);
double ffn_rate_nuebar(double, double, int, int);

double weak_R(double, double, int, int);
double weak_Q(double, double, int, int);

double neutron(double, double, double, int); 
double neutron_R(double, double); 
double neutron_Q(double, double); 
double proton(double, double, double, int);
double proton_R(double, double);
double proton_Q(double, double);
double capture(double, double, double, double); 
double decay(double, double, double, double);
double NSE_neutrino_spectrum(double, double, double, double, int); 
double NSE_neutrino_spectrum_tbl(double, double, double, double, int); 
void   NSE_neutrino_spectrum_table(double*, double, double, double, double, int, double, double ,int, int); 
