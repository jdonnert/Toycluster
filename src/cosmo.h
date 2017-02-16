/* cosmology */
extern struct Universe{
    // double Baryon_Fraction;  // Moved to Halo
    double h_100;
    double Omega_M;
    double Omega_L;
    double Omega_0;
    double K;
    double H0_cgs;
    double Overdensity_Parameter;    // Delta200 at cluster redshift
    double Rho_crit0;
    double (*a2t)(double);
} Cosmo;

double Overdensity_Parameter();
double Hubble_Parameter(const double z) ;
double Ez(const double z);
double Critical_Density(const double z);
double Omega_M(const double z);
