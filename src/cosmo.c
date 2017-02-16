#include "globals.h"

static double curvature_parameter(const double z);
double a2t_cgs(double a);

static const double a0 = 1;

void Set_cosmology()
{
    Cosmo.a2t = &a2t_cgs;
    Cosmo.h_100 = 0.7; // model parameters at z = 0
    Cosmo.Omega_M = 0.3;
    Cosmo.Omega_L = 0.7;

    Cosmo.Omega_0 = Cosmo.Omega_M + Cosmo.Omega_L;
    Cosmo.K = curvature_parameter( Cosmo.Omega_0 );

    Cosmo.H0_cgs = 100 * Cosmo.h_100 * 1e5 / 1000 /kpc2cgs; // cgs

    Cosmo.Rho_crit0 = 3./8./pi/Grav * p2(Cosmo.H0_cgs);

    printf("System at:   z = %g \n"
            "   H/100       = %g\n"
            "   Omega_M     = %g\n"
            "   rho_crit(0) = %g g/cm^3\n"
            "   rho_crit(z) = %g g/cm^3\n"
            "   mean mol. w.= %g\n"
            "   E(z)        = %g\n"
            "   Delta       = %g\n\n",
            Param.Redshift, Cosmo.h_100, Cosmo.Omega_M,
            Cosmo.Rho_crit0, Critical_Density(Param.Redshift),
            mean_mol_weight, Ez(Param.Redshift),
            Overdensity_Parameter());

    return;
}

double Omega_M(const double z)
{
    return Cosmo.Omega_M * p3(1+z) / p2(Ez(z));
}

double Critical_Density(const double z)
{
    return 3 * p2(Hubble_Parameter(z))/(8*pi*Grav);
}

double curvature_parameter(const double Omega_0)
{
    if (Omega_0 == 0)
        return 0;

    return (Omega_0) / fabs(Omega_0);
}


/* Mo, v.d.Bosch, White (3.75) */
double Hubble_Parameter(const double z)
{
    return Cosmo.H0_cgs*Ez(z); // in cgs
}

/* Mo, v.d.Bosch, White (2.62, 3.75), Boehringer 2012 (5) */
double Ez(const double z)
{
    return sqrt(Cosmo.Omega_L + (1-Cosmo.Omega_0) * p2(1+z)
            + Cosmo.Omega_M * p3(1+z));
}

/* Pierpaoli 2001, Boehringer+ 2012 */
const double cij[5][5] = {  // Pierpaoli+ 01 Table 1
    {546.67, -137.82, 94.083, -204.68,  111.51},
    {-1745.6, 627.22,  -1175.2, 2445.7,  -1341.7},
    {3928.8, -1519.3, 4015.8,  -8415.3, 4642.1},
    {-4384.8, 1748.7,  -5362.1, 11257.,  -6218.2},
    {1842.3, -765.53, 2507.7, -5210.7, 2867.5} };

double Overdensity_Parameter()
{
    const double x = Cosmo.Omega_M - 0.2;
    const double y = Cosmo.Omega_L;

    double result = 0;

    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
            result += cij[i][j] * pow(x,i) * pow(y,j);

    return Cosmo.Omega_M * result;
}

/* Concordance Cosmology */
double a2t_cgs(double a)
{
    const double km2cgs = 1e5;
    const double H0 = 100 * km2cgs / kpc2cgs/1000 * Cosmo.h_100;

    double t = 2.0/3.0 / (sqrt(Cosmo.Omega_M)*H0) // Mo+, eq (3.89)
        * asinh( pow( a/a0 * pow( Cosmo.Omega_L/Cosmo.Omega_M, 1.0/3.0), 3.0/2.0));

    return t;
}

double t2a_cgs(double t)
{
    const double km2cgs = 1e5;
    const double H0 = 100 * km2cgs / kpc2cgs/1000 * Cosmo.h_100;

    double a = a0 * pow( Cosmo.Omega_M/Cosmo.Omega_L, 1.0/3.0 )
        * pow( sinh( 1.5 * sqrt(Cosmo.Omega_L)*H0*t), 2.0/3.0 );

    return a;
}
