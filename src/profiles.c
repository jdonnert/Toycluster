#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include "globals.h"

#define NTABLE 256
#define NSAMPLE (32*NTABLE) // oversample !
#define GSL_SPLINE gsl_interp_cspline

static void setup_internal_energy_profile(const int i);
static void setup_gas_potential_profile(const int i);
static void    setup_dm_potential_profile(const int i);
static void setup_dm_distribution_function(const int iCluster);
static double hernquist_distribution_func(const int iCluster, const double E);

static void show_profiles(const int iCluster);

static struct int_param {
    double rho0;
    double beta;
    double rc;
    double rcut;
    double cuspy;
} ip;

/* This file contains all physical profiles that describe the cluster model.
 * All profiles follow the same scheme: define an integrant, make a table of
 * solutions using a GSL integrator, spline interpolate the table for speed
 * and accuracy. If the spline interpolation fails, its usually, because the
 * table contains data that is not strictly increasing in x. We include
 * analytical profiles for testing as well, however numerical stability is the
 * ultimate test of the model.
 * To setup the profiles for a given halo, Setup_Profiles has to be called in
 * the parent routine. */

void Setup_Profiles(const int i)
{
    Setup_DM_Mass_Profile(i);

    setup_dm_potential_profile(i);

    if ((Halo[i].bf > 0) && (!Halo[i].Is_Stripped)) {

        Setup_Gas_Mass_Profile(i);

        setup_gas_potential_profile(i);

        setup_internal_energy_profile(i);
    }

    setup_dm_distribution_function(i);

    show_profiles(i);

    return ;
}

/* As the NFW does not converge, we have to cut it off around R_Sample */

double DM_Density_Profile(const int i, const float r)
{
    double ra = r/Halo[i].Rs;

    double rho_nfw = Halo[i].Rho0_DM / (ra * p2(1+ra));

    const double rmax = Halo[i].R_Sample[1];

    return rho_nfw / (1+ p3(r/rmax)); // with cutoff

    //return Hernquist_Density_Profile(i, r);
}

/* Hernquist 1989, eq. 2 , eq. 17-19, Binney & Tremaine 2.64 */
double Hernquist_Density_Profile(const int i, const double r)
{
    double ra = r / Halo[i].A_hernq;
    double a = Halo[i].A_hernq;
    double Mdm = Halo[i].Mass[1] * Halo[i].MassCorrFac;

    return Halo[i].Rho0_DM / (ra * p3(1 + ra));
    //return Mdm * a/r/p3(a+r)/2/pi;
}

/* DM cumulative mass profile from an arbitrary DM density profile. */

static gsl_spline *DMMinv_Spline = NULL;
static gsl_interp_accel *DMMinv_Acc = NULL;
#pragma omp threadprivate(DMMinv_Spline, DMMinv_Acc)

static gsl_spline *DMM_Spline = NULL;
static gsl_interp_accel *DMM_Acc = NULL;
#pragma omp threadprivate(DMM_Spline, DMM_Acc)

double DM_Mass_Profile(const double r, const int i)
{
    double Mr = gsl_spline_eval(DMM_Spline, r, DMM_Acc);

    return Mr;
}

double Inverted_DM_Mass_Profile(double q, const int i)
{
    double m = q * Halo[i].Mass[1];

    return gsl_spline_eval(DMMinv_Spline, m, DMMinv_Acc);
}

static double dm_mr_integrant(double r, void * param)
{
    int i = *((int *) param);

    return 4.0*pi*r*r * DM_Density_Profile(i, r);
}

void Setup_DM_Mass_Profile(const int iCluster)
{
    double m_table[NTABLE] = { 0 };
    double r_table[NTABLE] = { 0 };

    double rmin = Zero;
    double rmax = Infinity;
    double log_dr = ( log10(rmax/rmin) ) / (NTABLE - 1);

    gsl_function gsl_F = { 0 };

    gsl_integration_workspace *gsl_workspace = NULL;
    gsl_workspace = gsl_integration_workspace_alloc(NTABLE);

    for (int i = 1; i < NTABLE; i++) {

        double error = 0;

        r_table[i] = rmin * pow(10, log_dr * i);

        gsl_F.function = &dm_mr_integrant;
        gsl_F.params = (void *) &iCluster;

        gsl_integration_qag(&gsl_F, 0, r_table[i], 0, 1e-5, NTABLE,
                GSL_INTEG_GAUSS61, gsl_workspace, &m_table[i], &error);
    }

    r_table[0] = m_table[0] = 0;

    for (int i = 1; i < NTABLE; i++)
        if (m_table[i-1] >= m_table[i]) // make strictly increasing
            m_table[i] = m_table[i-1] * (1 + 1e-4);

    #pragma omp parallel
    {

    DMM_Acc  = gsl_interp_accel_alloc();

    DMM_Spline = gsl_spline_alloc(GSL_SPLINE, NTABLE);
    gsl_spline_init(DMM_Spline, r_table, m_table, NTABLE);

    DMMinv_Acc  = gsl_interp_accel_alloc();

    DMMinv_Spline = gsl_spline_alloc(GSL_SPLINE, NTABLE);
    gsl_spline_init(DMMinv_Spline, m_table, r_table, NTABLE);

    } // omp parallel

    return ;
}

double DM_Mass_Profile_NFW(const double r, const int i)
{
    double rs = Halo[i].Rs;
    double rho0 = Halo[i].Rho0_DM; // (Wiki, Mo+ 2010)

    return 4*pi*rho0*p3(rs) * (log((rs+r)/rs) - r/(rs+r));
}

double DM_Mass_Profile_HQ(const double r, const int i)
{
    double ra = r/Halo[i].A_hernq;

    return Halo[i].Rho0_DM * 4*pi/2 * p3(Halo[i].A_hernq) * p2(ra)/p2(1+ra);
}

/* We need the precise potential up to very large distances. Because the
 * numerical integrator runs into precision issues for large r, we extra-
 * polate linearly the potential in log space beyond 1 Gpc */

static gsl_spline *DMPsi_Spline = NULL;
static gsl_interp_accel *DMPsi_Acc = NULL;
#pragma omp threadprivate(DMPsi_Spline, DMPsi_Acc)

static double dm_psi_integrant(const double r, void * param)
{
    int i = *((int *) param);

    return G/p2(r) * DM_Mass_Profile(r, i);
}

static void    setup_dm_potential_profile(const int iCluster)
{
    double psi_table[NTABLE] = { 0 };
    double r_table[NTABLE] = { 0 };

    double rmin = Zero;
    double rmax = Infinity;
    double log_dr = ( log10(rmax/rmin) ) / (NTABLE - 1);

    gsl_function gsl_F = { 0 };

    gsl_integration_workspace *gsl_workspace = NULL;
    gsl_workspace = gsl_integration_workspace_alloc(NTABLE);

    for (int i = 1; i < NTABLE; i++) {

        double error = 0;

        r_table[i] = rmin * pow(10, log_dr * i);

        int ip = iCluster;

        gsl_F.function = &dm_psi_integrant;
        gsl_F.params = (void *) &ip;

        gsl_integration_qag(&gsl_F, Zero, r_table[i], 0, 1e-5, NTABLE,
                GSL_INTEG_GAUSS61, gsl_workspace, &psi_table[i], &error);
    }
    psi_table[0] = psi_table[1];

    for (int i = 0; i < NTABLE; i++)
        psi_table[i] -= psi_table[NTABLE-1]; // == 0 @ infinity

    double dpsi_dr = 0, psi0 = 0;

    for (int i = 1; i < NTABLE-1; i++) { // extrapolate powerlaw

        if (r_table[i] < 2*Halo[iCluster].R_Sample[1])
            continue;

        if (dpsi_dr == 0) {

            dpsi_dr = (log10(-psi_table[i+1]) - log10(-psi_table[i-1]))
                    / (log10(   r_table[i+1]) - log10(   r_table[i-1]));

            psi0 =log10(-psi_table[i]) - dpsi_dr* log10(r_table[i]);
        }

        psi_table[i] = -pow(10, dpsi_dr*log10(r_table[i]) + psi0);
    }

    psi_table[NTABLE-1] = 0;

    #pragma omp parallel
    {

    DMPsi_Acc  = gsl_interp_accel_alloc();

    DMPsi_Spline = gsl_spline_alloc(GSL_SPLINE, NTABLE);
    gsl_spline_init(DMPsi_Spline, r_table, psi_table, NTABLE);

    } // omp parallel

    return ;
}

double DM_Potential_Profile(const int i, const float r)
{
    return gsl_spline_eval(DMPsi_Spline, r, DMPsi_Acc);
}

double DM_Potential_Profile_NFW(const int i, const double r)
{
    double rrs = r/Halo[i].Rs; // B&T 2.67

    return -G *  Halo[i].Rho0_DM * 4*pi*p2(Halo[i].Rs) * log(1+rrs)/(rrs);
}

double DM_Potential_Profile_HQ(const int i, const double r)
{
    double ra = r/Halo[i].A_hernq;

    return -G * Halo[i].Rho0_DM * 4*pi*p2(Halo[i].A_hernq) /2/ (1+ra);
}

/* Double beta profile at rc and rcut */

double Gas_Density_Profile(const double r, const int i)
{
    double rho0 = Halo[i].Rho0;
    double beta = Halo[i].Beta;
    double rc = Halo[i].Rcore;
    double rcut = Halo[i].Rcut;
    double Is_Cupsy = Halo[i].Have_Cuspy;

    double rho = rho0 * pow(1 + p2(r/rc), -3.0/2.0*beta)
                / (1 + p3(r/rcut));

#ifdef DOUBLE_BETA_COOL_CORES

    double rho0_cc = rho0 * Param.Rho0_Fac;
    double rc_cc = rc / Param.Rc_Fac;

    if (Is_Cuspy)
        rho += rho0_cc / (1 + p2(r/rc_cc)) / (1 + p3(r/rcut) * (r/rcut));

#endif // DOUBLE_BETA_COOL_CORES

    return rho;
}

/* Cumulative Mass profile from the gas density profile by numerical
 * integration and cubic spline interpolation. This has to be called once
 * before the mass profile of a halo is used, to set the spline variables. */

static gsl_spline *M_Spline = NULL;
static gsl_interp_accel *M_Acc = NULL;
#pragma omp threadprivate(M_Spline, M_Acc)

static gsl_spline *Minv_Spline = NULL;
static gsl_interp_accel *Minv_Acc = NULL;
#pragma omp threadprivate(Minv_Spline, Minv_Acc)

double Gas_Mass_Profile(const double r_in, const int i)
{
    double r = fmin(r_in, Halo[i].R_Sample[0]);

    return  gsl_spline_eval(M_Spline, r, M_Acc);
}

double Inverted_Gas_Mass_Profile(double M)
{
    return gsl_spline_eval(Minv_Spline, M, Minv_Acc);
}

double m_integrant(double r, void * param)
{
    int i = *((int *) param);

    return 4*pi*r*r*Gas_Density_Profile(r, i);
}

void Setup_Gas_Mass_Profile(const int j)
{
    double m_table[NTABLE] = { 0 };
    double r_table[NTABLE] = { 0 };
    double dr_table[NTABLE] = { 0 };

    double rmin = Zero;

    double rmax = Halo[j].R_Sample[0]*1.2; // include R_Sample

    double log_dr = ( log10(rmax/rmin) ) / (NTABLE - 1);

    gsl_function gsl_F = { 0 };

    gsl_integration_workspace *gsl_workspace = NULL;
    gsl_workspace = gsl_integration_workspace_alloc(NTABLE);

    for (int i = 1; i < NTABLE; i++) {

        double error = 0;

        r_table[i] = rmin * pow(10, log_dr * i);

        gsl_F.function = &m_integrant;
        gsl_F.params = (void *) &j;

        gsl_integration_qag(&gsl_F, 0, r_table[i], 0, 1e-5, NTABLE,
                GSL_INTEG_GAUSS61, gsl_workspace, &m_table[i], &error);

        if (m_table[i] < m_table[i-1])
            m_table[i] = m_table[i-1]; // integrator may fluctuate

        //printf("%g \n", m_table[i]);
    }

    m_table[0] = 0;

    #pragma omp parallel
    {

    M_Acc  = gsl_interp_accel_alloc();

    M_Spline = gsl_spline_alloc(GSL_SPLINE, NTABLE);
    gsl_spline_init(M_Spline, r_table, m_table, NTABLE);

    Minv_Acc  = gsl_interp_accel_alloc();

    Minv_Spline = gsl_spline_alloc(GSL_SPLINE, NTABLE);
    gsl_spline_init(Minv_Spline, m_table, r_table, NTABLE);

    } // omp parallel

    return ;
}

/* M(<= R) of a beta profile with beta=2/3 */

double Mass_Profile_23(const double r, const int i)
{
    double rho0 = Halo[i].Rho0;
    double rc = Halo[i].Rcore;
    double beta = Halo[i].Beta;
    double rcut = Halo[i].Rcut;
    int Is_Cuspy = Halo[i].Have_Cuspy;

    const double r2 = p2(r);
    const double rc2 = p2(rc);
    const double rcut2 = p2(rcut);

    double Mr = rho0 * rc2*rcut2*rcut/(8*(p2(rcut2)+p2(rc2))) * // fourth order
        ( sqrt2 *( (rc2-rcut2) *( log(rcut2 - sqrt2*rcut*r+r2)
                                - log(rcut2 + sqrt2*rcut*r+r2))
                   - 2 * (rc2 + rcut2) * atan(1 - sqrt2*r/rcut)
                   + 2 * (rc2 + rcut2) * atan(sqrt2 * r/rcut + 1))
          - 8 * rc * rcut * atan(r/rc));

#ifdef DOUBLE_BETA_COOL_CORES

    double rc_cc = rc / Param.Rc_Fac;
    double rc2_cc = p2(rc_cc);
    double rho0_cc = rho0 * Param.Rho0_Fac;

    double Mr_cc = 0;

    if (Is_Cuspy)
        Mr += rho0_cc * rc2_cc*rcut2*rcut/(8*(p2(rcut2)+p2(rc2_cc))) *
            ( sqrt2 *( (rc2-rcut2) *( log(rcut2 - sqrt2*rcut*r+r2)
                                    - log(rcut2 + sqrt2*rcut*r+r2))
                       - 2 * (rc2_cc + rcut2) * atan(1 - sqrt2*r/rcut)
                       + 2 * (rc2_cc + rcut2) * atan(sqrt2 * r/rcut + 1))
              - 8 * rc_cc * rcut * atan(r/rc));

#endif // DOUBLE_BETA_COOL_CORES

    return 4 * pi * Mr;
}

/* The grav. potential from the gas density. We solve Poisson's equation
 * numerically. Here psi is the -phi >= 0,
 * The potential is extrapolated away from the sampling radius for accuracy.
 * Because of the pole at r=0, we integrate from 1e-7 to help the quadrature
 * alg.*/

static gsl_spline *Psi_Spline = NULL;
static gsl_interp_accel *Psi_Acc = NULL;
#pragma omp threadprivate(Psi_Spline, Psi_Acc)

double Gas_Potential_Profile(const int i, const double r)
{
    double r_max = Halo[i].R_Sample[0];

    if (r < r_max)
        return gsl_spline_eval(Psi_Spline, r, Psi_Acc);

    double psi_max = gsl_spline_eval(Psi_Spline,r_max, Psi_Acc);

    return psi_max*r_max/r; // extrapolate, point mass
}

double psi_integrant(double r, void *param)
{
    int i = *((int *) param);

    if (r == 0)
        return 0;

    return G / (r*r) * Gas_Mass_Profile(r, i);
}

static void setup_gas_potential_profile(const int i)
{
    double error = 0;

    double psi_table[NTABLE] = { 0 };
    double r_table[NTABLE] = { 0 };

    double rmin = Zero;
    double rmax = Infinity;
    double log_dr = ( log10(rmax/rmin) ) / (NTABLE - 1);

    gsl_function gsl_F = { 0 };
    gsl_F.function = &psi_integrant;
    gsl_F.params = (void *) &i;

    gsl_integration_workspace *gsl_workspace = NULL;
    gsl_workspace = gsl_integration_workspace_alloc(NSAMPLE);

    for (int j = 0; j < NTABLE; j++) {

        r_table[j] = rmin * pow(10, log_dr * j);

        gsl_integration_qag(&gsl_F, 1e-7, r_table[j], 0, 1e-3, NSAMPLE,
              GSL_INTEG_GAUSS61, gsl_workspace, &psi_table[j], &error);
    }

    double gauge = 0;
    gsl_integration_qag(&gsl_F, 1e-7, Infinity, 0, 1e-3, NSAMPLE,
          GSL_INTEG_GAUSS61, gsl_workspace, &gauge, &error);

    for (int j = 0; j < NTABLE; j++) // psi = -phi > 0
        psi_table[j] = -1*(psi_table[j] - gauge);

    r_table[0] = 0;
    psi_table[0] = gauge; // make sure its correct

    #pragma omp parallel
    {

    Psi_Acc  = gsl_interp_accel_alloc();

    Psi_Spline = gsl_spline_alloc(GSL_SPLINE, NTABLE);
    gsl_spline_init(Psi_Spline, r_table, psi_table, NTABLE);

    } // omp parallel

    return ;
}

double Gas_Potential_Profile_23(const int i, const float r) // beta = 2/3
{
    if (r > 2*Halo[i].Rcut)
        return 0;

    double rc = Halo[i].Rcore;
    const double rcut = Halo[i].Rcut;

    const double r2 = r*r;
    double rc2 = rc*rc;
    const double rcut2 = rcut*rcut;

    double psi = -1 * rc2*rcut2/(8*(rc2*rc2 + rcut2*rcut2)*r)
                *(8*rc*rcut2*atan(r/rc) + 4*rc2*r*atan(r2/rcut2)
                + rcut*(2*sqrt2*(rc2 + rcut2)*atan(1 - (sqrt2*r)/rcut)
                    - 2*sqrt2*(rc2 + rcut2)* atan(1 + (sqrt2*r)/rcut)
                    + 4*rcut*r*log(rc2 + r2)
                    - sqrt2*rc2* log(rcut2 - sqrt2*rcut*r + r2)
                    + sqrt2*rcut2*log(rcut2 - sqrt2*rcut*r + r2)
                    + sqrt2*rc2*log(rcut2 + sqrt2*rcut*r + r2)
                    - sqrt2*rcut2*log(rcut2 + sqrt2*rcut*r + r2)
                    - 2*rcut*r*log(rcut2*rcut2 + r2*r2))) ;

    psi *= Halo[i].Rho0;

#ifdef DOUBLE_BETA_COOL_CORES

    rc /= rc / Param.Rc_Fac;
    rc2 = p2(rc);

    double psi_cc = -1 * rc2*rcut2/(8*(rc2*rc2 + rcut2*rcut2)*r)
                *(8*rc*rcut2*atan(r/rc) + 4*rc2*r*atan(r2/rcut2)
                + rcut*(2*sqrt2*(rc2 + rcut2)*atan(1 - (sqrt2*r)/rcut)
                    - 2*sqrt2*(rc2 + rcut2)* atan(1 + (sqrt2*r)/rcut)
                    + 4*rcut*r*log(rc2 + r2)
                    - sqrt2*rc2* log(rcut2 - sqrt2*rcut*r + r2)
                    + sqrt2*rcut2*log(rcut2 - sqrt2*rcut*r + r2)
                    + sqrt2*rc2*log(rcut2 + sqrt2*rcut*r + r2)
                    - sqrt2*rcut2*log(rcut2 + sqrt2*rcut*r + r2)
                    - 2*rcut*r*log(rcut2*rcut2 + r2*r2))) ;

    psi += psi_cc * Halo[i].Rho0 * Param.Rho0_Fac;

#endif // DOUBLE_BETA_COOL_CORES

    //double psi =  p3(rc)*( 0.5/rc*log(rc*rc+r*r) + atan(r/rc)/r )
//            * p2(p2(1 - r/rcut)) // no cutoff

    return 4*pi*G * psi;
}


/* Numerical solution for all kinds of gas densities. We spline interpolate
 * a table of solutions for speed. We have to solve the hydrostatic equation
 * (eq. 9 in Donnert 2014). */

static gsl_spline *U_Spline = NULL;
static gsl_interp_accel *U_Acc = NULL;
#pragma omp threadprivate(U_Spline, U_Acc)

double Internal_Energy_Profile(const int i, const double r)
{
    return gsl_spline_eval(U_Spline, r, U_Acc);;
}

static double u_integrant(double r, void *param) // Donnert 2014, eq. 9
{
    const int i = *((int *) param);

    double rho_gas = Gas_Density_Profile(r, i);
    double Mr_Gas = Gas_Mass_Profile(r, i);
    double Mr_DM = DM_Mass_Profile(r,i);

    return rho_gas /(r*r) * (Mr_Gas + Mr_DM);
}

static void setup_internal_energy_profile(const int i)
{
    double u_table[NTABLE] = { 0 },
           r_table[NTABLE] = { 0 };

    double rmin = Zero;
    double rmax = Infinity;
    double dr = ( log10(rmax/rmin) ) / (NTABLE-1);

#ifdef NO_RCUT_IN_T
    double rcut_old = Halo[i].Rcut;
    Halo[i].Rcut = 10 * Halo[i].R200;
#endif

    #pragma omp parallel
    {

    gsl_function gsl_F = { 0 };

    gsl_integration_workspace *gsl_workspace = NULL;
    gsl_workspace = gsl_integration_workspace_alloc(NTABLE);

    #pragma omp for
    for (int j = 1; j < NTABLE;  j++) {

        double error = 0;

        double r = rmin * pow(10, dr * j);

        r_table[j] = r;

        gsl_F.function = &u_integrant;
        gsl_F.params = (void *) &i;

        gsl_integration_qag(&gsl_F, r, rmax, 0, 1e-5, NTABLE,
                GSL_INTEG_GAUSS61, gsl_workspace, &u_table[j], &error);

        double rho_gas = Gas_Density_Profile(r, i);

        u_table[j] *= G/((adiabatic_index-1)*rho_gas); // Donnert 2014, eq. 9
    }

    gsl_integration_workspace_free(gsl_workspace);

    } // omp parallel

    u_table[0] = u_table[1];

    #pragma omp parallel
    {

    U_Acc = gsl_interp_accel_alloc();

    U_Spline = gsl_spline_alloc(GSL_SPLINE, NTABLE);
    gsl_spline_init(U_Spline, r_table, u_table, NTABLE);

    } // omp parallel

#ifdef NO_RCUT_IN_T
    Halo[i].Rcut = rcut_old;
#endif

    return ;
}

/* Standard analytical temperature profile from Donnert 2014.
 * To avoid negative temperatures we define rmax*sqrt3 as outer radius */

static double F1(const double r, const double rc, const double a)
{
    const double rc2 = rc*rc;
    const double a2 = a*a;

    double result = (a2-rc2)*atan(r/rc) - rc*(a2+rc2)/(a+r)
                + a*rc * log( (a+r)*(a+r) / (rc2 + r*r) );

    result *= rc / p2(a2 + rc2);

    return result;
}

static double F2(const double r, const double rc)
{
    return p2(atan(r/rc)) / (2*rc) + atan(r/rc)/r;
}

double Internal_Energy_Profile_Analytic(const int i, const double d)
{
    double rho0 = Halo[i].Rho0;
    double a = Halo[i].A_hernq;
    double rc = Halo[i].Rcore;
    double rmax = Param.Boxsize; // "open" T boundary
    double Mdm = Halo[i].Mass[1];

    double u = G / ( (adiabatic_index-1) ) * ( 1 + p2(d/rc) ) *
                ( Mdm * (F1(rmax, rc, a) - F1(d, rc, a))
                + 4*pi*rho0*p3(rc) * (F2(rmax, rc) - F2(d, rc) ) );
    return u;
}

/* Find f(E) for arbitrary spherically symmetric density-potential pairs by
 * numerical integration of the Eddington equation.
 * We interpolate rho(psi) with a cubic spline using the GSL library and
 * get the second derivative from the spline directly. This is a hard 1D
 * integral to get to floating point precision ! The Hernquist
 * f(E) is reproduced with a few 1e-3 accuracy ...
 * Kazantzidis+ 2004, Binney 1982, Binney & Tremaine pp. 298, Barnes 02 */

typedef struct {
    double E;
    gsl_spline *spline;
    gsl_interp_accel *acc;
} interpolation_parameters;

static interpolation_parameters fE_params;
#pragma omp threadprivate(fE_params)

double Distribution_Function(const double E)
{
    double fE =  gsl_spline_eval(fE_params.spline, E, fE_params.acc);

    return fE;
}

/* This is the absolute potential psi = -phi.  */
double Potential_Profile(const int i, const double r)
{
    double psi = fabs(DM_Potential_Profile(i, r)); // DM generated potential

    if (Halo[i].Npart[0]) // Gas generated potential
        psi += fabs(Gas_Potential_Profile(i,r));

    return psi;
}

/* Binney & Tremaine sect. 4.3.1, we take the second derivate from the
 * cubic spline */

static double eddington_integrant(double psi, void *params)
{
    const interpolation_parameters *p = params;

    if (psi >= p->E)
        return 0;

    const double dRhodPsi2 = gsl_spline_eval_deriv2(p->spline, psi, p->acc);

    double result = dRhodPsi2 / sqrt(p->E - psi);

    return result;
}

static void setup_dm_distribution_function(const int iCluster)
{
    const double rmin = Zero;
    const double rmax = Infinity;

    double rstep = log10(rmax/rmin) / NSAMPLE; // good up to one Gpc

    double psi[NSAMPLE] = { 0 },
           rho[NSAMPLE] = { 0 };

    for (int i = 0; i < NSAMPLE; i++) {

        double r = rmin * pow(10, i*rstep);

        rho[i] = DM_Density_Profile(iCluster, r);

        psi[i] = Potential_Profile(iCluster, r); // psi = - phi

        //printf("%d %g %g %g \n", i, r, rho[i], psi[i]);
    }

    psi[NSAMPLE-1] = rho[NSAMPLE-1] = 0; // make sure  E==0 is covered

    double E[NTABLE] = { 0 },
           fE[NTABLE] = { 0 };

    double x[NSAMPLE] = { 0 }, // holds inverted vals
           y[NSAMPLE] = { 0 };

    for (int i = 0; i < NSAMPLE; i++) { // invert profile for interpolation

        x[i] = psi[NSAMPLE-i-1];
        y[i] = rho[NSAMPLE-i-1];

        //printf("%d %g %g \n", i, x[i], y[i]);
    }

    #pragma omp parallel // make f(E) table
    {
    interpolation_parameters int_params;
    int_params.acc = gsl_interp_accel_alloc();
    int_params.spline = gsl_spline_alloc(gsl_interp_akima, NSAMPLE);

    gsl_spline_init(int_params.spline, x, y, NSAMPLE);

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(NSAMPLE);

    rstep = log10(rmax/rmin) / NTABLE;

    const double sqrt8 = sqrt(8.0);

    double err = 0;

    #pragma omp for
    for (int i = 0; i < NTABLE; i++) {

        double r = rmin * pow(10, i*rstep);

        E[i] = Potential_Profile(iCluster, r);

        gsl_function F = {&eddington_integrant, &int_params};
        int_params.E = E[i];

        double Emax = E[i];
        double Emin = Emax / 1024; // ~= 0, all of the integrant is close to Emax

        gsl_integration_qags(&F, Emin, Emax, 0, 1e-3, NSAMPLE, w, &fE[i], &err);

        fE[i] /= sqrt8 * pi * pi;

        //printf("%d %g %g %g %g %g %g \n", i, r, E[i], fE[i],
        //        hernquist_distribution_func(iCluster, E[i]), err/fE[i],
        //        eddington_integrant(0.5*E[i], &int_params));
    }

    E[NTABLE-1] = 0; // r == inf, f(E) -> 0
    fE[NTABLE-1] = 0;

    gsl_integration_workspace_free(w);

    gsl_interp_accel_free(int_params.acc);
    gsl_spline_free(int_params.spline);

    } // omp parallel

    #pragma omp parallel
    {
    fE_params.acc = gsl_interp_accel_alloc();
    fE_params.spline = gsl_spline_alloc(gsl_interp_cspline, NTABLE);
    }

    int nSpline = 0;

    for (int i = 0; i < NTABLE; i++) { // interpolate fE in log space

        x[i] = (E[NTABLE-i-1]);

        y[i] = (fE[NTABLE-i-1]);

        //printf("%d %g %g \n", i, x[i], y[i]);
    }

    #pragma omp parallel
    gsl_spline_init(fE_params.spline, x, y, NTABLE);

    /*for (int i = 0; i < NTABLE; i++) {

        double r = rmin * pow(10, i*rstep);

        double E = Potential_Profile(iCluster, r);

        printf("%d %g %g ", i, r, E);

        double solution =  hernquist_distribution_func(iCluster, E);
        double fe = distribution_function(E);

        printf("%g %g %g \n", solution, fe, fabs(fe-solution)/solution);
    } */

    return ;
}

static double hernquist_distribution_func(const int iCluster, const double E)
{
    const double a = Halo[iCluster].A_hernq;
    const double mass = Halo[iCluster].Mass[1];
    const double prefac = 1/(sqrt2 * p3(2*pi) * pow(G*mass*a, 3./2.));

    double q2 = a*E/(G*mass);

    double f_E = prefac * mass * sqrt(q2) / pow(1-q2,2)
                *( (1 - 2*q2) * (8*q2*q2 - 8*q2 - 3)
                 + (3*asin(sqrt(q2)))/sqrt(q2*(1-q2) ));
    return f_E;
}
static int iLast = -1;

static void show_profiles(const int iCluster)
{
    if (iCluster < iLast)
        return ;

    iLast = iCluster;

    char fname[128] = {""};

    sprintf(fname, "profiles_%03d.txt", iCluster);

    FILE *fp = fopen(fname, "w");

    fprintf(fp, "#r, rho_dm, mr_dm, psi_dm, rho_HQ,  mr_nfw, mr_hq,"
                " psi_nfw, psi_hq ");
    if (Halo[iCluster].Mass200[0] > 0)
        fprintf(fp, "rho_gas, mr_gas, psi_gas, u_gas, mr_23, psi_23, u_ana ");
    fprintf(fp, "E, fE, fE_hq \n");

    const double rmin = Zero;
    const double rmax = Infinity;

    double dr = ( log10(rmax/rmin) ) / (NTABLE-1);

    for (int i = 0; i < NTABLE; i++) {

        double r = rmin * pow(10, i*dr);

        double rho_dm = DM_Density_Profile(iCluster, r);
        double mr_dm = DM_Mass_Profile(r, iCluster);
        double psi_dm  = DM_Potential_Profile(iCluster,r);

        double rho_HQ = Hernquist_Density_Profile(iCluster, r);

        double mr_nfw = DM_Mass_Profile_NFW(r, iCluster);
        double mr_hq = DM_Mass_Profile_HQ(r, iCluster);
        double psi_nfw = DM_Potential_Profile_NFW(iCluster, r);
        double psi_hq = DM_Potential_Profile_HQ(iCluster, r);

        fprintf(fp,"%g %g %g %g %g %g %g %g %g ",
            r, rho_dm, mr_dm, psi_dm, rho_HQ,  mr_nfw, mr_hq, psi_nfw, psi_hq);

        if (Halo[iCluster].Mass200[0] > 0) {

            double rho_gas = Gas_Density_Profile(r, iCluster);
            double mr_gas = Gas_Mass_Profile(r, iCluster);
            double psi_gas = Gas_Potential_Profile(iCluster, r);
            double u_gas = Internal_Energy_Profile(iCluster, r);

            double mr_23 = Mass_Profile_23(r, iCluster);
            double psi_23 = Gas_Potential_Profile_23(iCluster, r);
            double u_ana = Internal_Energy_Profile_Analytic(iCluster, r);

            fprintf(fp,"%g %g %g %g %g %g %g ",
                rho_gas, mr_gas, psi_gas, u_gas, mr_23, psi_23, u_ana);
        }

        double E = Potential_Profile(iCluster, r);
        double fE = Distribution_Function(E);
        double fE_HQ = hernquist_distribution_func(iCluster, E);

        fprintf(fp, "%g %g %g ", E, fE, fE_HQ);

        fprintf(fp,"\n");
    }

    fclose(fp);

    return ;
}
