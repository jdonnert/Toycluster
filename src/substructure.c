#ifdef SUBSTRUCTURE

#include <gsl/gsl_sf_gamma.h>

#include "globals.h"

#define MIN_SUBHALO_MASS (10 * DESNNGB * (Param.Mpart[0] + Param.Mpart[1]))
#define MIN_DENSITY_CONTRAST 3

static void set_subhalo_masses(const double);
static void set_subhalo_particle_numbers();
static void set_subhalo_positions(const int);
static void set_subhalo_properties(const int);
static void    set_subhalo_pointers();
static void set_subhalo_bulkvel(const int i);
static bool reject_subhalo(const int);

static double subhalo_mass_fraction();
static double subhalo_mass_function(const double);
static double subhalo_number_density_profile(const double);
static double inverted_subhalo_number_density_profile(const double);
static double sampling_radius(const int, const double);
static double tidal_radius(const int, const double);
static double nfw_scale_radius(const double, const double, const double);
static double nfw_mass_profile(const double, const double, const double);

/*
 * Subhalos are treated as independent Halos inside r200 of SUBHOST
 */

void Setup_Substructure()
{
    Sub.First = 1;

    if (Param.Mass_Ratio != 0)
        Sub.First = 2;

    Sub.MassFraction = subhalo_mass_fraction();

    set_subhalo_masses(Sub.MassFraction);

    for (int i = Sub.First; i < Param.Nhalos; i++) {
	
        printf("."); fflush(stdout);

        do {

            set_subhalo_positions(i);

            set_subhalo_properties(i);

        } while (reject_subhalo(i));

#ifndef SLOW_SUBSTRUCTURE
        set_subhalo_bulkvel(i);
#endif
    }

    set_subhalo_particle_numbers();

    set_subhalo_pointers();

    printf("\nSubhalo Setup : \n"
        "   Total Mass DM   = %g \n"
        "   Mass Fraction   = %4.2g\n"
        "   Target Fraction = %g \n"
        "   Total Number    = %d / %d \n"
        "   Total Npart     = %d \n"
        "   Total Ngas      = %d \n"
        "   Total NDM       = %d \n",
        Sub.Mtotal, Sub.Mtotal / Halo[SUBHOST].Mtotal200, Sub.MassFraction,
        Sub.Nhalos, Param.Nhalos, Sub.Ntotal, Sub.Npart[0], Sub.Npart[1]);

#ifdef REPORTSUBHALOS
    for (int i = Sub.First; i < Param.Nhalos; i++)
        printf("Subhalo <%d> :\n"
                "   Npart         = %lld, %lld \n"
                "   Mass          = %g | %g %g \n"
                "   Mass200       = %g | %g %g \n"
                "   bf in rsample = %g \n"
                "   Mass Fraction = %g \n"
                "   DM  Mass      = %g \n"
                "   Gas Mass      = %g \n"
                "   c_nfw         = %g \n"
                "   rho0_DM         = %g \n"
                "   r_sample      = %g \n"
                "   R200          = %g \n"
                "   r_s           = %g \n"
                "   Hernquist a   = %g \n"
                "   core radius   = %g \n"
                "   rho0          = %g \n"
                "   MassCorrect.  = %g \n"
                "   x, y, z       = %g %g %g\n"
                "   vx,vy,vz      = %g %g %g\n",
                i, Halo[i].Npart[0], Halo[i].Npart[1],
                Halo[i].Mtotal,Halo[i].Mass[0], Halo[i].Mass[1],
                Halo[i].Mtotal200,Halo[i].Mass200[0], Halo[i].Mass200[1],
                Halo[i].Mass[0] / Halo[i].Mtotal,
                Halo[i].Mtotal200/Halo[0].Mtotal,
                Halo[i].Mass[1], Halo[i].Mass[0], Halo[i].C_nfw,
                Halo[i].Rho0_DM, Halo[i].R_Sample[0], Halo[i].R200,
                Halo[i].R200 / Halo[i].C_nfw,
                Halo[i].A_hernq, Halo[i].Rcore, Halo[i].Rho0,
                Halo[i].MassCorrFac,
                Halo[i].D_CoM[0],Halo[i].D_CoM[1],Halo[i].D_CoM[2],
                Halo[i].BulkVel[0],Halo[i].BulkVel[1],Halo[i].BulkVel[2]);
    fflush(stdout);
#endif

    return ;
}

/*
 * We sample the subhalo mass function from Giocoli 2010 MNRAS 408
 * this sets the DM mass in R_sample, not R200
 */

static void set_subhalo_masses(const double mass_fraction)
{
    const double mass_limit = Halo[SUBHOST].Mass200[1] * mass_fraction;

    const double qmax = subhalo_mass_function(MIN_SUBHALO_MASS)
        / MIN_SUBHALO_MASS;

    const double max_subhalo_mass = Sub.MassFraction * Halo[SUBHOST].Mass[1]/10;

    int i = Sub.First;

    while (Sub.Mtotal < mass_limit && (i < 70)) {

        double mDM = 0, q = 0;

        int j =0;

        for (j = 0; j < 10000; j++) { // rejection sampling

            mDM = MIN_SUBHALO_MASS + erand48(Omp.Seed) *
                (Halo[SUBHOST].Mass200[1] - MIN_SUBHALO_MASS);

            q = subhalo_mass_function(mDM) / mDM;

            double lower_bound = qmax * erand48(Omp.Seed);

            if (mass_limit - Sub.Mtotal < MIN_SUBHALO_MASS ){

                mDM = MIN_SUBHALO_MASS;

                break;
            }

            if (Sub.Mtotal + mDM > 1.05 * mass_limit)
                continue;

            if (mDM > max_subhalo_mass) // too big
                continue;

            if (q >= lower_bound)
                break;
        }

        if (j == 9999)
            mDM = MIN_SUBHALO_MASS;

#ifdef ADD_THIRD_SUBHALO
        if (i == Sub.First)
            mDM = Param.SubFirstMass;
#endif

        Halo[i].Mass[1] = mDM;

        Sub.Mtotal += Halo[i].Mass[1];

        Sub.Nhalos++;

        i++;

#ifdef THIRD_HALO_ONLY
        break;
#endif
    }

    Param.Nhalos += i-2;

    return ;
}

/* Sample mass dependent subhalo density profile (Gao+ 2002) */

static void set_subhalo_positions(int i)
{
    const double x_host = Halo[SUBHOST].D_CoM[0];
    const double y_host = Halo[SUBHOST].D_CoM[1];
    const double z_host = Halo[SUBHOST].D_CoM[2];

    double q = erand48(Omp.Seed);

	double x = 0, y = 0, z = 0;

    double r = 3 * Halo[SUBHOST].R200 * inverted_subhalo_number_density_profile(q);

	float theta = acos(2 *  erand48(Omp.Seed) - 1);
    float phi = 2*pi * erand48(Omp.Seed);

    x = r * sin(theta) * cos(phi);
    y = r * sin(theta) * sin(phi);
    z = r * cos(theta);

#ifdef ADD_THIRD_SUBHALO
    if (i == Sub.First) {

        x = Param.SubFirstPos[0] - x_host;
        y = Param.SubFirstPos[1] - y_host;
        z = Param.SubFirstPos[2] - z_host;
    }
#endif
    Halo[i].D_CoM[0] = (float) (x + x_host);
    Halo[i].D_CoM[1] = (float) (y + y_host);
    Halo[i].D_CoM[2] = (float) (z + z_host);

    return ;
}

/* here we impose restrictions on subhalos.
 * 1) they cannot overlap
 * 2) they have to be an substantial overdensity */

static bool reject_subhalo(const int i)
{
    bool resample = false;

    for (int j = Sub.First; j < i; j++) {

        double d[3] = { 0 };

        d[0] = Halo[i].D_CoM[0] - Halo[j].D_CoM[0];
        d[1] = Halo[i].D_CoM[1] - Halo[j].D_CoM[1];
        d[2] = Halo[i].D_CoM[2] - Halo[j].D_CoM[2];

        double r2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
        double size = Halo[i].R_Sample[0] + Halo[j].R_Sample[0];

        if (r2 < size*size)
            resample = true; // overlaps
    } // j

    double dx = Halo[i].D_CoM[0] - Halo[SUBHOST].D_CoM[0];
    double dy = Halo[i].D_CoM[1] - Halo[SUBHOST].D_CoM[1];
    double dz = Halo[i].D_CoM[2] - Halo[SUBHOST].D_CoM[2];

    double r = sqrt(dx*dx + dy*dy + dz*dz);

    double rho_host = DM_Density_Profile(SUBHOST, r);
    double rho_sub = DM_Density_Profile(i, 3 * Param.GravSofteningLength);

    if (rho_sub < rho_host*MIN_DENSITY_CONTRAST)
        resample = true;

    if (r > 3 * Halo[SUBHOST].R200)
        resample = true;

	if (SUBHOST == 1 && Halo[i].D_CoM[0] < 0)
        resample = true;

	if (SUBHOST == 0 && Halo[i].D_CoM[0] > 0)
        resample = true;

#ifdef ADD_THIRD_SUBHALO
    if (i == Sub.First)
        resample = false;
#endif

    return resample;
}

/*
 * fill Halo struct with profile parameters
 * note that we can determine the actual gas mass only at the
 * very end, because we need the mass in r200 + the mass profiles
*/

static void set_subhalo_properties(const int i)
{
	Halo[i].bf = 0.17;

    double dx = Halo[SUBHOST].D_CoM[0] - Halo[i].D_CoM[0];
    double dy = Halo[SUBHOST].D_CoM[1] - Halo[i].D_CoM[1];
    double dz = Halo[SUBHOST].D_CoM[2] - Halo[i].D_CoM[2];

    const double r_i = sqrt(dx*dx + dy*dy + dz*dz);

    double a =  Halo[SUBHOST].A_hernq / 10; // first guess (small)
    double r200 = Halo[SUBHOST].R200;
    double c_nfw = 0;
    double rsample = 0;

    int cnt = 0;

    for (;;) { // find a and R_Sample iteratively

        double last_a = a;

        rsample = fmax(sampling_radius(i, r_i), tidal_radius(i,r_i) );

        rsample = fmin(rsample, r200/4);

        Halo[i].C_nfw = c_nfw = Concentration_parameter(i);

        Halo[i].Rs = nfw_scale_radius(c_nfw, Halo[i].Mass[1], rsample);

        Halo[i].Rho0_DM = 1;
        Halo[i].Rho0_DM = Halo[i].Mass[1]
                                    / DM_Mass_Profile_NFW(rsample,i);

        a = Halo[i].Rs*sqrt(2*(log(1+c_nfw) - c_nfw/(1+c_nfw)));

        r200 = Halo[i].Rs * c_nfw;

#ifdef ADD_THIRD_SUBHALO
        rsample = r200;
#endif

        if ( fabs( (last_a - a)/a) < 1e-4)
            break;

        if (cnt++ > 100)
            break;
    }

    Halo[i].R_Sample[0] = Halo[i].R_Sample[1] = rsample;
    Halo[i].A_hernq = a;
    Halo[i].R200 = r200;
    Halo[i].C_nfw = c_nfw;

    const double r_strip = 0; //.1* Halo[0].R200;
                    // gas is assumed stripped beyond this radius

    Halo[i].Rcut = 0.6 * Halo[i].R_Sample[0];

    Halo[i].Mass200[1] = nfw_mass_profile(c_nfw, Halo[i].Rs, r200);

    if (r_i > r_strip)
        Halo[i].Mass200[0] = Halo[i].Mass200[1]
                            / (1/Halo[i].bf - 1);

#ifdef ADD_THIRD_SUBHALO
    if (i == Sub.First)
        Halo[i].Mass200[0] = Halo[i].Mass200[1]
            / (1/Cosmo.Baryon_Fraction - 1);
#endif

    Halo[i].Mtotal200 = Halo[i].Mass200[0] + Halo[i].Mass200[1];

    Halo[i].MassCorrFac = 1/(1 + 2*a/r200 + p2(a/r200));

    char tmp[CHARBUFSIZE] = {""};

    Halo[i].Beta = 2.0/3.0; // implicitely assumed

    double rc = Halo[i].Rcore = Gas_core_radius(i, tmp);

    Halo[i].Rho0 = Halo[i].Mass200[0]/(4*pi*p3(rc))/(r200/rc - atan(r200/rc));

    Halo[i].Mass[0] = 0;

    Halo[i].Is_Stripped  = true;

    if (r_i > r_strip) {

        Halo[i].Is_Stripped  = false;

        Halo[i].Mass[0] = Mass_Profile_23(Halo[i].R_Sample[0],i);
    }

#ifdef ADD_THIRD_SUBHALO
    if (i == Sub.First)
        Halo[i].Mass[0] = Gas_Mass_Profile(Halo[i].R_Sample[0], i);
#endif

    Halo[i].Mtotal = Halo[i].Mass[0] + Halo[i].Mass[1];

    return ;
}


static void set_subhalo_particle_numbers()
{
    const double mDM = Param.Mpart[1];
    const double mGas = Param.Mpart[0];

    for (int i = Sub.First; i < Param.Nhalos; i++) {

        int nDM = round(Halo[i].Mass[1] / mDM );
        int nGas = round(Halo[i].Mass[0] / mGas);

        if (mGas == 0) // DM only halo
            nGas = 0;

        int npart = nDM + nGas;

        Halo[i].Ntotal = npart;
        Halo[i].Npart[0] = nGas;
        Halo[i].Npart[1] = nDM;

        Sub.Ntotal += npart;
        Sub.Npart[0] += nGas;
        Sub.Npart[1] += nDM;

    }

    Halo[SUBHOST].Ntotal -= Sub.Ntotal; // correct <0> particle numbers
    Halo[SUBHOST].Npart[0] -= Sub.Npart[0];
    Halo[SUBHOST].Npart[1] -= Sub.Npart[1];

    return ;
}

static void    set_subhalo_pointers()
{
    int iGas = 0;
    int  iDM =  Param.Npart[0];

    for (int i = 0; i <SUBHOST+1; i++) {

        iGas += Halo[i].Npart[0];
        iDM += Halo[i].Npart[1];
    }

    for (int i = Sub.First; i < Param.Nhalos; i++) {

        Halo[i].Gas = &(P[iGas]);
        Halo[i].DM = &(P[iDM]);
        Halo[i].SphP = &(SphP[iGas]);

        iGas += Halo[i].Npart[0];
        iDM += Halo[i].Npart[1];
    }

    return ;
}

static double sampling_radius(const int i, const double d)
{
    const double rho_host = DM_Density_Profile(SUBHOST, d);

    const double m = Halo[i].Mass[1],
                   a = Halo[i].A_hernq;

    double left = 0,
           right = 10*Halo[0].R200,
           r = 0,
           delta = DBL_MAX;

    while (fabs(delta) > 1e-3) { // find root bisection

        r = left + 0.5 * (right - left);

        delta = (DM_Density_Profile(i, r) - rho_host)/rho_host;

        if (delta < 0)
            right = r;
        else
            left = r;
    }

    return r;
}

/* Tormen, Diaferio & Syer 1998, Hernquist 1990 */
static double tidal_radius(const int i, const double r)
{
    double m_sub = Halo[i].Mass[1];
    double m_host = Halo[SUBHOST].Mass200[1];
    double a = Halo[SUBHOST].A_hernq;

    double fac = (2*r*r / p2(a+r) * (1 - a*r*r / p3(r+a)));

    return r * pow(m_sub / (m_host * fac), 1.0/3.0);
}

/* Giocoli, Bartelmann et al. 2010, eq. 12, fig 2a */
static double subhalo_mass_function(const double m)
{
    const double cc = 1, Am = 9.33e-4, alpha = -0.9, beta = 12.2715;

    const double z = Param.Redshift;
    const double mSub = m * Unit.Mass / Msol2cgs;
    const double mHost = Halo[SUBHOST].Mass200[1] * Unit.Mass / Msol2cgs;

    const double x = mSub/mHost;

    return mHost * sqrt(1+z)*cc*Am* pow(mSub, alpha) * exp(-beta * p3(x));
}

/* Giocoli+ 2010 */
static double subhalo_mass_fraction()
{
#ifndef THIRD_HALO_ONLY
    return 0.22 * sqrt(1+Param.Redshift);
#else
    return Halo[SUBHOST].Mtotal200 / Param.SubFirstMass;
#endif
}

/* Gao+ 2004 */
static double subhalo_number_density_profile(const double r)
{
    const double  ac = 0.244 * Halo[SUBHOST].C_nfw, alpha = 2, beta = 2.75;

    return (1 + ac) * pow(r, beta) / (1 + ac*pow(r,alpha));
}

static double inverted_subhalo_number_density_profile(const double q)
{
    double left = 0, right = Halo[SUBHOST].R200, r = 0, delta = DBL_MAX;

    while (fabs(delta) > 1e-3) { // find root

        r = left + 0.5 * (right - left);

        delta = subhalo_number_density_profile(r) - q;

        if (delta > 0)
            right = r;
        else
            left = r;
    }

    return r;
}

static double nfw_scale_radius(const double c_nfw, const double M_t,
        const double r)
{ // Springel+ 2008 eq 7-9
    double left = 0, right = 10*Halo[SUBHOST].R_Sample[0],
           rs = 0, delta = DBL_MAX;

    while (fabs(delta) > 1e-3) { // find root

        rs = left + 0.5 * (right - left);

        delta = nfw_mass_profile(c_nfw, rs, r) - M_t;

        if (delta > 0)
            right = rs;
        else
            left = rs;
    }

    return rs;
}

static double  // wikipedia
nfw_mass_profile(const double c_nfw, const double rs, const double r)
{
    const double delta_c = Overdensity_Parameter();
    const double delta_s = delta_c /3*p3(c_nfw)
                         /(log(1+c_nfw) - c_nfw/(1+c_nfw));
    const double rho_crit = Critical_Density(Param.Redshift);
    const double rho_s = delta_s * Cosmo.Rho_crit0 / Unit.Density;

    return 4*pi*rho_s*p3(rs) * ( log((rs + r)/rs) - r/(rs+r) );
}

#define ANGLE_MAX ( 20 * DEG2RAD )
#define MAX_IMPACT_FACTOR Halo[0].R200;
#define ENERGY_ORBIT_FRACTION_SUBH 0.3

/* Velocity for bound Kepler orbits with random eccentricity,
 * plane orientation and anomaly (Binney & Tremaine pp 148) */
static void set_subhalo_bulkvel(const int i)
{
    const double G = Grav/p3(Unit.Length)*Unit.Mass*p2(Unit.Time);

    double dx = Halo[SUBHOST].D_CoM[0] - Halo[i].D_CoM[0];
    double dy = Halo[SUBHOST].D_CoM[1] - Halo[i].D_CoM[1];
    double dz = Halo[SUBHOST].D_CoM[2] - Halo[i].D_CoM[2];

    const double r = sqrt(dx*dx + dy*dy + dz*dz);

    double vel[3] = { 0 }; // velocity plane orientation

    vel[0] = erand48(Omp.Seed);
    vel[1] = erand48(Omp.Seed);
    vel[2] = erand48(Omp.Seed);

    double norm = sqrt( p2(vel[0]) + p2(vel[1]) + p2(vel[2]) );
    double impactfac = erand48(Omp.Seed) * MAX_IMPACT_FACTOR;

    vel[0] = Halo[i].D_CoM[0] - (Halo[SUBHOST].D_CoM[0] + impactfac*vel[0]/norm);
    vel[1] = Halo[i].D_CoM[1] - (Halo[SUBHOST].D_CoM[1] + impactfac*vel[1]/norm);
    vel[2] = Halo[i].D_CoM[2] - (Halo[SUBHOST].D_CoM[2] + impactfac*vel[2]/norm);

    norm = sqrt(p2(vel[0]) + p2(vel[1]) + p2(vel[2]));

    double v = ENERGY_ORBIT_FRACTION_SUBH *
        sqrt( 2*G * Halo[SUBHOST].Mtotal200 / r);

#ifdef ADD_THIRD_SUBHALO
    if (i == Sub.First) {

        v = norm = 1;

        vel[0] = -Param.SubFirstVel[0];
        vel[1] = -Param.SubFirstVel[1];
        vel[2] = -Param.SubFirstVel[2];
    }
#endif

    Halo[i].BulkVel[0] -= v * vel[0] / norm;
    Halo[i].BulkVel[1] -= v * vel[1] / norm;
    Halo[i].BulkVel[2] -= v * vel[2] / norm;

    return ;
}

#undef MIN_SUBHALO_MASS
#undef ANGLE_MAX
#undef ENERGY_ORBIT_FRACTION_SUBH

#endif // SUBSTRUCTURE
