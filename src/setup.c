#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#include "globals.h"

#define R200_RSAMPLE_RATIO 1.5

/* Set all relevant parameters for the collision.
 * We first find R200 from the baryon fraction in R200 and the total mass
 * in R200 given the profiles and assumptions on the scaling of the cluster
 * parameters c, beta, rho0 and rc with M200. */

void Setup()
{
    const double Xm = Param.Mass_Ratio;
    const double z = Param.Redshift;
    const double rho_crit =  Critical_Density(z);
    const double delta = 200;

    Halo[0].Mtotal200 = Param.Mtot200 / (1+Xm); // Halo Masses inside r200
    Halo[1].Mtotal200 = Param.Mtot200 - Halo[0].Mtotal200;

    if (Xm == 0)
        Param.Nhalos = 1;
    else
        Param.Nhalos = 2;

    for (int i = 0; i < Param.Nhalos; i++) {

        Halo[i].Mass200[1] = Halo[i].Mtotal200 / (1 + Halo[i].bf);
        Halo[i].Mass200[0] = Halo[i].Mtotal200 - Halo[i].Mass200[1];

        /* Kitayama & Suto 99, Boehringer+ 2012 (1), Lokas+ 2000 */
        Halo[i].R200 = pow(Halo[i].Mtotal200 * Unit.Mass
                / (delta*rho_crit*fourpithird),1./3.)/Unit.Length;

        double c_nfw = Halo[i].C_nfw = Concentration_parameter(i);

        Halo[i].Rs = Halo[i].R200 / Halo[i].C_nfw; // NFW scale radius

        /* Springel & Farrar 07 */
        Halo[i].A_hernq = Halo[i].Rs*sqrt(2*(log(1+c_nfw)-c_nfw/(1+c_nfw)));

        Halo[i].Rho0_DM = 1; // norm of DM profile
        Halo[i].Rho0_DM = Halo[i].Mass200[1] / DM_Mass_Profile_NFW(Halo[i].R200,i);
        // / DM_Mass_Profile_HQ(Halo[i].R200,i);
    }

    Param.Boxsize = floor(2*R200_TO_RMAX_RATIO * Halo[0].R200); // fit scalings

    char string[256];
    double mtot[2] = {0};

    for (int i = 0; i < Param.Nhalos; i++) { // Baryons and total mass

        Halo[i].R_Sample[0] = Halo[i].R200 * R200_RSAMPLE_RATIO;
        Halo[i].R_Sample[1] = Halo[i].R200 * R200_RSAMPLE_RATIO;
        Halo[i].Rcut = Halo[i].Rcut_R200_Ratio * Halo[i].R200;

        if (i == 0) { // 0 provides a background

            Halo[i].R_Sample[1] = Param.Boxsize/2;
            Halo[i].R_Sample[0] = Param.Boxsize/2 ;//* sqrt(3);
        }

        double rs_gas = Halo[i].R_Sample[0];
        double rs_dm = Halo[i].R_Sample[1];
        double a = Halo[i].A_hernq;
        double beta = Halo[i].Beta;
        double r200 = Halo[i].R200;
        double rcut = Halo[i].Rcut;

        Halo[i].Rcore = Gas_core_radius(i, string);

        double rc = Halo[i].Rcore;
        int cuspy = Halo[i].Have_Cuspy;

        /* DM Mass */

        Setup_DM_Mass_Profile(i);

        Halo[i].Mass[1] = DBL_MAX;
        Halo[i].Mass[1] = DM_Mass_Profile(Halo[i].R_Sample[1], i);

        Halo[i].MassCorrFac = DM_Mass_Profile(Infinity, i)
            / DM_Mass_Profile(Halo[i].R200, i);
        // correct for finite R_sample != infty
        /* Gas Mass */

        Halo[i].Mass[0] = 0;

        double rho0 = Halo[i].Rho0;

        if (Halo[i].bf > 0) {

            Halo[i].Rho0 = 1;

            Setup_Gas_Mass_Profile(i);

            Halo[i].Rho0 = Halo[i].Mass200[0] / Gas_Mass_Profile(r200, i);

            Setup_Gas_Mass_Profile(i);

            Halo[i].Mass[0] = Gas_Mass_Profile(rs_gas, i);
        }

        Halo[i].Mtotal = Halo[i].Mass[0] + Halo[i].Mass[1];

        Param.Mtotal += Halo[i].Mtotal;
        mtot[0] += Halo[i].Mass[0];
        mtot[1] +=  Halo[i].Mass[1];

        printf("Halo Setup : <%d>\n"
                "   Model             = %s\n"
                "   Sample Radius Gas = %g kpc\n"
                "   Sample Radius DM  = %g kpc\n"
                "   Mass in R_Sample  = %g 10^10 Msol\n"
                "   Mass in DM        = %g 10^10 Msol\n"
                "   Mass in Gas       = %g 10^10 Msol\n"
                "   Mass in R200      = %g 10^10 Msol\n"
                "   rho0_DM           = %g [gadget] \n"
                "   rs                = %g [kpc] \n"
                "   c_nfw             = %g \n"
                "   R200              = %g kpc\n"
                "   a_hernquist       = %g kpc\n"
                "   rho0_gas          = %g g/cm^3\n"
                "   rho0_gas          = %g [gadget]\n"
                "   rho0_gas          = %g [cm^-3]\n"
                "   beta              = %g \n"
                "   rc                = %g kpc\n"
                "   Rcut              = %g kpc\n",
                i,string,Halo[i].R_Sample[0], Halo[i].R_Sample[1],
                Halo[i].Mtotal, Halo[i].Mass[1], Halo[i].Mass[0],
                Halo[i].Mtotal200, Halo[i].Rho0_DM, Halo[i].Rs,
                Halo[i].C_nfw, Halo[i].R200*Unit.Length/kpc2cgs,
                Halo[i].A_hernq*Unit.Length/kpc2cgs,
                Density(Halo[i].Rho0), Halo[i].Rho0,
                Density(Halo[i].Rho0) / (0.6 * m_p), Halo[i].Beta,
                Halo[i].Rcore, Halo[i].Rcut);

#ifdef DOUBLE_BETA_COOL_CORES
        printf("   rho0_cc           = %g g/cm^3\n"
                "   rho0_cc           = %g [gadget]\n"
                "   rc_cc             = %g kpc\n",
                Density(Halo[i].Rho0*Param.Rho0_Fac),
                Halo[i].Rho0*Param.Rho0_Fac,Halo[i].Rcore/Param.Rc_Fac );
#endif // DOUBLE_BETA_COOL_CORES

        /* Calc & show effective Baryon Fraction in R500 */

        if (Halo[i].bf == 0) // DM only
            continue;

        if (Halo[i].Mtotal200 == 0)
            continue;

        /* R500 */
        double rho_crit =  Critical_Density(z);
        Halo[i].R500 = pow(Halo[i].Mtotal200*Unit.Mass
                / (500*rho_crit*fourpithird),1./3.)/Unit.Length;

        double r500 = Halo[i].R500 * Unit.Length;
        double Mdm = Halo[i].Mass[1] * Unit.Mass;

        rho0 = Density(Halo[i].Rho0);
        a = Halo[i].A_hernq * Unit.Length;
        rc = Halo[i].Rcore * Unit.Length;

        Halo[i].Bf_eff = 4*pi*p3(rc)*rho0 * (r500/rc - atan(r500/rc))
            / (Mdm * p2(r500) / p2(a+r500) );

        printf("   R500              = %g kpc\n"
                "   bf_200            = %g \n"
                "   bf_500            = %g \n",
                r500/Unit.Length, Halo[i].bf, Halo[i].Bf_eff);

        printf("\n");
    }

    /* Particle numbers are calculated from the global
     * masses, not from masses in r200 */

    int nDM  = 0.5 * Param.Ntotal;
    int nGas = 0.5 * Param.Ntotal;

    double mDM = mtot[1]/nDM,
		   mGas = mtot[0]/nGas;

    Param.Mpart[0] = mGas;
    Param.Mpart[1] = mDM;


    for (int i = 0; i<Param.Nhalos; i++) {

        Halo[i].Npart[0] = round(Halo[i].Mass[0]/mGas);
        Halo[i].Npart[1] = round(Halo[i].Mass[1]/mDM);

        Halo[i].Ntotal = Halo[i].Npart[0] + Halo[i].Npart[1];
    }

    if (!Halo[0].bf) { // DM only, redo masses and numbers

        Param.Mpart[0] = mGas = nGas = 0;
        Param.Mpart[1] = mDM = Param.Mtotal/Param.Ntotal;

        for (int i = 0; i < Param.Nhalos; i++) {

            Halo[i].Npart[1] = round(Halo[i].Mtotal / mDM);
            Halo[i].Npart[0] = 0;

            Halo[i].Ntotal = Halo[i].Npart[1];
        }
    }

    for (int i = 0; i < 6; i++)
        Param.Npart[i] = Halo[0].Npart[i]+Halo[1].Npart[i];

    printf("\nSystem Setup : \n"
            "   Mass Ratio      = %g\n"
            "   Boxsize         = %g kpc\n"
            "   Total Mass      = %g 10^10 Msol\n"
            "   Mass in Gas     = %g 10^10 Msol\n"
            "   Mass in DM      = %g 10^10 Msol\n"
            "   Mass Ratio      = %g\n"
            "   given bf_0      = %g\n"
            "   given bf_1      = %g\n"
            "   boxwide bf      = %g\n"
            "   # of Particles  = %lld\n"
            "   Sph Part Mass   = %g 10^10 Msol\n"
            "   DM Part Mass    = %g 10^10 Msol\n"
            "   Npart Parent    = %8lld, %8lld \n"
            "   Npart Bullet    = %8lld, %8lld \n"
            "   Npart Total     = %8lld, %8lld\n\n"
            , Param.Mass_Ratio, Param.Boxsize
            ,Param.Mtotal,mtot[0],mtot[1],Xm
            ,Halo[0].bf, Halo[1].bf, mtot[0]/mtot[1]
            ,Param.Ntotal,Param.Mpart[0],Param.Mpart[1],
            Halo[0].Npart[0], Halo[0].Npart[1],
            Halo[1].Npart[0], Halo[1].Npart[1],
            Param.Npart[0], Param.Npart[1]);

    /* allocate particles */
    size_t nBytes = Param.Ntotal * sizeof(*P);
    P = Malloc(nBytes);
    memset(P, 0, nBytes);

    nBytes = Param.Npart[0] * sizeof(*SphP);
    SphP = Malloc(nBytes);
    memset(SphP, 0, nBytes);

    /* set access pointers */
    Halo[0].Gas = &(P[0]);
    Halo[0].DM = &(P[nGas]);
    Halo[0].SphP = &(SphP[0]);

    if (Xm) {

        Halo[1].Gas = &(P[Halo[0].Npart[0]]);

        Halo[1].DM = &(P[nGas + Halo[0].Npart[1]]);

        Halo[1].SphP = &(SphP[Halo[0].Npart[0]]);
    }

    /* grav softening  from larger cluster */
    Param.GravSofteningLength = pow(
            p3(Halo[0].R_Sample[1])/Param.Ntotal, 1/3.) / 7;

    printf("\nGrav. Softening ~ %g kpc\n\n",
            Param.GravSofteningLength*Unit.Length/kpc2cgs);

    /* kinematics */

	double d_clusters;

    if (Xm) { // two clusters only

        d_clusters = 0.9 * (Halo[0].R200 + Halo[1].R200);

        Halo[0].D_CoM[0] = -1 * Halo[1].Mtotal200
            *d_clusters/Param.Mtot200;
        Halo[1].D_CoM[0] = d_clusters + Halo[0].D_CoM[0];

        /* if (Xm > 1) {

           Halo[1].D_CoM[0] = -1 * Halo[1].Mtotal200
         *d_clusters/Param.Mtot200;

         Halo[0].D_CoM[0] = d_clusters + Halo[0].D_CoM[0];
         }*/

        Halo[0].D_CoM[1] = -1 * Halo[1].Mtotal200
            * Param.Impact_Param/Param.Mtot200;

        Halo[1].D_CoM[1] = Param.Impact_Param + Halo[0].D_CoM[1];

// #ifndef GIVEPARAMS
        Param.VelMerger[0] = sqrt(2*G*Halo[1].Mtotal200
                /(d_clusters*(1+1/Xm)));
        Param.VelMerger[1] = -Param.Mtot200/Halo[1].Mtotal200
            * Param.VelMerger[0];

        Param.VelMerger[0] *= Param.Zero_Energy_Orbit_Fraction;
        Param.VelMerger[1] *= Param.Zero_Energy_Orbit_Fraction;
// #endif

        Halo[0].BulkVel[0] = Halo[0].BulkVel[1] = Halo[0].BulkVel[2] = 0;
        Halo[1].BulkVel[0] = Halo[1].BulkVel[1] = Halo[1].BulkVel[2] = 0;

#if !defined(PARABOLA) && !defined(COMET) // apply merger kinematics directly
        Halo[0].BulkVel[0] = Param.VelMerger[0];
        Halo[1].BulkVel[0] = Param.VelMerger[1];
#endif

        printf("Kinematics of Collision : \n"
                "   Zero-E fraction     = %g \n"
                "   Initial Distance    = %g kpc\n"
                "   CoM Distance of <0> = %g kpc\n"
                "   CoM Distance of <1> = %g kpc\n"
                "   CoM Velocity of <0> = %g km/s\n"
                "   CoM Velocity of <1> = %g km/s\n\n"
                "   Impact Parameter    = %g kpc\n"
                "   CoM Impact of <0>   = %g kpc\n"
                "   CoM Impact of <1>   = %g kpc\n\n"
                ,Param.Zero_Energy_Orbit_Fraction
                ,d_clusters, Halo[0].D_CoM[0], Halo[1].D_CoM[0]
                , Param.VelMerger[0], Param.VelMerger[1]
                ,Param.Impact_Param, Halo[0].D_CoM[1]
                ,Halo[1].D_CoM[1]);

    } else { // just one cluster

        d_clusters = 0;

        Halo[0].D_CoM[0] = 0;
        Halo[0].D_CoM[1] = 0;
        Halo[0].D_CoM[2] = 0;
        Param.VelMerger[0] = Param.VelMerger[1] = 0;
    }

#ifdef SUBSTRUCTURE
    printf("\nSubhalos hosted by cluster <%d> \n\n",SUBHOST );

    Setup_Substructure();
#endif

    return;
}

static double sph_kernel_wc2(const float r, const float h)
{
    double u = r/h;
    double t = fmax(1 - u, 0);

    return 21/2/pi/p3(h)*t*t*t*t * (1+4*u);
}

/* set orbit */
void Apply_kinematics()
{
    const float vx_host = Param.VelMerger[0];
    const float vx_infa = Param.VelMerger[1];

#ifdef PARABOLA // move origin to R200 touch point
    float dx =  - Halo[1].D_CoM[0] + Param.Boxsize/2 +Halo[1].R200;
    float dy =  - Halo[1].D_CoM[1] + Param.Boxsize/2;
    float dz =  - Halo[1].D_CoM[2] + Param.Boxsize/2;

#pragma omp parallel for
    for (size_t ipart = 0; ipart < Param.Ntotal; ipart++) {

        float x = P[ipart].Pos[0] +dx;
        float y = P[ipart].Pos[1] +dy;
        float z = P[ipart].Pos[2] +dz;

        if ( y*y + z*z < x*x && (x > 0) )  // infalling cluster
            P[ipart].Vel[0] += vx_infa;
        else                              // host cluster
            P[ipart].Vel[0] += vx_host;
    }
#endif // PARABOLA

#if defined(COMET)

    const double boxhalf = Param.Boxsize/2;

    const float x0 = Halo[1].D_CoM[0] + boxhalf;
    const float y0 = Halo[1].D_CoM[1] + boxhalf;
    const float z0 = Halo[1].D_CoM[2] + boxhalf;

    const float rVir2 = p2(Halo[1].R200);

    const double h = Halo[1].R_Sample[0];
    const double norm = sph_kernel_wc2(0, h);

#pragma omp parallel for
    for (size_t ipart = 0; ipart < Param.Ntotal; ipart++) {

        float dx = P[ipart].Pos[0] - x0;
        float dy = P[ipart].Pos[1] - y0;
        float dz = P[ipart].Pos[2] - z0;

        float r2_cyl = p2(dy) + p2(dz);
        float r2 = p2(dx) + p2(dy) + p2(dz);

        if( ((dx > 0) && (r2_cyl < rVir2)) || (r2 < rVir2) ) {

            if (dx < 0)    { // slow down front part

                double r = sqrt(dx*dx + dy*dy + dz*dz);

                double wk = 1; // sph_kernel_wc2(r,  h)/norm;

                P[ipart].Vel[0] += vx_infa * wk;

            } else { // std velocity

                P[ipart].Vel[0] += vx_infa;
            }

        } else {
            P[ipart].Vel[0] += vx_host;
        }
    }
#endif // COMET

    return;
}

/* center the simulation in the periodic box */
void Shift_Origin()
{
    const float boxsize = Param.Boxsize;
    const float boxHalf = boxsize / 2;

    printf("Shift Origin "); fflush(stdout);

    for (int i = 0; i < Param.Nhalos; i++) { // shift clusters away from 0

        printf("."); fflush(stdout);

        float dx = Halo[i].D_CoM[0];
        float dy = Halo[i].D_CoM[1];
        float dz = Halo[i].D_CoM[2];

        float vx = Halo[i].BulkVel[0];
        float vy = Halo[i].BulkVel[1];
        float vz = Halo[i].BulkVel[2];

#pragma omp parallel for // DM
        for (size_t ipart=0; ipart<Halo[i].Npart[1]; ipart++) {

            Halo[i].DM[ipart].Pos[0] += dx;
            Halo[i].DM[ipart].Pos[1] += dy;
            Halo[i].DM[ipart].Pos[2] += dz;

            Halo[i].DM[ipart].Vel[0] += vx;
            Halo[i].DM[ipart].Vel[1] += vy;
            Halo[i].DM[ipart].Vel[2] += vz;
        }

#pragma omp parallel for  // Gas
        for (size_t ipart=0; ipart<Halo[i].Npart[0]; ipart++) {

            Halo[i].Gas[ipart].Pos[0] += dx;
            Halo[i].Gas[ipart].Pos[1] += dy;
            Halo[i].Gas[ipart].Pos[2] += dz;

            Halo[i].Gas[ipart].Vel[0] += vx;
            Halo[i].Gas[ipart].Vel[1] += vy;
            Halo[i].Gas[ipart].Vel[2] += vz;
        }
    }

    printf("done \n");

#pragma omp parallel for // shift 0 to bottom left corner, wrap around box
    for (size_t ipart=0; ipart<Param.Ntotal; ipart++) {

        P[ipart].Pos[0] += boxHalf;
        P[ipart].Pos[1] += boxHalf;
        P[ipart].Pos[2] += boxHalf;

        while (P[ipart].Pos[0] > boxsize)
            P[ipart].Pos[0] -= boxsize;

        while (P[ipart].Pos[0] < 0)
            P[ipart].Pos[0] += boxsize;

        while (P[ipart].Pos[1] > boxsize)
            P[ipart].Pos[1] -= boxsize;

        while (P[ipart].Pos[1] < 0)
            P[ipart].Pos[1] += boxsize;

        while (P[ipart].Pos[2] > boxsize)
            P[ipart].Pos[2] -= boxsize;

        while (P[ipart].Pos[2] < 0)
            P[ipart].Pos[2] += boxsize;
    }

    return ;
}

/* NFW concentration parameter for main and subhalos */
double Concentration_parameter(const int i)
{
#ifdef GIVEPARAMS
    if (i < Sub.First)
        return Halo[i].C_nfw;
#else
    const double mass = Halo[i].Mtotal200 * Unit.Mass/Msol2cgs;

#ifdef NFWC_DUFFY08   // simulation fit in WMAP5 cosmology

    const double A = 5.74;
    const double B = -0.097;
    const double C = -0.47;
    const double Mpivot = 2e12 / Cosmo.h_100;

    double c_NFW = A * pow( mass/Mpivot, B) * pow( 1+Param.Redshift, C);

#endif // NFWC_DUFFY08

#ifdef NFWC_BUOTE07   // from observations in conc.cosm. (Buote+ 07)

    double c_NFW = 9*pow( mass / (1e14*Msol2cgs), -0.172);

#endif // NFWC_BUOTE07

#ifdef SUBSTRUCTURE

    if (i >= Sub.First) { // Giocoli 2010 eq 16

        double mass_sub = Halo[i].Mass[1] * Unit.Mass/Msol2cgs;

        const double aR = 0.237, c1 = 232.15, c2 = -181.74,
              a1 = 0.0146, a2 = 0.008;  // Pieri 2009, fit to Aquarius

        double dx = Halo[SUBHOST].D_CoM[0] - Halo[i].D_CoM[0];
        double dy = Halo[SUBHOST].D_CoM[1] - Halo[i].D_CoM[1];
        double dz = Halo[SUBHOST].D_CoM[2] - Halo[i].D_CoM[2];

        double d_vir = sqrt(dx*dx + dy*dy + dz*dz) / Halo[0].R200;

        c_NFW = pow(d_vir, -aR) * // Pieri 2009
            ( c1*pow(mass_sub, -a1) + c2*pow(mass_sub,-a2) );
        c_NFW /= 1+Param.Redshift;
    }

#endif  // SUBSTRUCTURE

    return c_NFW ; // * sqrt(Omega_M(Param.Redshift))
#endif  // GIVEPARAMS
}

/* r_core - beta model core radius model */
double Gas_core_radius(const int i, char *string)
{
    double rc = 0;

#ifdef GIVEPARAMS
    sprintf(string, "rc from parameter file ");

    if (i < Sub.First)
        return Halo[i].Rcore;
#endif // GIVEPARAMS


    if (Param.Cuspy & (1 << i) ) {

        sprintf(string, "Cool Core ");

        rc = Halo[i].Rs / 9;

        Halo[i].Have_Cuspy = 1;

#ifdef DOUBLE_BETA_COOL_CORES

        sprintf(string, "Double Beta Cool Core (%g, %g)",
                Param.Rho0_Fac, Param.Rc_Fac);

        rc = Halo[i].Rs / 3;
#endif
    } else {

        sprintf(string, "Disturbed ");

        rc = Halo[i].Rs / 3;

        Halo[i].Have_Cuspy = 0;
    }

    return rc;
}
