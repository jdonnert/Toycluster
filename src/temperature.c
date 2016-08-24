#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#include "globals.h"

static void setup_internal_energy_profile(const int i);

void Make_temperatures()
{
    const double boxhalf = 0.5 * Param.Boxsize;

    printf("Setting temperatures "); 
    
    for (int i = 0; i < Param.Nhalos; i++) {

		if (i < Sub.First)
        	printf("<%d> ",i);
		else
			printf(".");

		fflush(stdout);

		setup_internal_energy_profile(i);

		#pragma omp parallel for
        for (int ipart = 0; ipart < Halo[i].Npart[0]; ipart++) {

            float dx = Halo[i].Gas[ipart].Pos[0] - Halo[i].D_CoM[0] - boxhalf;
            float dy = Halo[i].Gas[ipart].Pos[1] - Halo[i].D_CoM[1] - boxhalf;
            float dz = Halo[i].Gas[ipart].Pos[2] - Halo[i].D_CoM[2] - boxhalf;

            double r = sqrt(dx*dx + dy*dy + dz*dz);

            double u = Internal_Energy_Profile(i, r);
			double u_ana = Internal_Energy_Profile_Analytic(i, r);
			
			Halo[i].SphP[ipart].U = u;
		}
    }
    
	printf("done\n\n"); fflush(stdout);

    return;
}

/* 
 * Standard analytical temperature profile from Donnert et al. 2016.
 * To avoid negative temperatures we define rmax*sqrt3 as outer radius
 */

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
    const double G = Grav/p3(Unit.Length)*Unit.Mass*p2(Unit.Time);
		
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

/* 
 * Numerical solution for all kinds of gas densities. We spline interpolate 
 * a table of solutions for speed. We have to solve the hydrostatic equation
 * (eq. 9 in Donnert 2014).
 */

#define TABLESIZE 1024

static gsl_spline *U_Spline = NULL;
static gsl_interp_accel *U_Acc = NULL;
#pragma omp threadprivate(U_Spline, U_Acc)

double Internal_Energy_Profile(const int i, const double r)
{
	return gsl_spline_eval(U_Spline, r, U_Acc);;
}

static double u_integrant(double r, void *param) // Donnert 2014, eq. 9
{
	int i = *((int *) param);

	double rho0 = Halo[i].Rho0;
	double rc = Halo[i].Rcore;
	double beta = Halo[i].Beta;
	double rcut = Halo[i].Rcut;
	int is_cuspy = Halo[i].Have_Cuspy;
	double a = Halo[i].A_hernq;
	double Mdm = Halo[i].Mass[1];

#ifdef NO_RCUT_IN_T
	rcut = 1e5;
#endif

	double rho_gas = Gas_density_profile(r, rho0, beta, rc, rcut, is_cuspy);
	double Mr_Gas = Mass_profile(r, i);
	double Mr_DM = Mdm * r*r/p2(r+a);

	return rho_gas /(r*r) * (Mr_Gas + Mr_DM);
}

static void setup_internal_energy_profile(const int i)
{
    double G = Grav/p3(Unit.Length)*Unit.Mass*p2(Unit.Time);

	gsl_function gsl_F = { 0 };

	double u_table[TABLESIZE] = { 0 }, r_table[TABLESIZE] = { 0 };

	double rmin = 0.1;
	double rmax = Param.Boxsize * sqrt(3);
	double dr = ( log10(rmax/rmin) ) / (TABLESIZE-1);

	Setup_Mass_Profile(i);

	#pragma omp parallel 
	{
	
	gsl_integration_workspace *gsl_workspace = NULL;
	gsl_workspace = gsl_integration_workspace_alloc(2*TABLESIZE);
	
	#pragma omp for
	for (int j = 1; j < TABLESIZE;  j++) {
	
		double error = 0;

		double r = rmin * pow(10, dr * j);
		
		r_table[j] = r;

		gsl_F.function = &u_integrant;
		gsl_F.params = (void *) &i;

		gsl_integration_qag(&gsl_F, r, rmax, 0, 1e-5, 2*TABLESIZE, 
				GSL_INTEG_GAUSS41, gsl_workspace, &u_table[j], &error);

		double rho0 = Halo[i].Rho0;
		double rc = Halo[i].Rcore;
		double beta = Halo[i].Beta;
		double rcut = Halo[i].Rcut;
		int is_cuspy = Halo[i].Have_Cuspy;
#ifdef NO_RCUT_IN_T
		rcut = 1e6;
#endif
		double rho_gas = Gas_density_profile(r, rho0, beta, rc, rcut, is_cuspy);

		u_table[j] *= G/((adiabatic_index-1)*rho_gas); // Donnert 2014, eq. 9
		
		//printf("%d %g %g %g %g \n", j,r, rho_gas, u_table[j], 
		//						u_integrant(r, (void *)&i ));
	}

	u_table[0] = u_table[1];
	//u_table[TABLESIZE-1] = 0;

	gsl_integration_workspace_free(gsl_workspace);
	
	U_Acc = gsl_interp_accel_alloc();

	U_Spline = gsl_spline_alloc(gsl_interp_cspline, TABLESIZE);
	gsl_spline_init(U_Spline, r_table, u_table, TABLESIZE);
	
	} // omp parallel

	return ;
}
