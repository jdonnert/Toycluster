#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include "globals.h"

#define NTABLE 512
#define NSAMPLE (NTABLE*4) // oversample integrand
#define RMIN 0.1 // smallest radius for which we have f(E)

static void calc_distribution_function_table(int);
static double distribution_function(const double);
static double eddington_integrant(double, void *);
static double potential_profile(const int, const float);
static double gas_potential_profile(const int i, const float r);
static double dm_density_profile(const int, const float);
static double dm_potential_profile(const int, const float);
static double hernquist_distribution_func(const int, const double);
static double sph_kernel_wc2(const float r, const float h);

#ifdef SUBSTRUCTURE
static void set_subhalo_bulk_velocities();
#endif

static double G = 0; 

typedef struct {
	double E;
	gsl_spline *spline;
	gsl_interp_accel *acc;
} interpolation_parameters;

static interpolation_parameters fE_params;
#pragma omp threadprivate(fE_params)

void Make_velocities()
{
    printf("\nSetting velocities ");fflush(stdout);

	const double boxhalf = Param.Boxsize/2;

	G = Grav / p3(Unit.Length) * Unit.Mass * p2(Unit.Time); // internal units
	
	for (int i = 0; i < Param.Nhalos; i++) {

		if (i < Sub.First)
        	printf("<%d> ", i); 
		else 
			printf(".");

		fflush(stdout);
		
		/* peculiar velocity */
		calc_distribution_function_table(i); 

		double M = Halo[i].Mtotal;

		double dCoM[3] = {Halo[i].D_CoM[0],Halo[i].D_CoM[1],Halo[i].D_CoM[2]};

		//#pragma omp parallel for schedule(dynamic) 
        for (size_t ipart = 0; ipart < Halo[i].Npart[1]; ipart++) { // DM

			double dx = Halo[i].DM[ipart].Pos[0] - dCoM[0] - boxhalf;
            double dy = Halo[i].DM[ipart].Pos[1] - dCoM[1] - boxhalf;
            double dz = Halo[i].DM[ipart].Pos[2] - dCoM[2] - boxhalf;

            double r = fmax(RMIN, sqrt(dx*dx + dy*dy + dz*dz)); 
			
			double pot = -potential_profile(i, r); // rejection sampling
            double vmax = sqrt(-2*pot);
            double emax = pot;
			double qmax = 4*pi * p2(vmax)/M * distribution_function(-emax);

            double v = 0;
            
			for (;;) { // Ascasibar+ 2005, Press+ 1992

                double lower_bound = qmax * erand48(Omp.Seed);

                v = vmax * erand48(Omp.Seed);

                double Etot = 0.5 * v*v + pot;  
				
				double q = 4*pi * p2(v)/M * distribution_function(-Etot);

				if (q >= lower_bound)
					break;
            }
    
            double theta = acos(2 *  erand48(Omp.Seed) - 1);
            double phi = 2*pi * erand48(Omp.Seed);
        
            Halo[i].DM[ipart].Vel[0] = (float) (v * sin(theta) * cos(phi));
            Halo[i].DM[ipart].Vel[1] = (float) (v * sin(theta) * sin(phi));
            Halo[i].DM[ipart].Vel[2] = (float) (v * cos(theta));

        } // for ipart

#if defined(SUBSTRUCTURE) && defined(SLOW_SUBSTRUCTURE)
		if(i == 0) // at main cluster only, because we take it f(E)
			set_subhalo_bulk_velocities();
#endif
		
		/* bulk velocity */
		for (size_t ipart = 0; ipart < Halo[i].Npart[1]; ipart++) { // DM 
            
            Halo[i].DM[ipart].Vel[0] += Halo[i].BulkVel[0]; 
            Halo[i].DM[ipart].Vel[1] += Halo[i].BulkVel[1]; 
            Halo[i].DM[ipart].Vel[2] += Halo[i].BulkVel[2]; 
        }

		if (i < Sub.First) { // Main Halos

	        for (size_t ipart = 0; ipart < Halo[i].Npart[0]; ipart++) { // gas 
    	        
        	    Halo[i].Gas[ipart].Vel[0] += Halo[i].BulkVel[0]; 
	            Halo[i].Gas[ipart].Vel[1] += Halo[i].BulkVel[1]; 
            	Halo[i].Gas[ipart].Vel[2] += Halo[i].BulkVel[2]; 
    	    }

		} else { // subhalos

			double h = Halo[i].R_Sample[0] * 1.1;
			double norm = sph_kernel_wc2(0,  h);
	
			#pragma omp parallel for schedule(dynamic)  
	        for (int ipart = 0; ipart < Halo[i].Npart[0]; ipart++) { // gas 
				
				double dx = Halo[i].Gas[ipart].Pos[0] - Halo[i].D_CoM[0] - boxhalf;
				double dy = Halo[i].Gas[ipart].Pos[1] - Halo[i].D_CoM[1] - boxhalf;
				double dz = Halo[i].Gas[ipart].Pos[2] - Halo[i].D_CoM[2] - boxhalf;

				double r = sqrt(dx*dx + dy*dy + dz*dz);

				double wk = sph_kernel_wc2(r,  h)/norm;
			
				Halo[i].Gas[ipart].Vel[0] += Halo[i].BulkVel[0] * wk;
				Halo[i].Gas[ipart].Vel[1] += Halo[i].BulkVel[1] * wk;
				Halo[i].Gas[ipart].Vel[2] += Halo[i].BulkVel[2] * wk;
			}
		}

        
	} // for i Halos

    printf("done \n\n");

    return;
}

static double sph_kernel_wc2(const float r, const float h)
{
	double u = r/h;
	double t = fmax(1 - u, 0);
	
	return 21/2/pi/p3(h)*t*t*t*t * (1+4*u);
}

/* return spline interpolation from tables fE and E */
static double distribution_function(const double E)
{
	return gsl_spline_eval(fE_params.spline, E, fE_params.acc);
}

/* 
 * Find f(E) for arbitrary spherically symmetric density-potential pairs by 
 * numerical integration of the Eddington equation.
 * We interpolate rho(psi) with an cubic spline using the GSL library and
 * get the second derivative from the spline directly. This is a hard 1D 
 * integral to get to floating point precision ! We disable the GSL error
 * handler, because target accuracy can't always be reached. The Hernquist
 * f(E) is reproduce with a few 1e-3 accuracy ...
 * Kazantzidis+ 2004, Binney 1982, Binney & Tremaine pp. 298, Barnes 02 
 */  

static void calc_distribution_function_table(int iCluster)
{
	const double rmin = RMIN; // zero wont work 
	const double rmax = 1e10 * Param.Boxsize; // large val for accuracy

	double rstep = log10(1e5 * rmax) / NSAMPLE; // good up to one Gpc

	double psi[NSAMPLE] = { 0 }, rho[NSAMPLE] = { 0 };
	
	#pragma omp parallel for  
	for (int i = 0; i < NSAMPLE; i++) {
		
		double r = rmin * pow(10, (i)*rstep) - 0.999;

		rho[i] = dm_density_profile(iCluster, r);

		psi[i] = potential_profile(iCluster, r); 
	}
	
	double E[NTABLE] = { 0 }, 
		   fE[NTABLE] = { 0 };

	double x[NSAMPLE] = { 0 }, // holds inverted vals
		   y[NSAMPLE] = { 0 }; 

	for (int i = 0; i < NSAMPLE; i++) { // for spline interpolation
	
		x[i] = psi[NSAMPLE-i-1];
		y[i] = rho[NSAMPLE-i-1];
	}

	gsl_error_handler_t * old_handler = gsl_set_error_handler_off();

	#pragma omp parallel // make f(E) table
	{

	gsl_integration_workspace *w = gsl_integration_workspace_alloc(NSAMPLE);

	interpolation_parameters int_params; 
	int_params.acc = gsl_interp_accel_alloc();
	int_params.spline = gsl_spline_alloc(gsl_interp_cspline, NSAMPLE);

	gsl_spline_init(int_params.spline, x, y, NSAMPLE);

	rstep = log10(rmax/rmin) / NTABLE; // only consider 10 Mpc radius
	
	const double sqrt8 = sqrt(8);

	#pragma omp for
	for (int i = 0; i < NTABLE; i++) {

		double r = rmin * pow(10, i*rstep);

		int_params.E = E[i] = potential_profile(iCluster, r);

		gsl_function F = {&eddington_integrant, &int_params};
		
		double error = 0;

		gsl_integration_qags(&F, 0, E[i], 0, 1e-4, NSAMPLE, w, &fE[i], &error); 

		fE[i] /= sqrt8 * pi * pi;
	
	//	printf("%d %g %g %g %g %g \n", i, r, E[i], fE[i], 
	//			hernquist_distribution_func(iCluster, E[i]), error/fE[i]);
	}

	gsl_integration_workspace_free(w);

	gsl_interp_accel_free(int_params.acc);
	gsl_spline_free(int_params.spline);

	} // omp parallel 

	gsl_set_error_handler(old_handler);

	#pragma omp parallel
	{
	fE_params.acc = gsl_interp_accel_alloc();
	fE_params.spline = gsl_spline_alloc(gsl_interp_cspline, NTABLE);
	}

	for (int i = 0; i < NTABLE; i++) {
		
		x[i] = E[NTABLE-i-1];
		
		y[i] = fE[NTABLE-i-1];
	}

	#pragma omp parallel
	gsl_spline_init(fE_params.spline, x, y, NTABLE);
	
	/*for (int i = 0; i < 2*NTABLE; i++) {
		
		double r = rmin * pow(10, i*rstep/2.01);

		double E = potential_profile(iCluster, r);

		printf("%d %g %g ", i, r, E);

		double solution =  hernquist_distribution_func(iCluster, E);
		double fe = distribution_function(E);

		printf("%g %g %g \n", solution, fe, fabs(fe-solution)/solution);
	}*/

	return ;
}

/* 
 * Binney & Tremaine sect. 4.3.1, we take the second derivate from the
 * cubic spline directly 
 */

static double eddington_integrant(double psi, void *params)
{
	const interpolation_parameters *p = params;

	if (psi == p->E)
		return 0;

	const double dRhodPsi2 = gsl_spline_eval_deriv2(p->spline, psi, p->acc);
	
	return dRhodPsi2 / sqrt(p->E - psi);
}

/* 
 * Hernquist 1989, eq. 2 
 */

static double dm_density_profile(const int i, const float r)
{
    const double a = Halo[i].A_hernq;
    const double m = Halo[i].Mass[1]; 

	return m/(2.0*pi) * a/r /p3(r+a) ;
}

static double potential_profile(const int i, const float r)
{
	const double rho0 = Halo[i].Rho0;
	const double rc = Halo[i].Rcore;
	const double rcut = Halo[i].Rcut;
	
	double psi = dm_potential_profile(i,r); // DM generated potential

	if (Halo[i].Npart[0]) // Gas generated potential
		psi += gas_potential_profile(i,r);  // exponential cut off near rcut

	return psi; 
}

/* 
 * For testing. Hernquist 1989, eq. 17-19  
 */

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

/* 
 * This is Psi = -Phi, i.e. Psi(r<inf) >= 0 
 */

double dm_potential_profile(const int i, const float r)
{
    const double a = Halo[i].A_hernq;
    const double mDM = Halo[i].Mass[1]; 
	
	double psi = G * mDM / (r+a);

	return psi;
}

double gas_potential_profile(const int i, const float r)
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
//			* p2(p2(1 - r/rcut)) // no cutoff

	return 4*pi*G * psi;
}

#ifdef SUBSTRUCTURE

/* 
 * Sample <0>'s f(E) with r relative to <0>'s centre for all subhalos
 * i.e. treat every subhalos like a particle of <0> following the 
 * Hernquist derived velocity distribution. This will give random 
 * bound orbits 
 */

void set_subhalo_bulk_velocities()
{
	const float d0[3] = {Halo[SUBHOST].D_CoM[0], Halo[SUBHOST].D_CoM[1], 
						Halo[SUBHOST].D_CoM[2]};

	printf("\n");

	for (int i = Sub.First; i < Param.Nhalos; i++) { 
		
		double M = Halo[i].Mtotal;
               
        float dx = Halo[i].D_CoM[0] - d0[0];
        float dy = Halo[i].D_CoM[1] - d0[1];
        float dz = Halo[i].D_CoM[2] - d0[2];
		
		double r = sqrt(dx*dx + dy*dy + dz*dz); 

		double pot = -potential_profile(i, r);
        double vmax = sqrt(-2*pot);
        double emax = -pot;
		double qmax = 4*pi * p2(vmax)/M * distribution_function(emax);

        double v = 0;

        for (;;) { 

            double lower_bound = qmax * erand48(Omp.Seed);

            v = vmax * erand48(Omp.Seed);

            double Etot = 0.5 * v*v + pot;  

	 		double q = 4*pi/M * p2(v) * distribution_function(-Etot);

			if (q >= lower_bound)
				break;
        }
    
        double theta = acos(2 *  erand48(Omp.Seed) - 1);
        double phi = 2*pi * erand48(Omp.Seed);
        
		double U = Internal_Energy_Profile(0, r);
		double cs = sqrt(U * adiabatic_index * (adiabatic_index-1));

		v *= Param.Zero_Energy_Orbit_Fraction;

        Halo[i].BulkVel[0] = (float) (v * sin(theta) * cos(phi));
        Halo[i].BulkVel[1] = (float) (v * sin(theta) * sin(phi));
        Halo[i].BulkVel[2] = (float) (v * cos(theta));

		printf("Sub=%d v=%g r=%g cs=%g Gas Stripped ?=%d\n",
				i, v, r/Halo[SUBHOST].R200, cs, Halo[i].Is_Stripped);
	}

	return ;
}
#endif // SUBSTRUCTURE

#undef NTABLE
#undef NSAMPLE
