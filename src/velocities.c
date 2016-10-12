#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include "globals.h"

#define NTABLE 512
#define NSAMPLE (NTABLE*2) // oversample integrand
#define RMIN 1 // smallest radius for which we have f(E)

static void calc_distribution_function_table(int);
static double distribution_function(const double);
static double eddington_integrant(double, void *);
static double potential_profile(const int, const float);

static double hernquist_distribution_func(const int, const double);
static double sph_kernel_wc2(const float r, const float h);

#ifdef SUBSTRUCTURE
static void set_subhalo_bulk_velocities();
#endif

void Make_velocities()
{
	const double boxhalf = Param.Boxsize/2;

    printf("\nSetting velocities "); 
	fflush(stdout);
	
	for (int i = 0; i < Param.Nhalos; i++) {

		if (i < Sub.First)
        	printf("<%d> ", i); 
		else
			printf(".");

		fflush(stdout);
		
		/* peculiar velocity */

		Setup_Profiles(i);

		calc_distribution_function_table(i); 
		
		double M = Halo[i].Mtotal;

		double dCoM[3] = {Halo[i].D_CoM[0],Halo[i].D_CoM[1],Halo[i].D_CoM[2]};

		#pragma omp parallel for schedule(dynamic) 
        for (size_t ipart = 0; ipart < Halo[i].Npart[1]; ipart++) { // DM

			double dx = Halo[i].DM[ipart].Pos[0] - dCoM[0] - boxhalf;
            double dy = Halo[i].DM[ipart].Pos[1] - dCoM[1] - boxhalf;
            double dz = Halo[i].DM[ipart].Pos[2] - dCoM[2] - boxhalf;

            double r = fmax(RMIN, sqrt(dx*dx + dy*dy + dz*dz)); 
			
			/* rejection sampling */

			double psi = potential_profile(i, r); // psi = -phi >= 0
            double vmax = sqrt(2*psi); // particle has to be bound
            double emax = psi;

			double qmax = 4*pi * p2(vmax)/M * distribution_function(emax);

            double v = 0;
            
			for (int i = 0; i < 9000; i++) { // Ascasibar+ 2005, Press+ 1992

                double lower_bound = qmax * erand48(Omp.Seed);

                v = vmax * erand48(Omp.Seed);

                double Etot = 0.5 * v*v - psi;  

				double q = 4*pi * p2(v)/M * distribution_function(-Etot);

				if (q >= lower_bound)
					break;

				v = 0;
            }

            double theta = acos(2 *  erand48(Omp.Seed) - 1);
            double phi = 2*pi * erand48(Omp.Seed);
        
            Halo[i].DM[ipart].Vel[0] = (float) (v * sin(theta) * cos(phi));
            Halo[i].DM[ipart].Vel[1] = (float) (v * sin(theta) * sin(phi));
            Halo[i].DM[ipart].Vel[2] = (float) (v * cos(theta));

        } // for ipart

#if defined(SUBSTRUCTURE) && defined(SLOW_SUBSTRUCTURE)
		if(i == SUBHOST) // at main cluster only, because we take its f(E)
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
				
				double dx = Halo[i].Gas[ipart].Pos[0] - Halo[i].D_CoM[0] 
					- boxhalf;
				double dy = Halo[i].Gas[ipart].Pos[1] - Halo[i].D_CoM[1] 
					- boxhalf;
				double dz = Halo[i].Gas[ipart].Pos[2] - Halo[i].D_CoM[2] 
					- boxhalf;

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


/* Find f(E) for arbitrary spherically symmetric density-potential pairs by 
 * numerical integration of the Eddington equation.
 * We interpolate rho(psi) with an cubic spline using the GSL library and
 * get the second derivative from the spline directly. This is a hard 1D 
 * integral to get to floating point precision ! We disable the GSL error
 * handler, because target accuracy can't always be reached. The Hernquist
 * f(E) is reproduce with a few 1e-3 accuracy ...
 * Kazantzidis+ 2004, Binney 1982, Binney & Tremaine pp. 298, Barnes 02 */  

typedef struct {
	double E;
	gsl_spline *spline;
	gsl_interp_accel *acc;
} interpolation_parameters;

static interpolation_parameters fE_params;
#pragma omp threadprivate(fE_params)

static double distribution_function(const double E)
{
	double log_E = log10(E);

	double log10_fE =  gsl_spline_eval(fE_params.spline, log_E, fE_params.acc);

	return pow(10, log10_fE);
}

static void calc_distribution_function_table(int iCluster)
{
	double rmin = Zero; // 0 wont work 
	double rmax = Infinity; // large val for accuracy

	double rstep = log10(rmax/rmin) / NSAMPLE; // good up to one Gpc

	double psi[NSAMPLE] = { 0 }, 
		   rho[NSAMPLE] = { 0 };

	#pragma omp parallel for  
	for (int i = 0; i < NSAMPLE; i++) {
		
		double r = rmin * pow(10, i*rstep);

		rho[i] = DM_Density_Profile(iCluster, r);

		psi[i] = potential_profile(iCluster, r); // psi = - phi

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

	gsl_integration_workspace *w = gsl_integration_workspace_alloc(NSAMPLE);

	interpolation_parameters int_params; 
	int_params.acc = gsl_interp_accel_alloc();
	int_params.spline = gsl_spline_alloc(gsl_interp_cspline, NSAMPLE);

	gsl_spline_init(int_params.spline, x, y, NSAMPLE);

	rstep = log10(rmax/rmin) / NTABLE; // only consider 10 Mpc radius
	
	const double sqrt8 = sqrt(8.0);

	double err = 0;

	#pragma omp for
	for (int i = 0; i < NTABLE; i++) {

		double r = rmin * pow(10, i*rstep);

		int_params.E = E[i] = potential_profile(iCluster, r);

		gsl_function F = {&eddington_integrant, &int_params};
		
		gsl_integration_qags(&F, 0, E[i], 0, 1e-3, NSAMPLE, w, &fE[i], &err); 
		
		fE[i] /= sqrt8 * pi * pi;
		
		//printf("%d %g %g %g %g %g %g \n", i, r, E[i], fE[i],
		//		hernquist_distribution_func(iCluster, E[i]), err/fE[i], 
		//		eddington_integrant(0.5*E[i], &int_params));
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
	fE_params.spline = gsl_spline_alloc(gsl_interp_akima, NTABLE);
	}
	
	int nSpline = 0;

	for (int i = 0; i < NTABLE; i++) { // interpolate fE in log space

		x[i] = log10(E[NTABLE-i-1]);
		
		y[i] = log10(fE[NTABLE-i-1]);
	}

	#pragma omp parallel
	gsl_spline_init(fE_params.spline, x, y, NTABLE);

	/*for (int i = 0; i < NTABLE; i++) {
		
		double r = rmin * pow(10, i*rstep);

		double E = potential_profile(iCluster, r);

		printf("%d %g %g ", i, r, E);

		double solution =  hernquist_distribution_func(iCluster, E);
		double fe = distribution_function(E);

		printf("%g %g %g \n", solution, fe, fabs(fe-solution)/solution);
	} */

	return ;
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


/* This is the absolute potential psi = -phi.  */
static double potential_profile(const int i, const float r)
{
	double psi = fabs(DM_Potential_Profile(i, r)); // DM generated potential
	
	if (Halo[i].Npart[0]) // Gas generated potential
		psi += fabs(Gas_Potential_Profile(i,r));

	return psi; 
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
