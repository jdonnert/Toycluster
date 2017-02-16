#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include "globals.h"

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

		double M = Halo[i].Mtotal;

		double dCoM[3] = {Halo[i].D_CoM[0],Halo[i].D_CoM[1],Halo[i].D_CoM[2]};

		#pragma omp parallel for schedule(dynamic) 
        for (size_t ipart = 0; ipart < Halo[i].Npart[1]; ipart++) { // DM
			
			double dx = Halo[i].DM[ipart].Pos[0] - dCoM[0] - boxhalf;
            double dy = Halo[i].DM[ipart].Pos[1] - dCoM[1] - boxhalf;
            double dz = Halo[i].DM[ipart].Pos[2] - dCoM[2] - boxhalf;

            double r = fmax(1, sqrt(dx*dx + dy*dy + dz*dz)); 
			
			/* rejection sampling */

			double psi = Potential_Profile(i, r); // psi = -phi >= 0
            double vmax = sqrt(2*psi); // particle has to be bound
            double emax = psi;

			double qmax = 4*pi * p2(vmax)/M * Distribution_Function(emax);

            double v = 0;
            
			for (int i = 0; i < 9000; i++) { // Ascasibar+ 2005, Press+ 1992

                double lower_bound = qmax * erand48(Omp.Seed);

                v = vmax * erand48(Omp.Seed);

                double Etot = 0.5 * v*v - psi;  

				double q = 4*pi * p2(v)/M * Distribution_Function(-Etot);

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
		if(i == SUBHOST) // at host cluster only, because we take its f(E)
			set_subhalo_bulk_velocities();
#endif
		
		/* bulk velocity */
		#pragma omp parallel for schedule(dynamic) 
		for (size_t ipart = 0; ipart < Halo[i].Npart[1]; ipart++) { // DM 
            
            Halo[i].DM[ipart].Vel[0] += Halo[i].BulkVel[0]; 
            Halo[i].DM[ipart].Vel[1] += Halo[i].BulkVel[1]; 
            Halo[i].DM[ipart].Vel[2] += Halo[i].BulkVel[2]; 
        }

		if (i < Sub.First) { // Main Halos

			#pragma omp parallel for schedule(dynamic) 
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

#ifdef SUBSTRUCTURE

/* Sample <0>'s f(E) with r relative to <0>'s centre for all subhalos
 * i.e. treat every subhalos like a particle of <0> following the 
 * Hernquist derived velocity distribution. This will give random 
 * bound orbits */

void set_subhalo_bulk_velocities()
{
	const double boxhalf = Param.Boxsize/2;
	const double M = Halo[SUBHOST].Mtotal;

	printf("\n");

	Setup_Profiles(SUBHOST);

	for (int i = Sub.First; i < Param.Nhalos; i++) { 
              
        float dx = Halo[i].D_CoM[0] - Halo[SUBHOST].D_CoM[0];
        float dy = Halo[i].D_CoM[1] - Halo[SUBHOST].D_CoM[1];
        float dz = Halo[i].D_CoM[2] - Halo[SUBHOST].D_CoM[2];
		
		double r = sqrt(dx*dx + dy*dy + dz*dz); 

		double psi = Potential_Profile(i, r);
        double vmax = sqrt(2*psi);
        double emax = psi;
		double qmax = 4*pi * p2(vmax)/M * Distribution_Function(emax);

        double v = 0;

        for (int j = 0; j < 100; j++) { 
            
			double lower_bound = qmax * erand48(Omp.Seed);

            v = vmax * erand48(Omp.Seed);

            double Etot = 0.5 * v*v - psi;  

	 		double q = 4*pi * p2(v)/M * Distribution_Function(-Etot);
			
			if (q >= lower_bound)
				break;

			v = 0;
        }
    
        double theta = acos(2 *  erand48(Omp.Seed) - 1);
        double phi = 2*pi * erand48(Omp.Seed);
        
		double U = Internal_Energy_Profile(0, r);
		double cs = sqrt(U * adiabatic_index * (adiabatic_index-1));

		v = fmin(0.5*cs,v);

		printf("Sub=%d v=%g r=%g cs=%g stripped =%d\n",
				i, v, r/Halo[SUBHOST].R200, cs, Halo[i].Is_Stripped);

        Halo[i].BulkVel[0] = (float) (v * sin(theta) * cos(phi));
        Halo[i].BulkVel[1] = (float) (v * sin(theta) * sin(phi));
        Halo[i].BulkVel[2] = (float) (v * cos(theta));
	}

	return ;
}
#endif // SUBSTRUCTURE
