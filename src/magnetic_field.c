#include "globals.h"
#include "tree.h"

#define BMAX 18e-6
#define KLOWCUT (2*pi / Param.Boxsize * 32)
#define KHIGHCUT (2*pi / 100)
#define SPECTRAL_INDEX (-11.0/3.0)

static void set_magnetic_vector_potential();
static void normalise_magnetic_field();

void Make_magnetic_field()
{
    printf("Magnetic field: \n"
            "   B0              = %g G\n"
            "   eta             = %g \n\n"
            ,Param.Bfld_Norm, Param.Bfld_Eta);

    set_magnetic_vector_potential();

	Bfld_from_rotA_SPH(); 

 	normalise_magnetic_field();

    return;
}

/* 
 * Bonafede 2010 scaling with central density 
 * this way we also can avoid smoothing of field 
 */

static void set_magnetic_vector_potential()
{
	const float boxhalf = 0.5 * Param.Boxsize;

		#pragma omp parallel for 
    	for (size_t ipart = 0; ipart < Param.Npart[0]; ipart++) {
		
			double A_max = 0;

			for (int i = 0; i < Param.Nhalos; i++) {

			if (Halo[i].Mass[0] == 0) // DM only halos
				continue;

			float dx = P[ipart].Pos[0] - Halo[i].D_CoM[0] - boxhalf,
				  dy = P[ipart].Pos[1] - Halo[i].D_CoM[1] - boxhalf,
				  dz = P[ipart].Pos[2] - Halo[i].D_CoM[2] - boxhalf;

			double r2 = dx*dx + dy*dy + dz*dz;

			double rho_i = Gas_density_profile(sqrt(r2), Halo[i].Rho0,
						Halo[i].Beta, Halo[i].Rcore, Halo[i].Rcut, 
						Halo[i].Have_Cuspy);

			double A = pow(rho_i/Halo[i].Rho0, Param.Bfld_Eta);

			if (A > A_max)
				A_max = A;
			}
		
			SphP[ipart].Apot[0] = (float) A_max;
			SphP[ipart].Apot[1] = (float) A_max;
	    	SphP[ipart].Apot[2] = (float) A_max;
		}

    return ;  
}

static void normalise_magnetic_field()
{
	const float boxhalf = 0.5 * Param.Boxsize;

 	double max_B2 = 0;

	#pragma omp parallel for 
	for (size_t ipart = 0; ipart < Param.Npart[0]; ipart++) {
		
		double bfld2 = p2(SphP[ipart].Bfld[0]) 
				+ p2(SphP[ipart].Bfld[1]) 
				+ p2(SphP[ipart].Bfld[2]);

		max_B2 = fmax(max_B2, bfld2);
	}

	double max_bfld = sqrt(max_B2);

   	double norm = Param.Bfld_Norm/max_bfld/ sqrt(3);

	printf("Bfld Norm = %g \n", norm );

	int cnt = 0;

	#pragma omp parallel for reduction(+:cnt) 
   	for (size_t ipart = 0; ipart < Param.Npart[0]; ipart++) {
		
		SphP[ipart].Bfld[0] *= norm;
		SphP[ipart].Bfld[1] *= norm;
		SphP[ipart].Bfld[2] *= norm;

		double B2 = p2(SphP[ipart].Bfld[0]) + p2(SphP[ipart].Bfld[1]) + 
					p2(SphP[ipart].Bfld[2]);
	
		float x = P[ipart].Pos[0] - boxhalf,
			  y = P[ipart].Pos[1] - boxhalf,
			  z = P[ipart].Pos[2] - boxhalf;
			
		int i = Halo_containing(ipart,x,y,z);
		
		double bmax = BMAX;

		if (i > 1) // subhaloes
			bmax = 2e-6;

		if ((B2 > p2(bmax))) {

			double B = sqrt(B2);

			SphP[ipart].Bfld[0] *= bmax/B;
			SphP[ipart].Bfld[1] *= bmax/B;
			SphP[ipart].Bfld[2] *= bmax/B;

			cnt++;
		}
	}

	printf("Bfld of %d particles limited to %g G\n", cnt, BMAX);

	return ;
}

