#include "globals.h"
#include "tree.h"
//#include "magnetic_field_turb.h"

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
#ifndef TURB_B_FIELD
	Bfld_from_rotA_SPH(); 
#else
	Bfld_from_turb_spectrum();
#endif

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

				double rho_i = Gas_Density_Profile(sqrt(r2), i);

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
/* Normalise BFLD inside 0.8 rc of halo 0 */

static void normalise_magnetic_field() // doesnt work correctly
{
	const float boxhalf = 0.5 * Param.Boxsize;

 	double mean_B = 0; // in 0.8rc of Halo[0]
	uint64_t cnt = 0;

	#pragma omp parallel for reduction(+:mean_B,cnt)
	for (size_t ipart = 0; ipart < Halo[0].Npart[0]; ipart++) {
		
		double dx = P[ipart].Pos[0] - Halo[0].D_CoM[0] - boxhalf,
			   dy = P[ipart].Pos[1] - Halo[0].D_CoM[1] - boxhalf,
			   dz = P[ipart].Pos[2] - Halo[0].D_CoM[2] - boxhalf;

		double r = sqrt(dx*dx + dy*dy + dz*dz);

		if (r > 0.8*Halo[0].Rcore)
			continue;

		double bfld = sqrt(p2(SphP[ipart].Bfld[0]) + p2(SphP[ipart].Bfld[1]) 
						 + p2(SphP[ipart].Bfld[2]));

		mean_B += bfld;
		cnt++;
	}

	mean_B /= cnt;

   	double norm = Param.Bfld_Norm/mean_B*sqrt(2);

	printf("Bfld Norm = %g \n", norm );

	cnt = 0;

	#pragma omp parallel for reduction(+:cnt) 
   	for (size_t ipart = 0; ipart < Param.Npart[0]; ipart++) {
		
		SphP[ipart].Bfld[0] *= norm;
		SphP[ipart].Bfld[1] *= norm;
		SphP[ipart].Bfld[2] *= norm;

	
		double 	x = P[ipart].Pos[0] - boxhalf,
			  	y = P[ipart].Pos[1] - boxhalf,
			 	z = P[ipart].Pos[2] - boxhalf;
			
		int i = Halo_containing(ipart,x,y,z);
		
		double bmax = BMAX;

		if (i >= Sub.First) // subhaloes
			bmax = 2e-6;

		double B2 = p2(SphP[ipart].Bfld[0]) + p2(SphP[ipart].Bfld[1]) + 
					p2(SphP[ipart].Bfld[2]);

		if ((B2 > p2(bmax))) {

			double B = sqrt(B2);

			SphP[ipart].Bfld[0] *= bmax/B;
			SphP[ipart].Bfld[1] *= bmax/B;
			SphP[ipart].Bfld[2] *= bmax/B;

			cnt++;
		}
	}

	printf("Bfld of %llu particles limited to %g G\n", cnt, BMAX);

	return ;
}

