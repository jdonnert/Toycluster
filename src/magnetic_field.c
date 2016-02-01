#include "globals.h"
#include "tree.h"

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

	for (int i = 0; i < Param.Nhalos; i++) {

		double rho0 =  Halo[i].Rho0;

		if (rho0 == 0)
			continue;
			
		#pragma omp parallel for 
    	for (size_t ipart = 0; ipart < Param.Npart[0]; ipart++) {
		
			float x = P[ipart].Pos[0] - boxhalf,
				  y = P[ipart].Pos[1] - boxhalf,
				  z = P[ipart].Pos[2] - boxhalf;

			int j = Halo_containing(0, x, y, z);

			if (j != i)
				continue;

			double rho_local = SphP[ipart].Rho;

			double rho_norm =  rho_local / rho0;
		
			double A = pow(rho_norm, Param.Bfld_Eta);
		
			SphP[ipart].Apot[0] = (float) A;
			SphP[ipart].Apot[1] = (float) A;
	    	SphP[ipart].Apot[2] = (float) A;
		}
	}

    return ;  
}

static void normalise_magnetic_field()
{
	const float boxhalf = 0.5 * Param.Boxsize;

 	double max_B2 = 0;

	for (int i = 0; i < Param.Nhalos; i++) {

		#pragma omp parallel for 
	    for (size_t ipart = 0; ipart < Param.Npart[0]; ipart++) {
		
			float x = P[ipart].Pos[0] - boxhalf,
				  y = P[ipart].Pos[1] - boxhalf,
				  z = P[ipart].Pos[2] - boxhalf;

			int j = Halo_containing(P[ipart].Type, x, y, z);

			if (j != i)
				continue;

			double bfld2 = p2(SphP[ipart].Bfld[0]) 
				+ p2(SphP[ipart].Bfld[1]) 
				+ p2(SphP[ipart].Bfld[2]);

			max_B2 = fmax(max_B2, bfld2);
		}

		double max_bfld = sqrt(max_B2);

    	double norm = 1.3*Param.Bfld_Norm/max_bfld;

		if (Param.Cuspy != 0)
			norm *= 6;

		#pragma omp parallel for 
    	for (size_t ipart = 0; ipart < Param.Npart[0]; ipart++) {
		
			float x = P[ipart].Pos[0] - boxhalf,
				  y = P[ipart].Pos[1] - boxhalf,
				  z = P[ipart].Pos[2] - boxhalf;

			int j = Halo_containing(P[ipart].Type, x, y, z);

			if (j != i)
				continue;

			SphP[ipart].Bfld[0] *= norm;
			SphP[ipart].Bfld[1] *= norm;
			SphP[ipart].Bfld[2] *= norm;
		}
	}

	return ;
}

