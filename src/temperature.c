#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#include "globals.h"

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

		Setup_Profiles(i);

		#pragma omp parallel for
        for (int ipart = 0; ipart < Halo[i].Npart[0]; ipart++) {

            float dx = Halo[i].Gas[ipart].Pos[0] - Halo[i].D_CoM[0] - boxhalf;
            float dy = Halo[i].Gas[ipart].Pos[1] - Halo[i].D_CoM[1] - boxhalf;
            float dz = Halo[i].Gas[ipart].Pos[2] - Halo[i].D_CoM[2] - boxhalf;

            double r = sqrt(dx*dx + dy*dy + dz*dz);

			Halo[i].SphP[ipart].U = Internal_Energy_Profile(i, r);
		}
    }
    
	printf("done\n\n"); fflush(stdout);

    return;
}
