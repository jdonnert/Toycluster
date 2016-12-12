#include "globals.h"

static void sample_Gas_particles();

void Make_positions()
{
    printf("Sampling positions ..."); fflush(stdout);

    sample_Gas_particles();

    printf(" done\n");

    return;
}

static void sample_Gas_particles()
{
    const double boxsize = Param.Boxsize;
    const double boxhalf = boxsize*0.5;

    // printf("Box.Npart[0]=%d\n",Box.Npart[0]);

#pragma omp parallel for
    for (int ipart = 0; ipart < Box.Npart[0]; ipart++) {
        for (;;) {
            double x = boxsize*erand48(Omp.Seed) - boxhalf;
            double y = boxsize*erand48(Omp.Seed) - boxhalf;
            double z = boxsize*erand48(Omp.Seed) - boxhalf;

            if ((x < -boxsize) || (x > boxsize))
                continue;
            if ((y < -boxsize) || (y > boxsize))
                continue;
            if ((z < -boxsize) || (z > boxsize))
                continue;

            // printf("ipart = %zu; (x, y, z) = (%g, %g, %g)\n", ipart, x, y, z);

            Box.Gas[ipart].Pos[0] = (float) x;
            Box.Gas[ipart].Pos[1] = (float) y;
            Box.Gas[ipart].Pos[2] = (float) z;

            Box.Gas[ipart].Type = 0;

            Box.Gas[ipart].Vel[0] = 0;
            Box.Gas[ipart].Vel[1] = 0;
            Box.Gas[ipart].Vel[2] = 0;

            Box.SphP[ipart].Rho = Gas_Density_Profile(x, y, z);
            break;
        }
        // printf("%g %g %g\n", Box.Gas[ipart].Pos[0], Box.Gas[ipart].Pos[1], Box.Gas[ipart].Pos[2]);
    }

    return;
}
