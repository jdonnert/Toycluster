#include "globals.h"

/* Set up a box with SPH particles, then WVT relax the particles */

int main(int argc, char *argv[])
{
    printf("--- This is %s, Version %s ---\n", PROG_NAME, VERSION);

#pragma omp parallel
    {

        Omp.ThreadID = omp_get_thread_num();
        Omp.NThreads = omp_get_num_threads();
        Omp.Seed[2] = 14041981 * (Omp.ThreadID + 1);
        erand48(Omp.Seed); // remove leading 0

        if (Omp.ThreadID == 0)
            printf("Running with %d Threads\n", Omp.NThreads);

    } // omp parallel

    Assert(argc == 2, "Usage : ./wvtbox $parameterfile\n");

    Read_param_file(argv[1]);

    printf("Parameters:\n  Ntotal  = %d\n  Boxsize = %g\n\n",
            Param.Ntotal, Param.Boxsize);

    Set_units();

    Setup();

    Make_positions();

    Make_IDs();

    Shift_Origin(); // from here onwards: 0 < Pos < Boxsize

    // Find_sph_quantities();

    // char wvt_filename[CHARBUFSIZE] = "test";
    // sprintf(wvt_filename, "test_%d", 99);
    // Write_positions(wvt_filename);

    /* for (int ipart = 0; ipart < Box.Npart[0]; ipart++) {
        printf("ipart=%d  |  %g %g %g %g\n", ipart, Box.Gas[ipart].Pos[0],
            Box.Gas[ipart].Pos[1], Box.Gas[ipart].Pos[2],
            sqrt(p2(Box.Gas[ipart].Pos[0])+p2(Box.Gas[ipart].Pos[1])
                +p2(Box.Gas[ipart].Pos[2])));
    } */

    Regularise_sph_particles();

//    Find_sph_quantities();
//
//    Make_velocities();
//
//    Write_output();

    return EXIT_SUCCESS ;
}

