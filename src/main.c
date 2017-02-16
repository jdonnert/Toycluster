#include "globals.h"

/* Make IC of a cluster collision
 * (Mastropietro+ 2008, Donnert 2014, Donnert+ 2017) */

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

    Assert(argc == 2, "Usage : ./Toycluster $parameterfile\n");

    Read_param_file(argv[1]);

    Set_units();

    Set_cosmology();

    Setup();

    Make_positions();

    Make_IDs();

    Shift_Origin(); // from here onwards: 0 < Pos < Boxsize

    Show_mass_in_r200();

    if (Halo[0].bf) { // make SPH ICM

        Regularise_sph_particles();

        Find_sph_quantities();

        Make_magnetic_field();

        Reassign_particles_to_halos();

        Show_mass_in_r200();

        Make_temperatures();
    }

    Make_velocities();

    Apply_kinematics();

    Write_output();

    return EXIT_SUCCESS ;
}
