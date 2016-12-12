#include "globals.h"

void Setup()
{
    int nGas = Param.Ntotal;

    Param.Mpart[0] = 0.0025;

    Box.Npart[0] = (int) nGas;

    for (int i = 0; i < 6; i++) {
        Param.Npart[i] = Box.Npart[i];
        printf("Param.Npart[%d] = %d\n", i, Box.Npart[i]);
        printf("Param.Mpart[%d] = %g\n", i, Param.Mpart[i]);
    }
    printf("\n");

    /* allocate particles */
    size_t nBytes = Param.Ntotal * sizeof(*P);
    P = Malloc(nBytes);
    memset(P, 0, nBytes);

    nBytes = Param.Npart[0] * sizeof(*SphP);
    SphP = Malloc(nBytes);
    memset(SphP, 0, nBytes);

    /* set access pointers */
    Box.Gas = &(P[0]);
    Box.SphP = &(SphP[0]);

    Box.D_CoM[0] = 0;
    Box.D_CoM[1] = 0;
    Box.D_CoM[2] = 0;

    Box.R_Sample[1] = Param.Boxsize/2;
    Box.R_Sample[0] = Param.Boxsize/2;

    return;
}

/* center the simulation in the periodic box */
void Shift_Origin()
{
    const float boxsize = Param.Boxsize;
    const float boxHalf = boxsize / 2;

    printf("Shift Origin ..."); fflush(stdout);

    float dx = Box.D_CoM[0];
    float dy = Box.D_CoM[1];
    float dz = Box.D_CoM[2];

#pragma omp parallel for  // Gas
    for (size_t ipart=0; ipart<Box.Npart[0]; ipart++) {

        Box.Gas[ipart].Pos[0] += dx;
        Box.Gas[ipart].Pos[1] += dy;
        Box.Gas[ipart].Pos[2] += dz;
    }

    printf(" done \n\n");

#pragma omp parallel for // shift 0 to bottom left corner, wrap around box
    for (size_t ipart=0; ipart<Param.Ntotal; ipart++) {

        P[ipart].Pos[0] += boxHalf;
        P[ipart].Pos[1] += boxHalf;
        P[ipart].Pos[2] += boxHalf;

        while (P[ipart].Pos[0] > boxsize)
            P[ipart].Pos[0] -= boxsize;

        while (P[ipart].Pos[0] < 0)
            P[ipart].Pos[0] += boxsize;

        while (P[ipart].Pos[1] > boxsize)
            P[ipart].Pos[1] -= boxsize;

        while (P[ipart].Pos[1] < 0)
            P[ipart].Pos[1] += boxsize;

        while (P[ipart].Pos[2] > boxsize)
            P[ipart].Pos[2] -= boxsize;

        while (P[ipart].Pos[2] < 0)
            P[ipart].Pos[2] += boxsize;
    }

    return ;
}

