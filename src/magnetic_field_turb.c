#ifdef TURB_B_FIELD

#include "globals.h"
#include "tree.h"

// #include <iostream>
// #include <fstream.h>

// #define pi 3.14159

//#define p2(a) (a * a)
//#define p3(a) (a * a * a)
#define length3(a) (sqrt(a[0] * a[0] + a[1] * a[1] + a[2] + a[2]))

static double Kmin = 0, Kmax = 0;


static void allocate_grids(const int, double **, complex **, double **);
static void fill_fourier_grid(complex **, double **, const int);
static void normalise_bfld_grid(const int, double **);
static void divergence_clean_fourier_grid(const int, complex **, double **);
static void grid2particles_NGP(double **, const int);

static inline double powerspectrum(double k);
static inline double magnetic_field_model(float x, float y, float z);
static inline int Idx(const int, int, int, int); /* For Grid.B only */


void Bfld_from_turb_spectrum()
{
    int i = 0;

    // needs to be scaled down, no idea why
    Param.Bfld_Norm *= 0.1;

    complex *Bk[3] = {NULL}; /* holds B in k-space */
    double *B[3] = {NULL};    /* holds B in real-space */
    double *KVec[3] = {NULL}; /* holds k vector */

    const int nGrid = 2 * ceil(Param.Boxsize / Param.Bfld_Scale);
    Kmin = 2 * pi / Param.Boxsize;
    Kmax = 0.5 * Kmin * nGrid; /* Nyquist frequency */

    printf("Magnetic field: (Bonafede et al. 2010) \n"
           "   B0              = %g G\n"
           "   N               = %d \n"
           "   max scale       = %g kpc \n"
           "   min scale       = %g kpc\n"
           "   kmax            = %g /kpc \n"
           "   kmin            = %g /kpc \n"
           "   spectral idx    = %g  \n",
           Param.Bfld_Norm, nGrid, Param.Boxsize, Param.Boxsize / nGrid * 2,
           Kmax, Kmin, Param.Spectral_Index);

    allocate_grids(nGrid, B, Bk, KVec);

    fftw_plan forward_plan[3], backward_plan[3]; /* Init FFTW3 */
    for (i = 0; i < 3; i++)
    {
        forward_plan[i] = fftw_plan_dft_r2c_3d(nGrid, nGrid, nGrid, B[i], Bk[i], FFTW_ESTIMATE);
        backward_plan[i] = fftw_plan_dft_c2r_3d(nGrid, nGrid, nGrid, Bk[i], B[i], FFTW_ESTIMATE);
    }

    printf("filling FFT grid\n");

    fill_fourier_grid(Bk, KVec, nGrid); /* where the magic happens */

    printf("FFT back transform\n");

    fftw_execute(backward_plan[0]);
    fftw_execute(backward_plan[1]);
    fftw_execute(backward_plan[2]);

    printf("normalizing Bfield\n");

    normalise_bfld_grid(nGrid, B);

    printf("FFT B-field for divergence cleaning\n");

    fftw_execute(forward_plan[0]);
    fftw_execute(forward_plan[1]);
    fftw_execute(forward_plan[2]);

    printf("divergence cleaning\n");

    divergence_clean_fourier_grid(nGrid, Bk, KVec); /* Ruszkowski+ 2006 */

    printf("FFT transform to real space\n");

    fftw_execute(backward_plan[0]);
    fftw_execute(backward_plan[1]);
    fftw_execute(backward_plan[2]);

    printf("interpolate particles\n");

    grid2particles_NGP(B, nGrid);

    for (i = 0; i < 3; i++)
    {
        fftw_destroy_plan(forward_plan[i]);
        fftw_destroy_plan(backward_plan[i]);
        fftw_free(Bk[i]);
        fftw_free(B[i]);
        free(KVec[i]);
    }

    printf("done \n\n");
    fflush(stdout);

    return;
}

static void allocate_grids(const int nGrid, double **B, complex **Bk,
                           double **KVec)
{
    int i = 0;
    const int nComplex = nGrid * nGrid * (nGrid / 2 + 1);

    size_t nBytes = nComplex * sizeof(**Bk);
    for (i = 0; i < 3; i++)
        Bk[i] = fftw_malloc(nBytes);

    nBytes = p3(nGrid) * sizeof(**B);
    for (i = 0; i < 3; i++)
        B[i] = fftw_malloc(nBytes);

    nBytes = nComplex * sizeof(**KVec);
    for (i = 0; i < 3; i++)
        KVec[i] = Malloc(nBytes);

    return;
}

/* here the difficulty is to correctly define k in three dimensions
 * and get the symmetry in the k=0 plane correctly, so the FFT of the grid
 * is indeed real and not complex. */
static void fill_fourier_grid(complex **Bk, double **KVec, const int nGrid)
{
    int i = 0, j = 0, k = 0, idxconj = 0, icomp = 0;
    complex bk[3] = {0 + I * 0}; /* complex bfld in k-space */
    double kVec[3] = {0};

    /* Init GSL random number generator */
    gsl_rng_env_setup();
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *rng = gsl_rng_alloc(T);

    /* sample Pk */
    for (i = 0; i < nGrid; i++)
    {
        for (j = 0; j < nGrid; j++)
        {
            for (k = 0; k < nGrid / 2 + 1; k++)
            {

                if (!i && !j && !k)
                    continue;

                /* construct symmetry index */
                int iconj = nGrid - i;
                if (!i)
                    iconj = 0;

                int jconj = nGrid - j;
                if (!j)
                    jconj = 0;

                int kconj = (nGrid / 2 + 1) - k;
                if (!k)
                    kconj = 0;

                /* construct k vector */
                if (i <= nGrid / 2)
                    kVec[0] = i * Kmin;
                else
                    kVec[0] = -iconj * Kmin;

                if (j <= nGrid / 2)
                    kVec[1] = j * Kmin;
                else
                    kVec[1] = -jconj * Kmin;

                if (k <= (nGrid / 2 + 1)) /* always true */
                    kVec[2] = k * Kmin;
                else
                    kVec[2] = -kconj * Kmin;

                size_t idx = (nGrid * i + j) * (nGrid / 2 + 1) + k;

                KVec[0][idx] = kVec[0]; /* save for div(B) cleaning */
                KVec[1][idx] = kVec[1];
                KVec[2][idx] = kVec[2];

                double kMag = length3(kVec);
                if (kMag > Kmax) /* sphere in k-space */
                    continue;

                /* Sample P(k) - Box-Mueller Method */
                double pk = powerspectrum(kMag);

                for (icomp = 0; icomp < 3; icomp++)
                {
                    double amp = sqrt(-log(gsl_rng_uniform_pos(rng)) * pk);

                    double phase = 2 * pi * gsl_rng_uniform(rng);

                    bk[icomp] = amp * (cos(phase) + I * sin(phase));
                }

                /* 3D symmetries to get real field after inverse FFT */
                if (k > 0)
                {
                    Bk[0][idx] = bk[0];
                    Bk[1][idx] = bk[1];
                    Bk[2][idx] = bk[2];
                    /* FFTW3 doesn't store the conjugated ones for k>0 */
                }
                else
                { /* k=0 plane is a bitch ! */
                    if (!i)
                    {
                        if (j >= nGrid / 2)
                            continue;

                        Bk[0][idx] = bk[0]; /* j < nGrid/2 here */
                        Bk[1][idx] = bk[1];
                        Bk[2][idx] = bk[2];

                        idxconj = (nGrid * i + jconj) * (nGrid / 2 + 1) + k;
                        Bk[0][idxconj] = conj(bk[0]);
                        Bk[1][idxconj] = conj(bk[1]);
                        Bk[2][idxconj] = conj(bk[2]);
                    }
                    else
                    { /* i != 0, nastyness doesn't end */
                        if (i >= nGrid / 2)
                            continue;

                        Bk[0][idx] = bk[0]; /* i < nGrid/2 here */
                        Bk[1][idx] = bk[1];
                        Bk[2][idx] = bk[2];

                        idxconj = (nGrid * iconj + jconj) * (nGrid / 2 + 1) + k;
                        Bk[0][idxconj] = conj(bk[0]);
                        Bk[1][idxconj] = conj(bk[1]);
                        Bk[2][idxconj] = conj(bk[2]);
                    }
                } /* ... thank god its over */
            }
        }
    }

    gsl_rng_free(rng);

    return;
}

/* we use the global mean of the P(k), because we are to lazy to
 * think how P(k) should be correctly normalised using Parseval's theorem. */
static void normalise_bfld_grid(const int nGrid, double **B)
{
    int i, j, k;

    const int nTotal = nGrid * nGrid * (nGrid / 2 + 1);

    const float cell2kpc = Param.Boxsize / nGrid;

    double global_mean = 0;
    for (i = 0; i < nTotal; i++)
        global_mean += sqrt(p2(B[0][i]) + p2(B[1][i]) + p2(B[2][i]));
    global_mean /= nTotal;

    const int fftw3_norm = nGrid * nGrid * nGrid;

    const double normInv = 1.0 / global_mean / fftw3_norm;

    
    for (i = 0; i < nGrid; i++)
    {
        for (j = 0; j < nGrid; j++)
        {
            for (k = 0; k < nGrid; k++)
            {
                float x = (i + 0.5 * (1 - nGrid)) * cell2kpc;
                float y = (j + 0.5 * (1 - nGrid)) * cell2kpc;
                float z = (k + 0.5 * (1 - nGrid)) * cell2kpc;

                double B_max = 0;
                double B_model = 0.0;

                for (int halo = 0; halo < Param.Nhalos; halo++)
                {

                    float dx = x - Halo[halo].D_CoM[0] - 0.5*Param.Boxsize,
                          dy = y - Halo[halo].D_CoM[1] - 0.5*Param.Boxsize,
                          dz = z - Halo[halo].D_CoM[2] - 0.5*Param.Boxsize;

                    double r2 = dx * dx + dy * dy + dz * dz;

                    double rho_i = Gas_Density_Profile(sqrt(r2), halo);

                    B_model = pow(rho_i / Halo[halo].Rho0, Param.Bfld_Eta);

                    if (B_model > B_max)
                        B_max = B_model;
                }

                double bNorm = B_model * normInv;

                int idx = Idx(nGrid, i, j, k);
                B[0][idx] *= bNorm;
                B[1][idx] *= bNorm;
                B[2][idx] *= bNorm;
            }
        }
    }

    return;
}

/* e.g. Balsara 1996, I guess one has to just try that one out. */
static void divergence_clean_fourier_grid(const int nGrid, complex **Bk,
                                          double **kVec)
{

    for (int i = 0; i < nGrid; i++)
    {
        for (int j = 0; j < nGrid; j++)
        {
            for (int k = 0; k < nGrid / 2 + 1; k++)
            {
                size_t idx = (nGrid * i + j) * (nGrid / 2 + 1) + k;

                complex Bx = Bk[0][idx];
                complex By = Bk[1][idx];
                complex Bz = Bk[2][idx];

                float kx = kVec[0][idx];
                float ky = kVec[1][idx];
                float kz = kVec[2][idx];

                float k2inv = 1.0 / (kx * kx + ky * ky + kz * kz);

                Bk[0][idx] = Bx * (1 - kx * kx * k2inv) - By * kx * ky * k2inv - Bz * kx * kz * k2inv;
                Bk[1][idx] = -Bx * ky * kx * k2inv + By * (1 - ky * ky * k2inv) - Bz * ky * kz * k2inv;
                Bk[2][idx] = -Bx * kz * kx * k2inv - By * kz * ky * k2inv - Bz * (1 - kz * kz * k2inv);
            }
        }
    }

    Bk[0][0] = Bk[1][0] = Bk[2][0] = 0;

    return;
}

/* Nearest Grid Point interpolation (Hockney & Eastwood) */
static void grid2particles_NGP(double **B, const int nGrid)
{

    double cellsize = Param.Boxsize / nGrid;

    for (int ipart = 0; ipart < Param.Npart[0]; ipart++)
    {
        float u = P[ipart].Pos[0] / cellsize; /* center */
        float v = P[ipart].Pos[1] / cellsize;
        float w = P[ipart].Pos[2] / cellsize;

        while (u < 0) /* ensure perdiodicity */
            u += nGrid;
        while (v < 0)
            v += nGrid;
        while (w < 0)
            w += nGrid;

        while (u >= nGrid)
            u -= nGrid;
        while (v >= nGrid)
            v -= nGrid;
        while (w >= nGrid)
            w -= nGrid;

        int i = floor(u); /* NGP */
        int j = floor(v);
        int k = floor(w);

        SphP[ipart].Bfld[0] = B[0][Idx(nGrid, i, j, k)];
        SphP[ipart].Bfld[1] = B[1][Idx(nGrid, i, j, k)];
        SphP[ipart].Bfld[2] = B[2][Idx(nGrid, i, j, k)];
    }

    return;
}

/* unnormed power-law, boring isn't it ?  */
static inline double powerspectrum(double k)
{
    return pow(k, Param.Spectral_Index);
}

/* Index of the real-space grid of doubles */
static inline int Idx(const int nGrid, int ix, int iy, int iz)
{
    return (nGrid * ix + iy) * nGrid + iz;
}

#endif