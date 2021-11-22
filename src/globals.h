#ifndef GLOBALS_H
#define GLOBALS_H

#define PROG_NAME "Toycluster_IC"
#define VERSION "2.0"

/* C std lib */
#include <stdlib.h>             // system
#include <stdio.h>
#include <stdint.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <complex.h>

/* GNU Scientifc Library */
#include <gsl/gsl_math.h>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_rng.h>

/* FFTW3 */
#include <complex.h>
#include <fftw3.h>

#include <omp.h>

/* Global Prototypes */
#include "macro.h"
#include "peano.h"
#include "cosmo.h"
#include "proto.h"

/* Code parameters */

#ifndef BETA
#define BETA (2.0/3.0) // Donnert 2014
#endif

#define CHARBUFSIZE 512        // For any char buffer !
#define MAXTAGS 300            // In parameter file

#ifdef SPH_CUBIC_SPLINE

#define DESNNGB  50        // SPH kernel weighted number of neighbours
#define NNGBDEV 0.05       // error tolerance in SPH kernel weighted neighb.
#define NGBMAX (DESNNGB*8)  // size of neighbour list

#else

#define DESNNGB  295        // SPH kernel weighted number of neighbours
#define NNGBDEV 0.05       // error tolerance in SPH kernel weighted neighb.
#define NGBMAX (DESNNGB*8)  // size of neighbour list

#endif // SPH_CUBIC_SPLINE


#define R200_TO_RMAX_RATIO 3.75 // this fits Y_M500 correlation
#define MAXHALOS 4096      // maximum number of subhalos
#define ZERO_ENERGY_ORBIT_FRACTION_SUB 1 // substructure vel fraction

/* mathematical constants */
#define pi             M_PI
#define sqrt2        M_SQRT2
#define sqrt3       1.73205080756887719
#define fourpithird 4.18879032135009765

/* physical constants cgs */
#define c            GSL_CONST_CGSM_SPEED_OF_LIGHT
#define k_B         GSL_CONST_CGSM_BOLTZMANN
#define m_p         GSL_CONST_CGSM_MASS_PROTON
#define m_e            GSL_CONST_CGSM_MASS_ELECTRON
#define Grav        GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT

/* unit conversions */
#define Msol2cgs    (1.98892e33)
#define kpc2cgs     (3.08568025e21)
#define K2eV        (1.5*8.617343e-5)
#define DEG2RAD        (pi / 180)

/* chemistry */
#define H_frac        0.76                    /* Hydrogen fraction */
#define He_frac        (1.0-H_frac)            /* Helium fraction */
#define u_mol        (4.0/(5.0*H_frac+3.0))    /* Mean mol. weight in hydr. mass */
#define n2ne         ((H_frac+0.5*He_frac)/(2.0*H_frac+0.75*He_frac))
#define yHelium        He_frac / (4.0 *H_frac)
#define mean_mol_weight (1.0+4.0*yHelium)/(1.0+3.0*yHelium+1.0)
#define adiabatic_index (5./3.)


extern struct OpenMP_infos{
    int NThreads;          // Number of openMP threads
    int ThreadID;          // Thread ID of this thread
    unsigned short Seed[3]; // random number seed: erand48(Omp.Seed)
} Omp;
#pragma omp threadprivate(Omp)

extern struct Parameters{
    char Output_File[CHARBUFSIZE];
    long long Ntotal;
    long long Npart[6];
    double Mtotal;
    double Mtot200;                    // total mass in R200 of system
    double Mass_Ratio;
    double Impact_Param;
    double Mpart[6];                // particle masses
    double Redshift;                // cluster redshift for scaling relations
    int Cuspy;                      // Cuspy Parameter, set in Binary
    double Bfld_Norm;               // B0
    double Bfld_Eta;                // Bfld scaling (Bonafede 10)
    double Boxsize;
    double VelMerger[2];            // Merging velocity
    int Nhalos;                     // Number of halos, incl. substructure
    double GravSofteningLength;
    double Zero_Energy_Orbit_Fraction;
#ifdef ADD_THIRD_SUBHALO
    double SubFirstMass;
    double SubFirstPos[3];
    double SubFirstVel[3];
#endif
#ifdef DOUBLE_BETA_COOL_CORES
    double Rho0_Fac;
    double Rc_Fac;
#endif
#ifdef TURB_B_FIELD
    double Bfld_Scale;
    double Spectral_Index;
#endif
} Param;

extern struct SubhaloData {
    int First;
    int Ntotal;
    int Npart[6];
    double Mtotal;
    double MassFraction;
    int Nhalos;
} Sub;

extern struct HaloProperties {
    long long Ntotal;               // total npart this cluster
    long long Npart[6];             // npart of particle type
    int Have_Cuspy;                 // Is this cluster cored or cuspy
    int Is_Stripped;
    double Mtotal;                  // Mass of this cluster
    double Mtotal200;               // Total Mass inside R200
    double Mass[6];                 // Mass in 6 particle types at R_Sample
    double Mass200[6];              // Mass inside R200
    double MassCorrFac;             // Correct DM profile for Rsample != infty
    double Rho0_DM;                 // Norm of NFW profile
    double C_nfw;                   // NFW profile concentration param
    double Rs;                      // NFW eq. scale radius
    double R200;                    // Virial Radius
    double R500;                    // Observational virial Radius
    double A_hernq;                 // Hernquist parameter <- C_nfw
    double Rho0;                    // Density normalisation
    double Beta;                    // Beta-model exponent
    double bf;                      // Halo-specific Baryon fraction
    double Rcut_R200_Ratio;         // Halo-specific rcut-to-r200 ratio
    double Rcore;                   // Beta-model core radius
    double Bf_eff;                  // Effective Baryon Fraction in r500
    double D_CoM[3];                // Distance from Center of Mass in x
    double BulkVel[3];              // Velocity rel to Center of Mass
    double R_Sample[2];             // DM Sampling size
    double Rcut;                    // SPH Sampling size
    double TempOffset;              // T of the parent cluster
    struct ParticleData *DM;        // DM Particle Data in P
    struct ParticleData *Gas;       // Gas Particle Data in P
    struct GasParticleData *SphP;   // Gas Particle Data in SphP
} Halo[MAXHALOS];

extern struct ParticleData{
    float Pos[3];
    float Vel[3];
    int32_t ID;
    int Type;
    peanoKey Key;
    int Tree_Parent;
} *P;

extern struct GasParticleData {
    float U;
    float Rho;
    float Hsml;
    float VarHsmlFac;
    float Bfld[3];
    float Apot[3];
    float ID;
    float Rho_Model;
    float Rs[3];
} *SphP;

/* code units */
extern struct Units{
    double Length;
    double Mass;
    double Vel;
    double Time;
    double Density;
    double Energy;
} Unit;

double G; // gravitational constant in code units

static const double Infinity = 1e25; // global boundaries for profiles
static const double Zero = 1;

#endif
