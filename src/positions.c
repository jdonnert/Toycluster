#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sort.h>

#include "globals.h"

#define NTABLE 2048    // length of interpolation table of M(r) 

/* Cum. mass profile is inverted from a table */
static void fill_mass_profile_table(int);
static double invert_mass_profile(double);
static void sort_particles(int *, const size_t );

static void sample_DM_particles(const int);
static void sample_Gas_particles(const int);

/* Cummulative Mass table interpolation */
static double dfdr_table[NTABLE], offset_table[NTABLE], 
       Mass_profile_table[NTABLE];

/*
 * Positions are sampled around 0 ! Haloes are moved into position later
 */

void Make_positions()
{
    printf("Sampling positions "); fflush(stdout);

    for (int i = 0; i < Param.Nhalos; i++) {

		if (i < Sub.First)
	        printf("<%d> ", i); 
		else
			printf(".");

		fflush(stdout);

		sample_DM_particles(i);
		
		sample_Gas_particles(i);
   	} 

    printf(" done\n");

    return;
}

static void sample_DM_particles(const int i)
{
	const double dCoM[3] = {Halo[i].D_CoM[0],Halo[i].D_CoM[1],Halo[i].D_CoM[2]};
	const double qmax = Halo[i].MassCorrFac; 

	#pragma omp parallel for
    for (int ipart = 0; ipart < Halo[i].Npart[1]; ipart++) { // DM

		for (;;) { // Hernquist halo M(<R) inverted

			float theta = acos(2 *  erand48(Omp.Seed) - 1);
           	float phi = 2*pi * erand48(Omp.Seed);

           	float sin_theta = sin(theta);
           	float cos_theta = cos(theta);

           	float sin_phi = sin(phi);
           	float cos_phi = cos(phi);

           	double sqrt_q =  sqrt(erand48(Omp.Seed) * qmax);  
           	double r = Halo[i].A_hernq * sqrt_q / (1-sqrt_q);

           	double x = r * sin_theta * cos_phi;
           	double y = r * sin_theta * sin_phi;
           	double z = r * cos_theta;

			if (i != Halo_containing(1, x+dCoM[0], y+dCoM[1], z+dCoM[2])) 
            	continue; // draw another one 
			
			Halo[i].DM[ipart].Pos[0] = (float) x;
           	Halo[i].DM[ipart].Pos[1] = (float) y;
           	Halo[i].DM[ipart].Pos[2] = (float) z;

			Halo[i].DM[ipart].Type = 1;

			break;
       	}
   	} // ipart

	return ;
}

static void sample_Gas_particles(const int i)
{
	const double dCoM[3] ={Halo[i].D_CoM[0],Halo[i].D_CoM[1],Halo[i].D_CoM[2]};
	const double boxhalf = Param.Boxsize/2;

	fill_mass_profile_table(i); 

	#pragma omp parallel for
   	for (size_t ipart = 0; ipart < Halo[i].Npart[0]; ipart++) {
		
		for (;;) { 

			double theta = acos(2 *  erand48(Omp.Seed) - 1);
        	double phi = 2*pi * erand48(Omp.Seed);

           	double q = erand48(Omp.Seed);  
           	double r = invert_mass_profile(q); // beta model, alpha = 2/3
   
           	double x = r * sin(theta) * cos(phi);
           	double y = r * sin(theta) * sin(phi);
            double z = r * cos(theta);
			
			if (i != Halo_containing(0, x+dCoM[0], y+dCoM[1], z+dCoM[2])) 
            	continue; // draw another one 

			if ((x < -boxhalf) || (x > boxhalf))	
				continue;
			if ((y < -boxhalf) || (y > boxhalf))	
				continue;
			if ((z < -boxhalf) || (z > boxhalf))	
				continue;

            Halo[i].Gas[ipart].Pos[0] = (float) x;
   	        Halo[i].Gas[ipart].Pos[1] = (float) y;
       	    Halo[i].Gas[ipart].Pos[2] = (float) z;

			Halo[i].Gas[ipart].Type = 0;

           	break;
   		}
   	} //ipart

	return ;
}

/* 
 * For gas particles we need to invert M(<r)/M
 * which is not analytical for this model.
 * instead we write M(<r)/M into a table and 
 * invert numerically with linear interpolation. 
 */

static double invert_mass_profile(const double q)
{
    int i = 0;

    for (i = 0; i < NTABLE-1; i++) // exit with i<=NTABLE-1 to avoid overflow
        if (Mass_profile_table[i+1] > q)
            break;

    return (q-offset_table[i])/dfdr_table[i];
}

static void fill_mass_profile_table(const int i)
{
	double step = log10(Halo[i].R_Sample[0])/NTABLE;

    double r[NTABLE] = { 0 };
	
	size_t nBytes = NTABLE * sizeof(double);

	memset(Mass_profile_table, 0, nBytes);
	memset(dfdr_table, 0, nBytes);
	memset(offset_table, 0, nBytes);

    for (int j = 0; j < NTABLE; j++) { // mass profile

        r[j] = pow(10,j*step) - 1;

	    Mass_profile_table[j] = Mass_profile(r[j], Halo[i].Rho0, 
					Halo[i].Rcore, Halo[i].Rcut, Halo[i].Have_Cuspy);
    }

	for (int j = 0; j < NTABLE; j++)
		Mass_profile_table[j] /= Mass_profile_table[NTABLE-1];
	
	for (int j = 1; j < NTABLE; j++) { // linear interpolation

        dfdr_table[j] = (Mass_profile_table[j] - Mass_profile_table[j-1])
            /(r[j] - r[j-1]);
        
        offset_table[j] =  Mass_profile_table[j] - dfdr_table[j]*r[j];
    }

    return;
}

/* 
 * Because sampling depends on boxsize
 * we print the mass in r200 
 * We show the main halo including substructure
 * if REPORTSUBHALOS is set, sampling of all subhalos is shown 
 */

void Show_mass_in_r200()
{
    const double mSph = Param.Mpart[0]*Unit.Mass/Msol2cgs;
    const double mDM = Param.Mpart[1]*Unit.Mass/Msol2cgs;
    const float boxhalf = Param.Boxsize / 2;

	int nShow = Sub.First;

#ifdef REPORTSUBHALOS 
	nShow = Param.Nhalos;
#endif

    for (int i = 0; i < nShow; i++) {

		int npart = Halo[i].Npart[0];

		if (i == 0)
			for (int j = Sub.First; j < Param.Nhalos; j++)
				npart += Halo[j].Npart[0];

        float rVir2 = p2(Halo[i].R200);

        size_t nSph = 0;

        for (size_t ipart = 0; ipart < npart; ipart++) {
            
            float dx = Halo[i].Gas[ipart].Pos[0] - Halo[i].D_CoM[0] - boxhalf;
            float dy = Halo[i].Gas[ipart].Pos[1] - Halo[i].D_CoM[1] - boxhalf;
            float dz = Halo[i].Gas[ipart].Pos[2] - Halo[i].D_CoM[2] - boxhalf;

            float r2 = dx*dx + dy*dy + dz*dz;
            
            if (r2 < rVir2)
                nSph++;
        }
    
		npart = Halo[i].Npart[1];

		if (!i)
			for (int j = Sub.First; j < Param.Nhalos; j++)
				npart += Halo[j].Npart[1];

        size_t nDM = 0;
        for (size_t ipart = 0; ipart < npart; ipart++) {

            float dx = Halo[i].DM[ipart].Pos[0] - Halo[i].D_CoM[0] - boxhalf;
            float dy = Halo[i].DM[ipart].Pos[1] - Halo[i].D_CoM[1] - boxhalf;
            float dz = Halo[i].DM[ipart].Pos[2] - Halo[i].D_CoM[2] - boxhalf;

            float r2 = dx*dx + dy*dy + dz*dz;

            if (r2 < rVir2)
                nDM++;
        }
        
        double mass200 = nSph*mSph + nDM*mDM;

        printf("\nSampling of Halo <%d> (r200 = %g kpc):\n"
        	"   Gas Mass in R200    = %g Msol \n"
            "   DM Mass in R200     = %g Msol \n"
            "   Total Mass in R200  = %g Msol \n"
            "   External Gas Mass   = %g Msol \n"
            "   External DM  Mass   = %g Msol \n"
            "   Total External Mass = %g Msol \n"
            "   Effective bf in r200= %g \n",
            i, sqrt(rVir2), nSph*mSph, nDM*mDM, mass200, 
			(Halo[i].Npart[0]-nSph)*mSph, (Halo[i].Npart[1]-nDM)*mDM,
            (Halo[i].Npart[0] - nSph)*mSph + (Halo[i].Npart[1] - nDM)*mDM, 
            nSph*mSph/ (nDM*mDM));
    }

    printf("\n");

    return;
}

/* center on CoM */
void center_positions()
{
    const double mGas = Param.Mpart[0];
    const double mDm = Param.Mpart[1];

    for (int i = 0; i < Param.Nhalos; i++) {

        double sum[3] = { 0 };

        for(size_t ipart=0; ipart<Halo[i].Npart[0]; ipart++) {

            sum[0] += mGas * Halo[i].Gas[ipart].Pos[0];
            sum[1] += mGas * Halo[i].Gas[ipart].Pos[1];
            sum[2] += mGas * Halo[i].Gas[ipart].Pos[2];
        }

        for(size_t ipart=0; ipart<Halo[i].Npart[1]; ipart++) {

            sum[0] += mDm * Halo[i].DM[ipart].Pos[0];
            sum[1] += mDm * Halo[i].DM[ipart].Pos[1];
            sum[2] += mDm * Halo[i].DM[ipart].Pos[2];
        }
        
        sum[0] /= mGas * Halo[i].Npart[0] + mDm * Halo[i].Npart[1];
        sum[1] /= mGas * Halo[i].Npart[0] + mDm * Halo[i].Npart[1];
        sum[2] /= mGas * Halo[i].Npart[0] + mDm * Halo[i].Npart[1];

        for(size_t ipart=0; ipart<Halo[i].Npart[0]; ipart++) {

            Halo[i].Gas[ipart].Pos[0] -= sum[0];
            Halo[i].Gas[ipart].Pos[1] -= sum[1];
            Halo[i].Gas[ipart].Pos[2] -= sum[2];
        }

        for(size_t ipart=0; ipart<Halo[i].Npart[1]; ipart++) {

            Halo[i].DM[ipart].Pos[0] -= sum[0];
            Halo[i].DM[ipart].Pos[1] -= sum[1];
            Halo[i].DM[ipart].Pos[2] -= sum[2];
        }
    }

    return;
}

void Reassign_particles_to_halos()
{
	const float boxhalf = 0.5 * Param.Boxsize;
	int *haloID = Malloc( Param.Npart[0]*sizeof(*haloID) );

	size_t npart[Param.Nhalos];
	memset(npart, 0, Param.Nhalos * sizeof(*npart));

   	for (size_t ipart = 0; ipart < Param.Npart[0]; ipart++) {
   	
		float x = P[ipart].Pos[0] - boxhalf,
			  y = P[ipart].Pos[1] - boxhalf,
			  z = P[ipart].Pos[2] - boxhalf;

		int i = Halo_containing(P[ipart].Type, x, y, z);
		
		haloID[ipart] = i;

		npart[i]++;
	}

	sort_particles(haloID, Param.Npart[0]);

	Free(haloID);

	Sub.Ntotal = Sub.Npart[1]; // update particle numbers
	Sub.Npart[0] = 0;

	for (int i = 0; i < Param.Nhalos; i++) { 

		Halo[i].Npart[0] = npart[i];

		Halo[i].Ntotal = Halo[i].Npart[0] + Halo[i].Npart[1];

		if (i >= Sub.First) {
		
			Sub.Ntotal += npart[i];
			
			Sub.Npart[0] += npart[i];
		}
	}
	
	int iGas = Halo[0].Npart[0]; // update pointers

	for (int i = 1; i < Param.Nhalos; i++) {

		Halo[i].Gas = &(P[iGas]);
		Halo[i].SphP = &(SphP[iGas]);

		iGas += Halo[i].Npart[0];
	}

	printf("Particle Distribution after Relaxation :\n"
			"   Main     %8lld   %8lld   %8lld  \n",
			Halo[0].Ntotal, Halo[0].Npart[0], Halo[0].Npart[1]);
	if (Param.Mass_Ratio)
		printf("   Bullet   %8lld   %8lld   %8lld  \n", 
				Halo[1].Ntotal, Halo[1].Npart[0], Halo[1].Npart[1]);
	
#ifdef SUBSTRUCTURE
			printf("   Subhalos %8d   %8d   %8d \n",
				Sub.Ntotal, Sub.Npart[0], Sub.Npart[1]);
#endif

	return ;
}

/* check for collision between position xyz and clusters
 * via maximum of density for gas particles and sampling radius for DM */
int Halo_containing(const int type, const float x, const float y, const float z)
{
	const double boxsize = Param.Boxsize;

	if ((x > boxsize || y > boxsize || z > boxsize)) 
		return -1;

	int i = 0;

	if (type > 0) { // DM

		float r = sqrt(p2(x - Halo[1].D_CoM[0]) + p2(y - Halo[1].D_CoM[1])  
				 + p2(z - Halo[1].D_CoM[2]));

		if ( (r < Halo[1].R_Sample[1]) && (x > 0) )
			i = 1;

		for (int j = Sub.First; j < Param.Nhalos; j++) { // 0 is std

   	   		float r = sqrt(p2(x - Halo[j].D_CoM[0]) + p2(y - Halo[j].D_CoM[1])  
				 + p2(z - Halo[j].D_CoM[2]));

			if (r < Halo[j].R_Sample[1]) {

				i = j;  // we are in the sampling radius

				break;
			}
		}

	} else { // SPH
	
		double rho_max = 0;

		for (int j = 0; j < Param.Nhalos; j++) { // 0 is std

			if (Halo[j].Is_Stripped)
				continue; // stripped halos don't contain SPH particles !

   	   		float r = sqrt(p2(x - Halo[j].D_CoM[0]) + p2(y - Halo[j].D_CoM[1])  
				 + p2(z - Halo[j].D_CoM[2]));

			double rho_gas = Gas_density_profile(r, Halo[j].Rho0, 
							Halo[j].Rcore, Halo[j].Rcut,Halo[j].Have_Cuspy);

       		if ( (rho_gas > rho_max) && (r < Halo[j].R_Sample[0]) ) {
			
				i = j;

				rho_max = rho_gas;
			}
		}
	}
	
	return i;
}

/* testing */
int compare_int(const void * a, const void *b) 
{
	const int *x = (const int*)a;
	const int *y = (const int*)b;

    return (*x > *y) - (*x < *y);
}


/* memory efficient out of place sorting of both particle structures 
 * Knowing where idx starts in memory we can reconstruct ipart */
static void sort_particles(int *ids, const size_t nPart) 
{
    size_t *idx = Malloc(nPart*sizeof(*idx));

	Qsort_Index(Omp.NThreads, idx, ids, nPart, sizeof(*ids), &compare_int);
	
    for (size_t i = 0; i < nPart; i++) {

        if (idx[i] == i)
            continue;

		size_t dest = i;

        struct ParticleData Ptmp = P[dest];
        struct GasParticleData SphPtmp = SphP[dest];
		size_t src = idx[i];

        for (;;) {

			memcpy(&P[dest], &P[src], sizeof(*P));
            memcpy(&SphP[dest], &SphP[src], sizeof(*SphP));
            idx[dest] = dest;

			dest = src;
			src = idx[dest];

            if (src == i) 
                break;
        }

		memcpy(&P[dest], &Ptmp, sizeof(*P));
        memcpy(&SphP[dest], &SphPtmp, sizeof(*SphP));
        idx[dest] = dest;
    }
 
    Free(idx);

    return;
}
