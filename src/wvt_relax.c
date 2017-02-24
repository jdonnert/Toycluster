#include "globals.h"
#include "tree.h"

#define WVTNNGB DESNNGB // 145 for WC2 that equals WC6 with 295

int Find_ngb_simple(const int ipart,  const float hsml, int *ngblist);

static float global_density_model(const int ipart);
static inline float sph_kernel_M4(const float r, const float h);
static inline double sph_kernel_WC2(const float r, const float h);
static inline double sph_kernel_WC6(const float r, const float h);
static inline float gravity_kernel(const float r, const float h);

/* Settle SPH particle with weighted Voronoi tesselations (Diehl+ 2012).
 * Here hsml is not the SPH smoothing length, but is related to a local 
 * metric defined ultimately by the density model.
 * Relaxation is done in units of the boxsize, hence the box volume is 1 
 * We do not stop on the SPH sampling error, but on the distribution of particle
 * displacements relative to the local mean particle density. This is because 
 * the minimum energy state does not coincide with the lowest sampling error. */

void Regularise_sph_particles()
{
	const int maxiter = 128;
	const double mps_frac = 5; 		// move this fraction of the mean particle sep
	const double step_red = 0.95; 	// force convergence at this rate
	const double bin_limits[3] = { -1, 5, -1 }; // displacement limits in 100%, 10%, 1%

    const int nPart = Param.Npart[0];
    const double boxsize = Param.Boxsize;

	const double npart2percent = 100.0/nPart;

    printf("Starting iterative SPH regularisation \n"
			"   max %d iterations, mpsfrac=%g, force convergence at %g \n"
			"   bin limits: %g %g %g\n",
			maxiter, mps_frac, step_red, bin_limits[0], bin_limits[1], bin_limits[2]); 
	fflush(stdout);

    float *hsml = NULL;
    hsml = Malloc(nPart * sizeof(*hsml));
    
    float *displ[3] = { NULL }; 
   
    displ[0] = Malloc(nPart * sizeof(**displ));
    displ[1] = Malloc(nPart * sizeof(**displ));
    displ[2] = Malloc(nPart * sizeof(**displ));

	double rho_mean = nPart * Param.Mpart[0] / p3(boxsize);
	double step_mean = boxsize/pow(nPart, 1.0/3.0) / mps_frac;

	double errLast = DBL_MAX, errLastTree = DBL_MAX;
	double errDiff = DBL_MAX;

	double last_cnt = DBL_MAX;

    int it = -1;

    for (;;) {
			
		Sort_Particles_By_Peano_Key();	
	
		Build_Tree();	

		Find_sph_quantities();
	
        double vSphSum = 0; // total volume defined by hsml
 
		#pragma omp parallel for shared(hsml) reduction(+:vSphSum)
        for (int ipart = 0; ipart < nPart; ipart++) { // find hsml

            float rho = global_density_model(ipart);

			SphP[ipart].Rho_Model= rho;

            hsml[ipart] = pow(WVTNNGB * Param.Mpart[0]/rho/fourpithird, 1./3.);
            
            vSphSum += p3(hsml[ipart]);
        }

        float norm_hsml = pow(WVTNNGB/vSphSum/fourpithird , 1.0/3.0);
    
		#pragma omp parallel for
        for (int ipart = 0; ipart < nPart; ipart++) 
            hsml[ipart] *= norm_hsml;

		#pragma omp parallel for shared(displ, hsml, P) schedule(dynamic, nPart/Omp.NThreads/256)
        for (int ipart = 0; ipart < nPart; ipart++) { 

            displ[0][ipart] = displ[1][ipart] = displ[2][ipart] = 0;

            int ngblist[NGBMAX] = { 0 };

            //int ngbcnt = Find_ngb_simple(ipart, hsml[ipart]*boxsize, ngblist);
            int ngbcnt = Find_ngb_tree(ipart, hsml[ipart]*boxsize, ngblist);

			for (int i = 0; i < ngbcnt; i++) { // neighbour loop

				int jpart = ngblist[i];

                if (ipart == jpart)
                    continue;

                float dx = (P[ipart].Pos[0] - P[jpart].Pos[0])/boxsize;
	    		float dy = (P[ipart].Pos[1] - P[jpart].Pos[1])/boxsize;
    		    float dz = (P[ipart].Pos[2] - P[jpart].Pos[2])/boxsize;
			
                dx = dx > 0.5 ? dx-1 : dx; // find closest image
                dy = dy > 0.5 ? dy-1 : dy;
                dz = dz > 0.5 ? dz-1 : dz;

                dx = dx < -0.5 ? dx+1 : dx;
                dy = dy < -0.5 ? dy+1 : dy;
                dz = dz < -0.5 ? dz+1 : dz; 

                float r2 = (dx*dx + dy*dy + dz*dz);
                
                float h = 0.5 * (hsml[ipart] + hsml[jpart]);

    		    if (r2 > p2(h)) 
                    continue ;

    		    float r = sqrt(r2);

		        float wk = sph_kernel_WC6(r, h);
                
				/* scale mean step size with local density */

				double dens_contrast = pow(SphP[ipart].Rho_Model/rho_mean, 1/3);
				double step = step_mean / dens_contrast;
				
				displ[0][ipart] += step * hsml[ipart] * wk * dx/r;
                displ[1][ipart] += step * hsml[ipart] * wk * dy/r;
                displ[2][ipart] += step * hsml[ipart] * wk * dz/r;
            }
        }

        int cnt_100 = 0, cnt_10 = 0, cnt_1 = 0 ;

		#pragma omp parallel for reduction(+:cnt_100,cnt_10,cnt_1)
        for (int ipart = 0; ipart < nPart; ipart++) { // move particles

			float rho = global_density_model(ipart);

            float d = sqrt(p2(displ[0][ipart])
                    + p2( displ[1][ipart]) + p2( displ[2][ipart]));

            float d_mps = pow(Param.Mpart[0] / rho / DESNNGB, 1.0/3.0);

            if (d > 1 * d_mps) // simple distribution function of displs
                cnt_100++;
            if (d > 0.1 * d_mps)
                cnt_10++;
            if (d > 0.01 * d_mps)
                cnt_1++;

            P[ipart].Pos[0] += displ[0][ipart];
            P[ipart].Pos[1] += displ[1][ipart];
            P[ipart].Pos[2] += displ[2][ipart];

            while (P[ipart].Pos[0] < 0) // keep it in the box
                P[ipart].Pos[0] += boxsize;

            while (P[ipart].Pos[0] > boxsize)
                P[ipart].Pos[0] -= boxsize;
        
            while (P[ipart].Pos[1] < 0)
                P[ipart].Pos[1] += boxsize;

            while (P[ipart].Pos[1] > boxsize)
                P[ipart].Pos[1] -= boxsize;
            
            while (P[ipart].Pos[2] < 0)
                P[ipart].Pos[2] += boxsize;

            while (P[ipart].Pos[2] > boxsize)
                P[ipart].Pos[2] -= boxsize;
        }
	
		int nIn = 0;

        double  errMax = 0, errMean = 0; 

		#pragma omp parallel for reduction(+:errMean,nIn) reduction(max:errMax)
        for (int  ipart = 0; ipart < nPart; ipart++) { // get error

        	float rho = global_density_model(ipart);

            float err = fabs(SphP[ipart].Rho-rho) / rho;

            errMax = fmax(err, errMax);

            errMean += err;

			nIn++;
        }

        errMean /= nIn;

		errDiff = (errLast - errMean) / errMean;

		double bins[3] = { cnt_100*npart2percent, cnt_10*npart2percent, cnt_1*npart2percent };

		printf("   #%04d: Delta %4g%% > 1; %4g%% > 1/10; %4g%% > 1/100 of d_mps\n" 
			   "          Error max=%3g; mean=%03g; diff=%03g step_mean=%g\n",
				it, bins[0], bins[1], bins[2], errMax, errMean,errDiff, step_mean); 

		errLast = errMean;

		if (cnt_10 > last_cnt)  // force convergence if distribution doesnt tighten
            step_mean *= step_red;

		last_cnt = cnt_10;

		if ((bins[0] < bin_limits[0]) ||
			(bins[1] < bin_limits[1]) ||
			(bins[2] < bin_limits[2]) )
				break;

		if (it++ >= maxiter)
			break;
    }

    Free(hsml); Free(displ[0]); Free(displ[1]); Free(displ[2]);

    printf("\ndone\n\n"); fflush(stdout);
    
    return ;
}

static float global_density_model(const int ipart)
{
    const double boxhalf = Param.Boxsize*0.5;

	const double px = P[ipart].Pos[0], 
		  		 py = P[ipart].Pos[1], 
		  		 pz = P[ipart].Pos[2];

    double rho = 0;  

    for (int i = 0; i < Param.Nhalos; i++) {

		if (Halo[i].Mass[0] == 0) // DM only halos
			continue;

        double dx = px - Halo[i].D_CoM[0] - boxhalf;
        double dy = py - Halo[i].D_CoM[1] - boxhalf;
        double dz = pz - Halo[i].D_CoM[2] - boxhalf;

        double r2 = dx*dx + dy*dy + dz*dz;

		double rho_i = Gas_Density_Profile(sqrt(r2), i);

		rho = fmax(rho_i, rho);
    }

    return rho;
}
    
static inline double sph_kernel_WC2(const float r, const float h)
{   
	const float u= r/h;
    const float t = 1-u;

    return 21/(2*pi)*t*t*t*t*(1+4*u);
}

static inline float gravity_kernel(const float r, const float h)
{
    const float epsilon = 0.1;
    const float offset = h / (h + epsilon);
    const float val = h / (r + epsilon) - offset;

    return val * val;
}

static inline double sph_kernel_WC6(const float r, const float h)
{   
	const double u = r/h;
    const double t = 1-u;

    return 1365.0/(64*pi) *t*t*t*t*t*t*t*t*(1+8*u + 25*u*u + 32*u*u*u);
}

static inline float sph_kernel_M4(const float r, const float h) // cubic spline
{
	double wk = 0;
   	double u = r/h;
 
	if(u < 0.5) 
		wk = (2.546479089470 + 15.278874536822 * (u - 1) * u * u);
	else
		wk = 5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u);
	
	return wk/p3(h);
}

int Find_ngb_simple(const int ipart,  const float hsml, int *ngblist)
{
    const float boxhalf = Param.Boxsize/2;
    const float boxsize = Param.Boxsize;

    int ngbcnt = 0;

    for (int i = 0; i < NGBMAX; i++)
        ngblist[i] = 0;
    
	for (int jpart=0; jpart < Param.Npart[0]; jpart++) {
        
        float dx = fabs(P[ipart].Pos[0] - P[jpart].Pos[0]);
        float dy = fabs(P[ipart].Pos[1] - P[jpart].Pos[1]);
        float dz = fabs(P[ipart].Pos[2] - P[jpart].Pos[2]);

        if (dx > boxhalf)	// find closest image 
			dx -= boxsize;

    	if (dy > boxhalf)
	    	dy -= boxsize;

		if (dz > boxhalf)
			dz -= boxsize;
		
    	float r2 = dx*dx + dy*dy + dz*dz;

        if (r2 < hsml*hsml) 
			ngblist[ngbcnt++] = jpart;

		if (ngbcnt == NGBMAX)
			break;
    }

    return ngbcnt ;
}

