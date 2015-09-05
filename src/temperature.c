#include "globals.h"

static double F1(const double, const double, const double), 
              F2(const double, const double);

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

#pragma omp parallel
        for (size_t ipart = 0; ipart < Halo[i].Npart[0]; ipart++) {

            float dx = Halo[i].Gas[ipart].Pos[0] - Halo[i].D_CoM[0] - boxhalf;
            float dy = Halo[i].Gas[ipart].Pos[1] - Halo[i].D_CoM[1] - boxhalf;
            float dz = Halo[i].Gas[ipart].Pos[2] - Halo[i].D_CoM[2] - boxhalf;

            double r = sqrt(dx*dx + dy*dy + dz*dz);

            double u = Internal_Energy_Profile(i, r);

			double u_main = 0;

			if (i > 0 ) {
	         
            	dx = Halo[i].Gas[ipart].Pos[0] - Halo[0].D_CoM[0] - boxhalf;
            	dy = Halo[i].Gas[ipart].Pos[1] - Halo[0].D_CoM[1] - boxhalf;
            	dz = Halo[i].Gas[ipart].Pos[2] - Halo[0].D_CoM[2] - boxhalf;
				
            	r = sqrt(dx*dx + dy*dy + dz*dz);

			//	if (r < Halo[0].R200)
				u_main = Internal_Energy_Profile(0, r);
			}
			
			//u = fmax(u, u_main);

			Halo[i].SphP[ipart].U = u;
		}
    }

    printf("done\n\n"); fflush(stdout);

    return;
}

/* 
* to avoid negative temperatures we define rmax*sqrt3 as outer radius
*/

double Internal_Energy_Profile(const int i, const double d)
{
    const double G = Grav/p3(Unit.Length)*Unit.Mass*p2(Unit.Time);
		
    const double rho0 = Halo[i].Rho0; 
    const double a = Halo[i].A_hernq;
    const double rc = Halo[i].Rcore;
    double rmax = Param.Boxsize; // "open" T boundary
    double Mdm = 1.10 * Halo[i].Mass[1];

	double u = G / ( (adiabatic_index-1) ) * ( 1 + p2(d/rc) ) *
                ( Mdm * (F1(rmax, rc, a) - F1(d, rc, a))
                + 4*pi*rho0*p3(rc) * (F2(rmax, rc) - F2(d, rc) ) );
	return u;
}

static double F1(const double r, const double rc, const double a)
{
    const double rc2 = rc*rc;
    const double a2 = a*a;

    double result = (a2-rc2)*atan(r/rc) - rc*(a2+rc2)/(a+r) 
                + a*rc * log( (a+r)*(a+r) / (rc2 + r*r) );

    result *= rc / p2(a2 + rc2);

    return result;
}

static double F2(const double r, const double rc)
{
    return p2(atan(r/rc)) / (2*rc) + atan(r/rc)/r;
}
