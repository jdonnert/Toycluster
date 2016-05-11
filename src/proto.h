#include "macro.h"
#include "sort.h"
#include "peano.h"
#include "cosmo.h"

void Read_param_file(char *);
void Set_units();
void Set_cosmology();
void Setup();
void Make_positions();
void Make_IDs();
void Read_positions();
void Center_positions();
void Make_velocities();
void Make_temperatures();
void Make_magnetic_field();
void Find_sph_quantities();
void Apply_kinematics();
void Show_mass_in_r200();
void Wvt_relax();
void Shift_particles();
void Write_output();
void Bfld_from_rotA_SPH();
void Shift_Origin();
void Regularise_sph_particles();
void Show_mass_in_r200();
void Setup_Substructure();
void Reassign_particles_to_halos();
extern void Smooth_SPH_quantities();

double Internal_Energy_Profile(const int,const double);


double Concentration_parameter(const int);
double Gas_core_radius(const int, char *);
double Hernquist_density_profile(const double, const double, const double);
double Beta_density_profile(const double, const double, const double, 
		const double);
double Mass_profile(const double, const double, const double, 
		const double, const bool); 
double Gas_density_profile(const double, const double, const double, 
		const double, const bool); 
int Halo_containing(const int, const float,const float,const float);
float Global_density_model(const int);


/* Helper Monkeys */
void *Malloc_info(const char* func, const char* file, const int line, 
        size_t size);
void *Realloc_info(const char* func, const char* file, const int line, 
        void *ptr, size_t size);
void Free_info(const char* func, const char* file, const int line, void *ptr);
void Assert_Info(const char *func, const char *file, int line, int64_t expr, 
        const char *errmsg, ...);
double U2T(double U);
double T2U(double T);
double Density(float rho);

/* Cosmo */
double Redshift2Time(const double);

/* From system libs */
double erand48(unsigned short *);
