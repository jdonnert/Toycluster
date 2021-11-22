
/* profiles */

void Setup_Profiles(const int);
void Setup_DM_Mass_Profile(const int j);
void Setup_Gas_Mass_Profile(const int j);

double DM_Density_Profile(const int, const float);
double Potential_Profile(const int, const double);
double DM_Potential_Profile(const int, const float);
double DM_Mass_Profile(const double, const int); 
double DM_Mass_Profile_NFW(const double, const int); 
double DM_Mass_Profile_HQ(const double r, const int i);
double Inverted_DM_Mass_Profile(const double, const int); 
double Hernquist_Density_Profile(const int, const double);
double DM_Potential_Profile_HQ(const int i, const double r);
double DM_Potential_Profile_NFW(const int i, const double r);
double Distribution_Function(const double E);

double Gas_Density_Profile(const double, const int); 
double Gas_Mass_Profile(const double, const int); 
double Inverted_Gas_Mass_Profile(double); 
double Gas_Potential_Profile(const int i, const double r);
double Mass_Profile_23(const double, const int); 
double Gas_Potential_Profile_23(const int i, const float r);

double Internal_Energy_Profile(const int,const double);
double Internal_Energy_Profile_Analytic(const int i, const double d);

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
void Bfld_from_turb_spectrum();
void Shift_Origin();
void Regularise_sph_particles();
void Show_mass_in_r200();
void Setup_Substructure();
void Reassign_particles_to_halos();
void Smooth_SPH_quantities();


int Halo_containing(const int, const float,const float,const float);
double Concentration_parameter(const int);
double Gas_core_radius(const int, char *);


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
