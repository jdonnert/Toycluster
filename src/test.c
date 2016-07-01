


double Mass_profile(const double r, const double rho0, const double rc,
        const double rcut, const double beta, const bool Is_Cuspy)
{
    const double r2 = p2(r);
    const double rc2 = p2(rc);
    const double rcut2 = p2(rcut);

#ifdef GIVEPARAMS
    printf("r=%g\trho0=%g\trc=%g\trcut=%g\tbeta=%.5f",
            r, rho0, rc, rcut, beta);

    //  Probably a mighty bad idea ...
    gsl_error_handler_t * old_handler = gsl_set_error_handler_off();
    double Mr = rho0 * p3(r)/3 * gsl_sf_hyperg_2F1(1.5, 1.5*beta, 2.5, -p2(r/rc));
    gsl_set_error_handler(old_handler);
    printf("\tMr=%g\n", 4*pi*Mr);

    return 4 * pi * Mr;
