#ifndef LATAN_MASS_H_
#define LATAN_MASS_H_

#include <latan/latan_globals.h>
#include <latan/latan_fit.h>

__BEGIN_DECLS

typedef struct
{
    size_t start;
    size_t end;
    double mean;
    double sig;
} plat;

latan_errno effmass(mat *res, mat *mprop, const int parity);
latan_errno effmass_PCAC(mat *res, mat *mprop_AP, mat *mprop_PP);
plat *search_plat(size_t *nplat, mat *data, mat *sigdata,\
                  const size_t ntmax, const double nsig, const double tol);
latan_errno fit_data_mass_fit_tune(fit_data *d, mat *fit_init, mat *prop,\
                                   mat *em, mat *sigem,            \
                                   const int parity);

__END_DECLS

#endif