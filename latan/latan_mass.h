#ifndef LATAN_MASS_H_
#define LATAN_MASS_H_

#include <latan/latan_globals.h>

__BEGIN_DECLS

typedef struct
{
	size_t start;
	size_t end;
	double mean;
	double sig;
} plat;

latan_errno effmass(mat res, const mat mprop, const int parity);
latan_errno effmass_PCAC(mat res, const mat mprop_AP, const mat mprop_PP);
plat* search_plat(size_t *nplat, const mat data, const mat sigdata,\
				  const size_t ntmax, const double nsig, const double tol);
								   

__END_DECLS

#endif