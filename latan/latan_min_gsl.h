#ifndef LATAN_MIN_GSL_H_
#define LATAN_MIN_GSL_H_

#include <latan/latan_globals.h>
#include <latan/latan_minimizer.h>

__BEGIN_DECLS

latan_errno minimize_gsl(mat x, double* f_min, min_func* f, void* param);

__END_DECLS

#endif