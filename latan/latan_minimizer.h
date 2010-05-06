#ifndef LATAN_MINIMIZER_H_
#define LATAN_MINIMIZER_H_

#include <latan/latan_globals.h>

__BEGIN_DECLS

typedef double min_func(const mat var, void* param);

latan_errno minimize(mat var, min_func* f, void* param);

__END_DECLS

#endif