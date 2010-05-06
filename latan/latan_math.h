#ifndef LATAN_MATH_H_
#define LATAN_MATH_H_

#include <latan/latan_globals.h>

#define SQ(x) ((x)*(x))
#define DRATIO(a,b) (((double)(a))/((double)(b)))
#define C_PI 3.1415926535897932384626433832795028841970

__BEGIN_DECLS

unsigned int binomial(const unsigned int n, const unsigned int p);
latan_errno finite_diff(mat ddat, const mat dat);

__END_DECLS

#endif