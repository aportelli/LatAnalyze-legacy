#ifndef LATAN_BLAS_H_
#define	LATAN_BLAS_H_

#include <latan/latan_globals.h>

/* level 1 */
latan_errno latan_blas_ddot(const mat *x, const mat *y, double *res);

/* level 3 */
latan_errno latan_blas_dgemm(const char opA, const char opB, double alpha,   \
                             const mat *A, const mat *B, double beta, mat *C);
latan_errno latan_blas_dsymm(const char side, const char uploA, double alpha,\
                             const mat *A, const mat *B, double beta, mat *C);

#endif
