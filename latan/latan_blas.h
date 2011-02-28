#ifndef LATAN_BLAS_H_
#define	LATAN_BLAS_H_

#include <latan/latan_globals.h>

/* level 1 */
latan_errno latan_blas_ddot(const mat *x, const mat *y, double *res);

/* level 2 */
latan_errno latan_blas_dgemv(const char opA, const double alpha,            \
                             const mat *A, const mat *x, const double beta, \
                             mat *y);
latan_errno latan_blas_dsymv(const char uploA, const double alpha,          \
                             const mat *A, const mat *x, const double beta, \
                             mat *y);

/* level 3 */
latan_errno latan_blas_dgemm(const char opA, const char opB,                \
                             const double alpha, const mat *A, const mat *B,\
                             const double beta, mat *C);
latan_errno latan_blas_dsymm(const char side, const char uploA,             \
                             const double alpha, const mat *A, const mat *B,\
                             const double beta, mat *C);

#endif
