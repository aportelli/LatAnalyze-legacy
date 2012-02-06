/* latan_blas.h, part of LatAnalyze library
 *
 * Copyright (C) 2010, 2011, 2012 Antonin Portelli
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

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
