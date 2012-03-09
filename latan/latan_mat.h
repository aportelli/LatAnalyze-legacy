/* latan_mat.h, part of LatAnalyze library
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

#ifndef LATAN_MAT_H_
#define LATAN_MAT_H_

#include <latan/latan_globals.h>
#include <gsl/gsl_matrix.h>

__BEGIN_DECLS

/* definition of the matrix type */
typedef enum
{
    MAT_GEN = 0,      /* general matrix    */
    MAT_SYM = 1 << 0, /* symmetric matrix  */
    MAT_POS = 1 << 1  /* positive spectrum */
} mat_flag;

typedef struct 
{
    gsl_matrix *data_cpu;
    unsigned int prop_flag;
} mat;

/* loop */
#define FOR_VAL(m,i,j)\
for (i=0;i<nrow(m);i++)\
for (j=0;j<ncol(m);j++)

/* functions */
/** allocation **/
mat *mat_create(const size_t init_nrow, const size_t init_ncol);
#define mat_create_from_dim(n) mat_create(nrow(n),ncol(n))
#define mat_create_from_trdim(n) mat_create(ncol(n),nrow(n))
mat *mat_create_from_mat(const mat *n);
mat *mat_create_from_ar(const double *ar, const size_t init_nrow,\
                        const size_t init_ncol);
mat **mat_ar_create(const size_t nmat, const size_t init_nrow,\
                    const size_t init_ncol);
#define mat_ar_create_from_dim(nmat,n) mat_ar_create(nmat,nrow(n),ncol(n))
void mat_destroy(mat *m);
void mat_ar_destroy(mat **m, const size_t nmat);

/** access **/
size_t nrow(const mat *m);
size_t ncol(const mat *m);
size_t nel(const mat *m);
double mat_get(const mat *m, const size_t i, const size_t j);
double mat_get_rm(const mat* m, const size_t ind);
void mat_set(mat *m, const size_t i, const size_t j, const double val);
void mat_set_rm(mat *m, const size_t ind, const double val);
latan_errno mat_get_subm(mat *m, const mat *n, const size_t k1,           \
                         const size_t l1, const size_t k2, const size_t l2);
latan_errno mat_set_subm(mat *m, const mat *n, const size_t k1,           \
                         const size_t l1, const size_t k2, const size_t l2);
latan_errno mat_get_diag(mat *diag, const mat *m);
latan_errno mat_set_diag(mat *m, const mat *diag);
latan_errno mat_set_step(mat *m, const double x0, const double step);
latan_errno mat_set_from_ar(mat *m, const double *ar);
#define mat_inc(m,i,j,val) mat_set(m,i,j,mat_get(m,i,j)+val)
#define mat_pp(m,i,j) mat_inc(m,i,j,1.0)
double mat_get_min(const mat *m);
double mat_get_max(const mat *m);
void mat_reset_assump(mat *m);
void mat_assume(mat *m, const mat_flag flag);

/** tests **/
bool mat_is_samedim(const mat *m, const mat *n);
bool mat_is_square(const mat *m);
bool mat_is_row_vector(const mat *m);
bool mat_is_col_vector(const mat *m);
bool mat_is_vector(const mat *m);
bool mat_is_assumed(const mat *m, const mat_flag flag);

/** sort **/
void mat_get_sind(size_t *sind, const mat *m);

/** operations **/
void mat_zero(mat *m);
latan_errno mat_cst(mat *m, const double x);
void mat_rand_u(mat *m, const double a, const double b);
void mat_id(mat *m);
latan_errno mat_cp(mat *m, const mat *n);
latan_errno mat_eqadd(mat *m, const mat *n);
latan_errno mat_add(mat *m, const mat *n, const mat *o);
#define mat_eqadds(m,s) mat_adds(m,m,s)
latan_errno mat_adds(mat *m, const mat *n, const double s);
latan_errno mat_eqsub(mat *m, const mat *n);
latan_errno mat_sub(mat *m, const mat *n, const mat *o);
latan_errno mat_mul(mat *m, const mat *n, const char opn, const mat *o,\
                    const char opo);
latan_errno mat_eqtranspose(mat *m);
latan_errno mat_transpose(mat *m, const mat *n);
latan_errno mat_eqmulp(mat *m, const mat *n);
latan_errno mat_mulp(mat *m, const mat *n, const mat *o);
#define mat_eqinvp(m) mat_invp(m,m)
latan_errno mat_invp(mat *m, const mat *n);
latan_errno mat_eqmuls(mat *m, const double s);
latan_errno mat_muls(mat *m, const mat *n, const double s);
latan_errno mat_eqdivp(mat *m, const mat *n);
latan_errno mat_divp(mat *m, const mat *n, const mat *o);
#define mat_eqabs(m) mat_abs(m,m)
latan_errno mat_abs(mat *m, const mat *n);
#define mat_eqsqrt(m) mat_sqrt(m,m)
latan_errno mat_sqrt(mat *m, const mat *n);
#define mat_eqexpp(m) mat_expp(m,m)
latan_errno mat_expp(mat *m, const mat *n);

/** linear algebra **/
#define mat_eqinv_LU(m) mat_inv_LU(m,m)
latan_errno mat_inv_LU(mat *m, const mat *n);
#define mat_eqinv_symChol(m) mat_inv_symChol(m,m)
latan_errno mat_inv_symChol(mat *m, const mat *n);
#define mat_eqpseudoinv(m) mat_pseudoinv(m,m);
latan_errno mat_pseudoinv(mat *m, const mat *n);

__END_DECLS

#endif
