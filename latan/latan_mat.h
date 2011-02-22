/* latan_mat.h, part of LatAnalyze library
 *
 * Copyright (C) 2010, 2011 Antonin Portelli
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
enum
{
    CPU_ALLOCATED = 0x01,\
    GPU_ALLOCATED = 0x02,\
    CPU_LAST      = 0x04,\
    GPU_LAST      = 0x08,\
    SYNCED        = 0x10
};

typedef struct 
{
    gsl_matrix *data_cpu;
    double *data_gpu;
    int mem_flag;
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
mat *mat_create_from_mat(mat *n);
mat *mat_create_from_ar(const double *ar, const size_t init_nrow,\
                       const size_t init_ncol);
mat **mat_ar_create(const size_t nmat, const size_t init_nrow,\
                   const size_t init_ncol);
#define mat_ar_create_from_dim(nmat,n) mat_ar_create(nmat,nrow(n),ncol(n))
void mat_destroy(mat *m);
void mat_ar_destroy(mat **m, const size_t nmat);

/** CPU/GPU transfer **/
extern latan_errno (*mat_on_cpu)(mat *m);
extern latan_errno (*mat_on_gpu)(mat *m);

/** flag management **/
#define MAT_CPU_LAST(m)\
{\
    (m)->mem_flag -= (m)->mem_flag&(GPU_LAST|SYNCED);\
    (m)->mem_flag |= CPU_LAST;\
}
#define MAT_GPU_LAST(m)\
{\
    (m)->mem_flag -= (m)->mem_flag&(CPU_LAST|SYNCED);\
    (m)->mem_flag |= GPU_LAST;\
}
#define MAT_SYNCED(m)\
{\
    (m)->mem_flag -= (m)->mem_flag&(CPU_LAST|GPU_LAST);\
    (m)->mem_flag |= SYNCED;\
}

/** access **/
size_t nrow(mat *m);
size_t ncol(mat *m);
double mat_get(mat *m, const size_t i, const size_t j);
void mat_set(mat *m, const size_t i, const size_t j, const double val);
latan_errno mat_get_subm(mat *m, mat *n, const size_t k1, const size_t l1, \
                         const size_t k2, const size_t l2);
latan_errno mat_set_subm(mat *m, mat *n, const size_t k1, const size_t l1, \
                         const size_t k2, const size_t l2);
#define MAT_PT_SUBM(m,n,k1,l1,k2,l2)\
{\
    gsl_matrix_view _nview;\
    mat_on_cpu(n);\
    if ((m)->mem_flag & CPU_ALLOCATED)\
    {\
        LATAN_WARNING("memory leak : modifying a pointer on an allocated memory zone",\
                      LATAN_EFAULT);\
        (m)->mem_flag -= CPU_ALLOCATED;\
    }\
    _nview        = gsl_matrix_submatrix((n)->data_cpu,k1,l1,(k2)-(k1)+1,\
                                         (l2)-(l1)+1);                   \
    (m)->data_cpu = &(_nview.matrix);\
    \
    MAT_CPU_LAST(m);\
}
latan_errno mat_get_diag(mat *diag, mat *m);
latan_errno mat_set_diag(mat *m, mat *diag);
latan_errno mat_set_step(mat *m, const double x0, const double step);
latan_errno mat_set_from_ar(mat *m, const double *ar);
#define mat_inc(m,i,j,val) mat_set(m,i,j,mat_get(m,i,j)+val)
#define mat_pp(m,i,j) mat_inc(m,i,j,1.0)
double mat_get_min(mat *m);
double mat_get_max(mat *m);

/** tests **/
bool mat_is_samedim(mat *m, mat *n);
bool mat_is_square(mat *m);
bool mat_is_row_vector(mat *m);
bool mat_is_col_vector(mat *m);
bool mat_is_vector(mat *m);

/** operations **/
void mat_zero(mat *m);
void mat_cst(mat *m, const double x);
void mat_rand_u(mat *m, const double a, const double b);
void mat_id(mat *m);
latan_errno mat_cp(mat *m, mat *n);
latan_errno mat_eqadd(mat *m, mat *n);
latan_errno mat_add(mat *m, mat *n, mat *o);
#define mat_eqadds(m,s) mat_adds(m,m,s);
latan_errno mat_adds(mat *m, mat *n, const double s);
latan_errno mat_eqsub(mat *m, mat *n);
latan_errno mat_sub(mat *m, mat *n, mat *o);
#define mat_eqmul_l(m,n) mat_mul(m,n,m);
#define mat_eqmul_r(m,n) mat_mul(m,m,n);
#define mat_eqmul_l_t(m,n) mat_mul_tn(m,n,m);
#define mat_eqmul_r_t(m,n) mat_mul_nt(m,m,n);
latan_errno mat_mul_nn(mat *m, mat *n, mat *o);
latan_errno mat_mul_nt(mat *m, mat *n, mat *o);
latan_errno mat_mul_tn(mat *m, mat *n, mat *o);
latan_errno mat_mul_tt(mat *m, mat *n, mat *o);
#define mat_mul(m,n,o) mat_mul_nn(m,n,o);
latan_errno mat_dvmul(double *res, mat *m, mat *n);
latan_errno mat_eqtranspose(mat *m);
latan_errno mat_transpose(mat *m, mat *n);
latan_errno mat_eqmulp(mat *m, mat *n);
latan_errno mat_mulp(mat *m, mat *n, mat *o);
latan_errno mat_eqmuls(mat *m, const double s);
latan_errno mat_muls(mat *m, mat *n, const double s);
latan_errno mat_eqdivp(mat *m, mat *n);
latan_errno mat_divp(mat *m, mat *n, mat *o);
#define mat_eqabs(m) mat_abs(m,m)
latan_errno mat_abs(mat *m, mat *n);
#define mat_eqsqrt(m) mat_sqrt(m,m)
latan_errno mat_sqrt(mat *m, mat *n);

/** linear algebra **/
#define mat_eqinv(m) mat_inv(m,m);
latan_errno mat_inv(mat *m, mat *n);

__END_DECLS

#endif