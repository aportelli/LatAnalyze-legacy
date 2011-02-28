/* latan_mat.c, part of LatAnalyze library
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

#include <latan/latan_mat.h>
#include <latan/latan_includes.h>
#include <latan/latan_blas.h>
#include <latan/latan_math.h>
#include <latan/latan_rand.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>

/*                              allocation                                  */
/****************************************************************************/
mat *mat_create(const size_t init_nrow, const size_t init_ncol)
{
    gsl_error_handler_t *error_handler;
    mat *m;
    
    MALLOC_ERRVAL(m,mat *,1,NULL);
    if ((init_nrow > 0)&&(init_ncol > 0))
    {
        error_handler = gsl_set_error_handler_off();
        m->data_cpu   = gsl_matrix_alloc(init_nrow,init_ncol);
        if (m == NULL)
        {
            LATAN_ERROR_VAL("memory allocation failed",LATAN_ENOMEM,NULL);
        }
        gsl_set_error_handler(error_handler);
    }
    
    return m;
}

mat *mat_create_from_mat(const mat *n)
{
    mat *m;
    
    m = mat_create_from_dim(n);
    
    mat_cp(m,n);
    
    return m;
}

mat *mat_create_from_ar(const double *ar, const size_t init_nrow,\
                        const size_t init_ncol)
{
    mat *m;
    
    m = mat_create(init_nrow,init_ncol);
    
    mat_set_from_ar(m,ar);
    
    return m;
}

mat **mat_ar_create(const size_t nmat, const size_t init_nrow,\
                    const size_t init_ncol)
{
    size_t i;
    mat **m;
    
    MALLOC_NOERRET(m,mat**,nmat);
    for (i=0;i<nmat;i++)
    {
        m[i] = mat_create(init_nrow,init_ncol);
    }
    
    return m;
}

void mat_destroy(mat *m)
{
    gsl_matrix_free(m->data_cpu);
    FREE(m);
}

void mat_ar_destroy(mat **m, const size_t nmat)
{
    size_t i;
    
    for (i=0;i<nmat;i++)
    {
        mat_destroy(m[i]);
    }
    FREE(m);
}

/*                              access                                      */
/****************************************************************************/
size_t nrow(const mat *m)
{
    return m->data_cpu->size1;
}

size_t ncol(const mat *m)
{
    return m->data_cpu->size2;
}

double mat_get(const mat *m, const size_t i, const size_t j)
{
    if ((i>=nrow(m))||(j>=ncol(m)))
    {
        LATAN_ERROR_NORET("index out of range",LATAN_EBADLEN);
    }

    return gsl_matrix_get(m->data_cpu,i,j);
}

void mat_set(mat *m, const size_t i, const size_t j, const double val)
{
    if ((i>=nrow(m))||(j>=ncol(m)))
    {
        LATAN_ERROR_VOID("index out of range",LATAN_EBADLEN);
    }
    
    gsl_matrix_set(m->data_cpu,i,j,val);
}

latan_errno mat_get_subm(mat *m, const mat *n, const size_t k1,           \
                         const size_t l1, const size_t k2, const size_t l2)
{
    latan_errno status;
    gsl_matrix_const_view nview = gsl_matrix_const_submatrix(n->data_cpu,      \
                                                             (size_t)(k1),     \
                                                             (size_t)(l1),     \
                                                             (size_t)(k2-k1+1),\
                                                             (size_t)(l2-l1+1));
    
    if ((k2-k1+1>nrow(n))||(k2-k1+1<1)||(l2-l1+1>ncol(n))||(l2-l1+1<1))
    {
        LATAN_ERROR("invalid sub-matrix dimensions",LATAN_EBADLEN);
    }
    if ((k2-k1+1!=nrow(m))||(l2-l1+1!=ncol(m)))
    {
        LATAN_ERROR("sub-matrix and destination matrix dimensions do not match"\
                    ,LATAN_EBADLEN);
    }
    
    status = gsl_matrix_memcpy(m->data_cpu,&(nview.matrix));
    
    return status;
}

latan_errno mat_set_subm(mat *m, const mat *n, const size_t k1,          \
                         const size_t l1,const size_t k2, const size_t l2)
{
    latan_errno status;
    gsl_matrix_view mview;
    
    if ((k2-k1+1>nrow(m))||(k2-k1+1<1)||(l2-l1+1>ncol(m))||(l2-l1+1<1))
    {
        LATAN_ERROR("invalid sub-matrix dimensions",LATAN_EBADLEN);
    }
    if ((k2-k1+1!=nrow(n))||(l2-l1+1!=ncol(n)))
    {
        LATAN_ERROR("destination sub-matrix and matrix dimensions do not match"\
                    ,LATAN_EBADLEN);
    }
    
    mview  = gsl_matrix_submatrix(m->data_cpu,(size_t)(k1),(size_t)(l1),\
                                  (size_t)(k2-k1+1),(size_t)(l2-l1+1));
    status = gsl_matrix_memcpy(&(mview.matrix),n->data_cpu);
    
    return status;
}

latan_errno mat_get_diag(mat *diag, const mat *m)
{
    size_t min_dim;
    size_t i;
    
    min_dim = MIN(nrow(m),ncol(m));
    if (min_dim != nrow(diag))
    {
        LATAN_ERROR("trying to set the diagonal of a matrix with invalid dimensions",\
                    LATAN_EBADLEN);
    }
    
    for (i=0;i<min_dim;i++)
    {
        mat_set(diag,i,0,mat_get(m,i,i));
    }
    
    return LATAN_SUCCESS;
}

latan_errno mat_set_diag(mat *m, const mat *diag)
{
    size_t min_dim;
    size_t i;
    
    min_dim = MIN(nrow(m),ncol(m));
    if (min_dim != nrow(diag))
    {
        LATAN_ERROR("trying to set the diagonal of a matrix with invalid dimensions",\
                    LATAN_EBADLEN);
    }
    
    for (i=0;i<min_dim;i++)
    {
        mat_set(m,i,i,mat_get(diag,i,0));
    }

    return LATAN_SUCCESS;
}

latan_errno mat_set_step(mat *m, const double x0, const double step)
{
    size_t i;
    
    if (ncol(m) > 1)
    {
        LATAN_ERROR("step matrix output is not a column vector",LATAN_EBADLEN);
    }
    
    for (i=0;i<nrow(m);i++)
    {
        mat_set(m,i,0,x0+((double)(i))*step);
    }
    
    return LATAN_SUCCESS;
}

latan_errno mat_set_from_ar(mat *m, const double *ar)
{
    gsl_matrix_const_view ar_view = gsl_matrix_const_view_array(ar,\
                                                                nrow(m),\
                                                                ncol(m));
    latan_errno status;
    
    status = gsl_matrix_memcpy(m->data_cpu,&(ar_view.matrix));
    
    return status;
}

double mat_get_min(const mat *m)
{
    return gsl_matrix_min(m->data_cpu);
}

double mat_get_max(const mat *m)
{
    return gsl_matrix_max(m->data_cpu);
}

/*                              tests                                       */
/****************************************************************************/
bool mat_is_samedim(const mat *m, const mat *n)
{
    return ((nrow(m) == nrow(n))&&(ncol(m) == ncol(n)));
}

bool mat_is_square(const mat *m)
{
    return (nrow(m) == ncol(m));
}

bool mat_is_row_vector(const mat *m)
{
    return (nrow(m) == 1);
}

bool mat_is_col_vector(const mat *m)
{
    return (ncol(m) == 1);
}

bool mat_is_vector(const mat *m)
{
    return ((nrow(m) == 1)||(ncol(m) == 1));
}

/*                              operations                                  */
/****************************************************************************/
void mat_zero(mat *m)
{
    gsl_matrix_set_zero(m->data_cpu);
}

void mat_cst(mat *m, const double x)
{
    gsl_matrix_set_all(m->data_cpu,x);
}

void mat_rand_u(mat *m, const double a, const double b)
{
    size_t i,j;

    FOR_VAL(m,i,j)
    {
        mat_set(m,i,j,rand_u(a,b));
    }
}

void mat_id(mat *m)
{
    gsl_matrix_set_identity(m->data_cpu);
}

latan_errno mat_cp(mat *m, const mat *n)
{
    latan_errno status;
    
    status = LATAN_SUCCESS;
    
    if (!mat_is_samedim(m,n))
    {
        LATAN_ERROR("matrix copy with dimension mismatch",LATAN_EBADLEN);
    }
    
    LATAN_UPDATE_STATUS(status,gsl_matrix_memcpy(m->data_cpu,n->data_cpu));
   
    return status;
}

latan_errno mat_eqadd(mat *m, const mat *n)
{
    latan_errno status;
    size_t i;
    
    status = LATAN_SUCCESS;
    
    if (!mat_is_samedim(m,n))
    {
        LATAN_ERROR("operation between matrices with dimension mismatch",\
                    LATAN_EBADLEN);
    }

    if ((m->data_cpu->tda == ncol(m))&&(n->data_cpu->tda == ncol(n)))
    {
        cblas_daxpy((int)(ncol(m)*nrow(m)),1.0,n->data_cpu->data,1,\
                    m->data_cpu->data,1);
    }
    else
    {
        for (i=0;i<nrow(m);i++)
        {
            cblas_daxpy((int)(ncol(m)),1.0,                    \
                        n->data_cpu->data+i*n->data_cpu->tda,1,\
                        m->data_cpu->data+i*m->data_cpu->tda,1);
        }
    }
    
    return status;
}

latan_errno mat_add(mat *m, const mat *n, const mat *o)
{
    latan_errno status;
    
    status = LATAN_SUCCESS;
    
    LATAN_UPDATE_STATUS(status,mat_cp(m,n));
    LATAN_UPDATE_STATUS(status,mat_eqadd(m,o));
    
    return status;
}

latan_errno mat_adds(mat *m, const mat *n, const double s)
{
    latan_errno status;
    mat *cst;

    cst = mat_create_from_dim(n);

    mat_cst(cst,s);
    status = mat_add(m,n,cst);

    mat_destroy(cst);

    return status;
}

latan_errno mat_eqsub(mat *m, const mat *n)
{
    latan_errno status;
    size_t i;

    status = LATAN_SUCCESS;

    if (!mat_is_samedim(m,n))
    {
        LATAN_ERROR("operation between matrices with dimension mismatch",\
                    LATAN_EBADLEN);
    }

    if ((m->data_cpu->tda == ncol(m))&&(n->data_cpu->tda == ncol(n)))
    {
        cblas_daxpy((int)(ncol(m)*nrow(m)),-1.0,n->data_cpu->data,1,\
                    m->data_cpu->data,1);
    }
    else
    {
        for (i=0;i<nrow(m);i++)
        {
            cblas_daxpy((int)(ncol(m)),-1.0,                    \
                        n->data_cpu->data+i*n->data_cpu->tda,1, \
                        m->data_cpu->data+i*m->data_cpu->tda,1);
        }
    }

    return status;
}

latan_errno mat_sub(mat *m, const mat *n, const mat *o)
{
    latan_errno status;
    
    status = LATAN_SUCCESS;
    
    LATAN_UPDATE_STATUS(status,mat_cp(m,n));
    LATAN_UPDATE_STATUS(status,mat_eqsub(m,o));
    
    return status;
}

latan_errno mat_mul(mat *m, const mat *n, const char opn, const mat *o,\
                    const char opo)
{
    latan_errno status;
    bool have_duplicate_arg;
    mat *pt,*buf;
    double dbuf;
    
    status             = LATAN_SUCCESS;
    have_duplicate_arg = ((m == n)||(m == o));
    buf                = NULL;
    
    if (have_duplicate_arg)
    {
        buf = mat_create_from_dim(m);
        pt  = buf;
    }
    else
    {
        pt = m;
    }

    if ((nrow(m) == 1)&&(ncol(m) == 1))
    {
        LATAN_UPDATE_STATUS(status,latan_blas_ddot(n,o,&dbuf));
        mat_set(m,0,0,dbuf);
    }
    else
    {
        LATAN_UPDATE_STATUS(status,latan_blas_dgemm(opn,opo,1.0,n,o,0.0,pt));
    }
    
    if (have_duplicate_arg)
    {
        LATAN_UPDATE_STATUS(status,mat_cp(m,buf));
        mat_destroy(buf);
    }
    
    return status;
}

latan_errno mat_eqtranspose(mat *m)
{
    latan_errno status;
    
    if (!mat_is_square(m))
    {
        LATAN_ERROR("cannot auto-transpose a non-square matrix",LATAN_ENOTSQR);
    }
    
    status = gsl_matrix_transpose(m->data_cpu);
    
    return status;
}

latan_errno mat_transpose(mat *m, const mat *n)
{
    latan_errno status;
    
    if ((nrow(m) != ncol(n))||(ncol(m) != nrow(n)))
    {
        LATAN_ERROR("operation between matrices with dimension mismatch",\
                    LATAN_EBADLEN);
    }
    
    status = gsl_matrix_transpose_memcpy(m->data_cpu,n->data_cpu);
    
    return status;
}

latan_errno mat_inv(mat *m, const mat *n)
{
    latan_errno status;
    int signum;
    mat *LU;
    size_t i;
    gsl_permutation *perm;
    gsl_error_handler_t *error_handler;
    
    status = LATAN_SUCCESS;
    
    if (!mat_is_square(m))
    {
        LATAN_ERROR("cannot invert a non-square matrix",LATAN_ENOTSQR);
    }
    
    LU = mat_create_from_mat(n);
    error_handler = gsl_set_error_handler(&latan_error);
    perm = gsl_permutation_alloc(nrow(n));
    gsl_set_error_handler(error_handler);
    
    LATAN_UPDATE_STATUS(status,gsl_linalg_LU_decomp(LU->data_cpu,perm,&signum));
    for (i=0;i<nrow(LU);i++)
    {
        if (mat_get(LU,i,i) == 0.0)
        {
            LATAN_ERROR("trying to invert a singular matrix",LATAN_EDOM);
        }
    }
    LATAN_UPDATE_STATUS(status,gsl_linalg_LU_invert(LU->data_cpu,perm,\
                                                    m->data_cpu));
    
    mat_destroy(LU);
    gsl_permutation_free(perm);
    
    return status;
}

latan_errno mat_eqmulp(mat *m, const mat *n)
{
    latan_errno status;
    
    status = LATAN_SUCCESS;
    
    if (!(mat_is_samedim(m,n)))
    {
        LATAN_ERROR("operation between matrices with dimension mismatch",\
                    LATAN_EBADLEN);
    }
    
    LATAN_UPDATE_STATUS(status,gsl_matrix_mul_elements(m->data_cpu,n->data_cpu));
    
    return status;
}

latan_errno mat_mulp(mat *m, const mat *n, const mat *o)
{
    latan_errno status;
    
    status = LATAN_SUCCESS;
    
    LATAN_UPDATE_STATUS(status,mat_cp(m,n));
    LATAN_UPDATE_STATUS(status,mat_eqmulp(m,o));
    
    return status;
}

latan_errno mat_eqmuls(mat *m, const double s)
{
    latan_errno status;
    
    status = LATAN_SUCCESS;
    
    LATAN_UPDATE_STATUS(status,gsl_matrix_scale(m->data_cpu,s));
    
    return status;
}

latan_errno mat_muls(mat *m, const mat *n, const double s)
{
    latan_errno status;
    
    status = LATAN_SUCCESS;
    
    LATAN_UPDATE_STATUS(status,mat_cp(m,n));
    LATAN_UPDATE_STATUS(status,mat_eqmuls(m,s));
    
    return status;
}

latan_errno mat_eqdivp(mat *m, const mat *n)
{
    latan_errno status;
    
    status = LATAN_SUCCESS;
    
    if (!(mat_is_samedim(m,n)))
    {
        LATAN_ERROR("operation between matrices with dimension mismatch",\
                    LATAN_EBADLEN);
    }
    
    LATAN_UPDATE_STATUS(status,gsl_matrix_div_elements(m->data_cpu,n->data_cpu));
    
    return status;
}

latan_errno mat_divp(mat *m, const mat *n, const mat *o)
{
    latan_errno status;
    
    status = LATAN_SUCCESS;

    LATAN_UPDATE_STATUS(status,mat_cp(m,n));
    LATAN_UPDATE_STATUS(status,mat_eqdivp(m,o));

    return status;
}

latan_errno mat_abs(mat *m, const mat *n)
{
    size_t i,j;
    
    if (!mat_is_samedim(m,n))
    {
        LATAN_ERROR("operation between matrices with dimension mismatch",\
                    LATAN_EBADLEN);
    }
     
    FOR_VAL(m,i,j)
    {
        mat_set(m,i,j,fabs(mat_get(n,i,j)));
    }

    return LATAN_SUCCESS;
}

latan_errno mat_sqrt(mat *m, const mat *n)
{
    size_t i,j;
    
    if (!mat_is_samedim(m,n))
    {
        LATAN_ERROR("operation between matrices with dimension mismatch",\
                    LATAN_EBADLEN);
    }
      
    FOR_VAL(m,i,j)
    {
        mat_set(m,i,j,sqrt(fabs(mat_get(n,i,j))));
    }

    return LATAN_SUCCESS;
}
