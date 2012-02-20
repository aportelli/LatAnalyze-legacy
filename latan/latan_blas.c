/* latan_blas.c, part of LatAnalyze library
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

#include <latan/latan_blas.h>
#include <latan/latan_includes.h>
#include <latan/latan_io.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*                             flag parsers                                 */
/****************************************************************************/
static latan_errno parse_op(CBLAS_TRANSPOSE_t *op_no, const char op)
{
    switch(op)
    {
        case 'n':
        case 'N':
            *op_no = CblasNoTrans;
            break;
        case 't':
        case 'T':
            *op_no = CblasTrans;
            break;
        case 'h':
        case 'H':
            *op_no = CblasConjTrans;
            break;
        default:
            LATAN_ERROR("wrong matrix transposition flag",LATAN_EINVAL);
            break;
    }

    return LATAN_SUCCESS;
}

static latan_errno parse_side(CBLAS_SIDE_t *side_no, const char side)
{
    switch(side)
    {
        case 'l':
        case 'L':
            *side_no = CblasLeft;
            break;
        case 'r':
        case 'R':
            *side_no = CblasRight;
            break;
        default:
            LATAN_ERROR("wrong matrix side flag",LATAN_EINVAL);
            break;
    }

    return LATAN_SUCCESS;
}

static latan_errno parse_uplo(CBLAS_UPLO_t *uplo_no, const char uplo)
{
    switch (uplo)
    {
        case 'u':
        case 'U':
            *uplo_no = CblasUpper;
            break;
        case 'l':
        case 'L':
            *uplo_no = CblasLower;
            break;
        default:
            LATAN_ERROR("wrong matrix upper/lower flag",LATAN_EINVAL);
            break;
    }

    return LATAN_SUCCESS;
}

/*                              level 1                                     */
/****************************************************************************/
latan_errno latan_blas_ddot(const mat *x, const mat *y, double *res)
{
    latan_errno status;
    gsl_vector_view x_vview,y_vview;
    
    if (mat_is_row_vector(x))
    {
        x_vview = gsl_matrix_row(x->data_cpu,0);
    }
    else if (mat_is_col_vector(x))
    {
        x_vview = gsl_matrix_column(x->data_cpu,0);
    }
    else
    {
        LATAN_ERROR("vector operation on matrix",LATAN_EBADLEN);
    }
    if (mat_is_row_vector(y))
    {
        y_vview = gsl_matrix_row(y->data_cpu,0);
    }
    else if (mat_is_col_vector(y))
    {
        y_vview = gsl_matrix_column(y->data_cpu,0);
    }
    else
    {
        LATAN_ERROR("vector operation on matrix",LATAN_EBADLEN);
    }
    if (x_vview.vector.size != y_vview.vector.size)
    {
        LATAN_ERROR("operation between vectors with dimension mismatch",\
                    LATAN_EBADLEN);
    }
    status = gsl_blas_ddot(&(x_vview.vector),&(y_vview.vector),res);

    return status;
}

double latan_blas_dnrm2(const mat *x)
{
    gsl_vector_view x_vview;
    
    if (mat_is_row_vector(x))
    {
        x_vview = gsl_matrix_row(x->data_cpu,0);
    }
    else if (mat_is_col_vector(x))
    {
        x_vview = gsl_matrix_column(x->data_cpu,0);
    }
    else
    {
        LATAN_ERROR_VAL("vector operation on matrix",LATAN_EBADLEN,GSL_NAN);
    }
    
    return gsl_blas_dnrm2(&(x_vview.vector));
}

/*                              level 2                                     */
/****************************************************************************/
latan_errno latan_blas_dgemv(const char opA, const double alpha,            \
                             const mat *A, const mat *x, const double beta, \
                             mat *y)
{
    latan_errno status;
    gsl_vector_view x_vview,y_vview;
    CBLAS_TRANSPOSE_t opA_no;

    status = LATAN_SUCCESS;
    opA_no = CblasNoTrans;

    if (mat_is_row_vector(x))
    {
        x_vview = gsl_matrix_row(x->data_cpu,0);
    }
    else if (mat_is_col_vector(x))
    {
        x_vview = gsl_matrix_column(x->data_cpu,0);
    }
    else
    {
        LATAN_ERROR("vector operation on matrix",LATAN_EBADLEN);
    }
    if (mat_is_row_vector(y))
    {
        y_vview = gsl_matrix_row(y->data_cpu,0);
    }
    else if (mat_is_col_vector(y))
    {
        y_vview = gsl_matrix_column(y->data_cpu,0);
    }
    else
    {
        LATAN_ERROR("vector operation on matrix",LATAN_EBADLEN);
    }
    if (x_vview.vector.size != y_vview.vector.size)
    {
        LATAN_ERROR("operation between vectors with dimension mismatch",\
                    LATAN_EBADLEN);
    }
    if (ncol(A) != x_vview.vector.size)
    {
        LATAN_ERROR("operation between matrix and vector with dimension mismatch",\
                    LATAN_EBADLEN);
    }
    USTAT(parse_op(&opA_no,opA));
    USTAT(gsl_blas_dgemv(opA_no,alpha,A->data_cpu,&(x_vview.vector),beta,\
                         &(y_vview.vector)));

    return status;
}

latan_errno latan_blas_dsymv(const char uploA, const double alpha,          \
                             const mat *A, const mat *x, const double beta, \
                             mat *y)
{
    latan_errno status;
    gsl_vector_view x_vview,y_vview;
    CBLAS_UPLO_t uploA_no;

    status   = LATAN_SUCCESS;
    uploA_no = CblasUpper;
    
    if (mat_is_row_vector(x))
    {
        x_vview = gsl_matrix_row(x->data_cpu,0);
    }
    else if (mat_is_col_vector(x))
    {
        x_vview = gsl_matrix_column(x->data_cpu,0);
    }
    else
    {
        LATAN_ERROR("vector operation on matrix",LATAN_EBADLEN);
    }
    if (mat_is_row_vector(y))
    {
        y_vview = gsl_matrix_row(y->data_cpu,0);
    }
    else if (mat_is_col_vector(y))
    {
        y_vview = gsl_matrix_column(y->data_cpu,0);
    }
    else
    {
        LATAN_ERROR("vector operation on matrix",LATAN_EBADLEN);
    }
    if (x_vview.vector.size != y_vview.vector.size)
    {
        LATAN_ERROR("operation between vectors with dimension mismatch",\
                    LATAN_EBADLEN);
    }
    if (ncol(A) != x_vview.vector.size)
    {
        LATAN_ERROR("operation between matrix and vector with dimension mismatch",\
                    LATAN_EBADLEN);
    }
    if (!mat_is_square(A))
    {
        LATAN_ERROR("symmetric matrix is not square",LATAN_ENOTSQR);
    }
    if (!mat_is_assumed_sym(A))
    {
        LATAN_WARNING("matrix is not declared symmetric",LATAN_EINVAL);
    }
    USTAT(parse_uplo(&uploA_no,uploA));
    USTAT(gsl_blas_dsymv(uploA_no,alpha,A->data_cpu,&(x_vview.vector),beta,\
                         &(y_vview.vector)));

    return status;
}

/*                              level 3                                     */
/****************************************************************************/
latan_errno latan_blas_dgemm(const char opA, const char opB,                \
                             const double alpha, const mat *A, const mat *B,\
                             const double beta, mat *C)
{
    latan_errno status;
    CBLAS_TRANSPOSE_t opA_no,opB_no;
    size_t nropA,ncopA,nropB,ncopB;

    status = LATAN_SUCCESS;
    opA_no = CblasNoTrans;
    opB_no = CblasNoTrans;

    USTAT(parse_op(&opA_no,opA));
    USTAT(parse_op(&opB_no,opB));
    nropA = (opA_no == CblasNoTrans) ? nrow(A) : ncol(A);
    ncopA = (opA_no == CblasNoTrans) ? ncol(A) : nrow(A);
    nropB = (opB_no == CblasNoTrans) ? nrow(B) : ncol(B);
    ncopB = (opB_no == CblasNoTrans) ? ncol(B) : nrow(B);
    if ((nrow(C) != nropA)||(ncol(C) != ncopB)||(ncopA != nropB))
    {
        LATAN_ERROR("operation between matrices with dimension mismatch",\
                    LATAN_EBADLEN);
    }
    USTAT(gsl_blas_dgemm(opA_no,opB_no,alpha,A->data_cpu,B->data_cpu,beta,\
                         C->data_cpu));

    return status;
}

latan_errno latan_blas_dsymm(const char side, const char uploA,             \
                             const double alpha, const mat *A, const mat *B,\
                             const double beta, mat *C)
{
    latan_errno status;
    CBLAS_SIDE_t side_no;
    CBLAS_UPLO_t uploA_no;
    const mat *L,*R;

    status   = LATAN_SUCCESS;
    side_no  = CblasLeft;
    uploA_no = CblasUpper;

    USTAT(parse_side(&side_no,side));
    USTAT(parse_uplo(&uploA_no,uploA));
    L = (side_no == CblasLeft) ? A : B;
    R = (side_no == CblasLeft) ? B : A;
    if ((nrow(C) != nrow(L))||(ncol(C) != ncol(R))||(ncol(L) != nrow(R)))
    {
        LATAN_ERROR("operation between matrices with dimension mismatch",\
                    LATAN_EBADLEN);
    }
    if (!mat_is_square(A))
    {
        LATAN_ERROR("symmetric matrix is not square",LATAN_ENOTSQR);
    }
    if (!mat_is_assumed_sym(A))
    {
        LATAN_WARNING("matrix is not declared symmetric",LATAN_EINVAL);
    }
    USTAT(gsl_blas_dsymm(side_no,uploA_no,alpha,A->data_cpu,B->data_cpu,beta,\
                         C->data_cpu));

    return status;
}


