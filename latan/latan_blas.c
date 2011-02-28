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
        case 'h':
        case 'H':
            *op_no = CblasTrans;
            break;
        default:
            LATAN_ERROR("wrong matrix transposition flag",LATAN_EINVAL);
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

/*                              level 3                                     */
/****************************************************************************/
latan_errno latan_blas_dgemm(const char opA, const char opB, double alpha,\
                             const mat *A, const mat *B, double beta, mat *C)
{
    latan_errno status;
    CBLAS_TRANSPOSE_t opA_no,opB_no;
    size_t nropA,ncopA,nropB,ncopB;

    status = LATAN_SUCCESS;

    LATAN_UPDATE_STATUS(status,parse_op(&opA_no,opA));
    LATAN_UPDATE_STATUS(status,parse_op(&opB_no,opB));
    nropA = (opA_no == CblasNoTrans) ? nrow(A) : ncol(A);
    ncopA = (opA_no == CblasNoTrans) ? ncol(A) : nrow(A);
    nropB = (opB_no == CblasNoTrans) ? nrow(B) : ncol(B);
    ncopB = (opB_no == CblasNoTrans) ? ncol(B) : nrow(B);
    if ((nrow(C) != nropA)||(ncol(C) != ncopB)||(ncopA != nropB))
    {
        LATAN_ERROR("operation between matrices with dimension mismatch",\
                    LATAN_EBADLEN);
    }
    LATAN_UPDATE_STATUS(status,gsl_blas_dgemm(opA_no,opB_no,alpha,A->data_cpu,\
                        B->data_cpu,beta,C->data_cpu));

    return status;
}

latan_errno latan_blas_dsymm(const char side, const char uploA, double alpha,\
                             const mat *A, const mat *B, double beta, mat *C)
{
    latan_errno status;

    status = LATAN_SUCCESS;

    return status;
}


