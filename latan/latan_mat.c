#include <latan/latan_mat.h>
#include <latan/latan_includes.h>
#include <latan/latan_math.h>
#include <latan/latan_rand.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

/*								allocation									*/
/****************************************************************************/
mat mat_create(const size_t init_nrow, const size_t init_ncol)
{
	gsl_error_handler_t *error_handler;
	mat m;
	
	error_handler = gsl_set_error_handler(&latan_error);
	m = gsl_matrix_alloc(init_nrow,init_ncol);
	gsl_set_error_handler(error_handler);
	
	return m;
}

mat mat_create_from_mat(const mat n)
{
	mat m;
	
	m = mat_create_from_dim(n);
	
	mat_cp(m,n);
	
	return m;
}

mat mat_create_from_ar(const double *ar, const size_t init_nrow,\
					   const size_t init_ncol)
{
	mat m;
	
	m = mat_create(init_nrow,init_ncol);
	
	mat_set_from_ar(m,ar);
	
	return m;
}

mat *mat_ar_create(const size_t nmat, const size_t init_nrow,\
				   const size_t init_ncol)
{
	size_t i;
	mat *m;
	
	MALLOC_NOERRET(m,mat*,nmat);
	for (i=0;i<nmat;i++)
	{
		m[i] = mat_create(init_nrow,init_ncol);
	}
	
	return m;
}

void mat_destroy(mat m)
{
	gsl_matrix_free(m);
}

void mat_ar_destroy(mat *m, const size_t nmat)
{
	size_t i;
	
	for (i=0;i<nmat;i++)
	{
		mat_destroy(m[i]);
	}
	FREE(m);
}

/*								access										*/
/****************************************************************************/
size_t nrow(const mat m)
{
	return m->size1;
}

size_t ncol(const mat m)
{
	return m->size2;
}

double mat_get(const mat m, const size_t i, const size_t j)
{
	if ((i>=nrow(m))||(j>=ncol(m)))
	{
		LATAN_ERROR_NORET("index out of range",LATAN_EBADLEN);
	}

	return gsl_matrix_get(m,i,j);
}

void mat_set(mat m, const size_t i, const size_t j, const double val)
{
	if ((i>=nrow(m))||(j>=ncol(m)))
	{
		LATAN_ERROR_VOID("index out of range",LATAN_EBADLEN);
	}

	gsl_matrix_set(m,i,j,val);
}

latan_errno mat_get_subm(mat m, const mat n, const size_t k1, const size_t l1, \
						 const size_t k2, const size_t l2)
{
	gsl_matrix_view nview;
	
	if ((k2-k1+1>nrow(n))||(k2-k1+1<1)||(l2-l1+1>ncol(n))||(l2-l1+1<1))
	{
		LATAN_ERROR("invalid sub-matrix dimensions",LATAN_EBADLEN);
	}
	if ((k2-k1+1!=nrow(m))||(l2-l1+1!=ncol(m)))
	{
		LATAN_ERROR("sub-matrix and destination matrix dimensions do not match"\
					,LATAN_EBADLEN);
	}
	
	nview =	gsl_matrix_submatrix(n,(size_t)(k1),(size_t)(l1),			\
								 (size_t)(k2-k1+1),(size_t)(l2-l1+1));
	mat_cp(m,&(nview.matrix));
	
	return LATAN_SUCCESS;
}

latan_errno mat_set_subm(mat m, const mat n, const size_t k1, const size_t l1, \
						 const size_t k2, const size_t l2)
{
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
	
	mview =	gsl_matrix_submatrix(m,(size_t)(k1),(size_t)(l1),			\
								 (size_t)(k2-k1+1),(size_t)(l2-l1+1));
	mat_cp(&(mview.matrix),n);
	
	return LATAN_SUCCESS;
}

latan_errno mat_get_diag(mat diag, const mat m)
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

latan_errno mat_set_diag(mat m, const mat diag)
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

latan_errno mat_set_step(mat m, const double x0, const double step)
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

latan_errno mat_set_from_ar(mat m, const double *ar)
{
	gsl_matrix_const_view ar_view = gsl_matrix_const_view_array(ar,\
																nrow(m),\
																ncol(m));
	latan_errno status;
	
	status = gsl_matrix_memcpy(m,&(ar_view.matrix));
	
	return status;
}

double mat_get_min(const mat m)
{
	return gsl_matrix_min(m);
}

double mat_get_max(const mat m)
{
	return gsl_matrix_max(m);
}

/*								tests   									*/
/****************************************************************************/
bool mat_is_samedim(const mat m, const mat n)
{
	return ((nrow(m) == nrow(n))&&(ncol(m) == ncol(n)));
}

bool mat_is_square(const mat m)
{
	return (nrow(m) == ncol(m));
}

bool mat_is_row_vector(mat m)
{
	return (ncol(m) == 1);
}

bool mat_is_col_vector(mat m)
{
	return (nrow(m) == 1);
}

bool mat_is_vector(mat m)
{
	return ((nrow(m) == 1)||(ncol(m) == 1));
}

/*								operations									*/
/****************************************************************************/
void mat_zero(mat m)
{
	gsl_matrix_set_zero(m);
}

void mat_cst(mat m, const double x)
{
	gsl_matrix_set_all(m,x);
}

void mat_rand_u(mat m, const double a, const double b)
{
	size_t i,j;
	
	FOR_VAL(m,i,j)
	{
		mat_set(m,i,j,rand_u(a,b));
	}
}

void mat_id(mat m)
{
	gsl_matrix_set_identity(m);
}

latan_errno mat_cp(mat m, const mat n)
{
	latan_errno status;
	
	status = LATAN_SUCCESS;
	
	if (!mat_is_samedim(m,n))
	{
		LATAN_ERROR("matrix copy with dimension mismatch",LATAN_EBADLEN);
	}
	
	LATAN_UPDATE_STATUS(status,gsl_matrix_memcpy(m,n));
	
	return status;
}

latan_errno mat_eqadd(mat m, const mat n)
{
	latan_errno status;
	
	status = LATAN_SUCCESS;
	
	if (!(mat_is_samedim(m,n)))
	{
		LATAN_ERROR("operation between matrices with dimension mismatch",\
					LATAN_EBADLEN);
	}
	
	LATAN_UPDATE_STATUS(status,gsl_matrix_add(m,n));
	
	return status;
}

latan_errno mat_add(mat m, const mat n, const mat o)
{
	latan_errno status;
	
	status = LATAN_SUCCESS;
	
	LATAN_UPDATE_STATUS(status,mat_cp(m,n));
	LATAN_UPDATE_STATUS(status,mat_eqadd(m,o));
	
	return status;
}


latan_errno mat_eqsub(mat m, const mat n)
{
	latan_errno status;
	
	status = LATAN_SUCCESS;
	
	if (!(mat_is_samedim(m,n)))
	{
		LATAN_ERROR("operation between matrices with dimension mismatch",\
					LATAN_EBADLEN);
	}
	
	LATAN_UPDATE_STATUS(status,gsl_matrix_sub(m,n));
	
	return status;
}

latan_errno mat_sub(mat m, const mat n, const mat o)
{
	latan_errno status;
	
	status = LATAN_SUCCESS;
	
	LATAN_UPDATE_STATUS(status,mat_cp(m,n));
	LATAN_UPDATE_STATUS(status,mat_eqsub(m,o));
	
	return status;
}

latan_errno mat_mul_nn(mat m, const mat n, const mat o)
{
	latan_errno status;
	bool have_duplicate_arg;
	mat pt,buf;
	
	status             = LATAN_SUCCESS;
	have_duplicate_arg = ((m == n)||(m == o));
	
	if ((nrow(m) != nrow(n))||(ncol(m) != ncol(o))||(ncol(n) != nrow(o)))
	{
		LATAN_ERROR("operation between matrices with dimension mismatch",\
					LATAN_EBADLEN);
	}
	
	if (have_duplicate_arg)
	{
		buf = mat_create_from_dim(m);
		pt  = buf;
	}
	else
	{
		pt = m;
	}
	LATAN_UPDATE_STATUS(status,gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,\
											  1.0,n,o,0.0,pt));
	if (have_duplicate_arg)
	{
		LATAN_UPDATE_STATUS(status,mat_cp(m,buf));
		mat_destroy(buf);
	}
	
	return status;
}

latan_errno mat_mul_nt(mat m, const mat n, const mat o)
{
	latan_errno status;
	bool have_duplicate_arg;
	mat pt,buf;
	
	status             = LATAN_SUCCESS;
	have_duplicate_arg = ((m == n)||(m == o));
	
	if ((nrow(m) != nrow(n))||(ncol(m) != nrow(o))||(ncol(n) != ncol(o)))
	{
		LATAN_ERROR("operation between matrices with dimension mismatch",\
					LATAN_EBADLEN);
	}
	
	if (have_duplicate_arg)
	{
		buf = mat_create_from_dim(m);
		pt  = buf;
	}
	else
	{
		pt = m;
	}
	LATAN_UPDATE_STATUS(status,gsl_blas_dgemm(CblasNoTrans,CblasTrans,\
											  1.0,n,o,0.0,pt));
	if (have_duplicate_arg)
	{
		LATAN_UPDATE_STATUS(status,mat_cp(m,buf));
		mat_destroy(buf);
	}
	
	return status;
}

latan_errno mat_mul_tn(mat m, const mat n, const mat o)
{
	latan_errno status;
	bool have_duplicate_arg;
	mat pt,buf;
	
	status             = LATAN_SUCCESS;
	have_duplicate_arg = ((m == n)||(m == o));
	
	if ((nrow(m) != ncol(n))||(ncol(m) != ncol(o))||(nrow(n) != nrow(o)))
	{
		LATAN_ERROR("operation between matrices with dimension mismatch",\
					LATAN_EBADLEN);
	}
	
	if (have_duplicate_arg)
	{
		buf = mat_create_from_dim(m);
		pt  = buf;
	}
	else
	{
		pt = m;
	}
	LATAN_UPDATE_STATUS(status,gsl_blas_dgemm(CblasTrans,CblasNoTrans,\
											  1.0,n,o,0.0,pt));
	if (have_duplicate_arg)
	{
		LATAN_UPDATE_STATUS(status,mat_cp(m,buf));
		mat_destroy(buf);
	}
	
	return status;
}

latan_errno mat_mul_tt(mat m, const mat n, const mat o)
{
	latan_errno status;
	bool have_duplicate_arg;
	mat pt,buf;
	
	status             = LATAN_SUCCESS;
	have_duplicate_arg = ((m == n)||(m == o));
	
	if ((nrow(m) != ncol(n))||(ncol(m) != nrow(o))||(nrow(n) != ncol(o)))
	{
		LATAN_ERROR("operation between matrices with dimension mismatch",\
					LATAN_EBADLEN);
	}
	
	if (have_duplicate_arg)
	{
		buf = mat_create_from_dim(m);
		pt  = buf;
	}
	else
	{
		pt = m;
	}
	LATAN_UPDATE_STATUS(status,gsl_blas_dgemm(CblasTrans,CblasTrans,\
											  1.0,n,o,0.0,pt));
	if (have_duplicate_arg)
	{
		LATAN_UPDATE_STATUS(status,mat_cp(m,buf));
		mat_destroy(buf);
	}
	
	return status;
}

latan_errno mat_eqtranspose(mat m)
{
	latan_errno status;
	
	if (!mat_is_square(m))
	{
		LATAN_ERROR("cannot auto-transpose a non-square matrix",LATAN_ENOTSQR);
	}
	
	status = gsl_matrix_transpose(m);
	
	return status;
}

latan_errno mat_transpose(mat m, const mat n)
{
	latan_errno status;
	
	if ((nrow(m) != ncol(n))||(ncol(m) != nrow(n)))
	{
		LATAN_ERROR("operation between matrices with dimension mismatch",\
					LATAN_EBADLEN);
	}
	
	status = gsl_matrix_transpose_memcpy(m,n);
	
	return status;
}

latan_errno mat_inv(mat m, const mat n)
{
	latan_errno status;
	int signum;
	mat LU;
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
	
	LATAN_UPDATE_STATUS(status,gsl_linalg_LU_decomp(LU,perm,&signum));
	for (i=0;i<nrow(LU);i++)
	{
		if (mat_get(LU,i,i) == 0.0)
		{
			LATAN_ERROR("trying to invert a singular matrix",LATAN_EDOM);
		}
	}
	LATAN_UPDATE_STATUS(status,gsl_linalg_LU_invert(LU,perm,m));
	
	mat_destroy(LU);
	gsl_permutation_free(perm);
	
	return status;
}

latan_errno mat_eqmulp(mat m, const mat n)
{
	latan_errno status;
	
	status = LATAN_SUCCESS;
	
	if (!(mat_is_samedim(m,n)))
	{
		LATAN_ERROR("operation between matrices with dimension mismatch",\
					LATAN_EBADLEN);
	}
	
	LATAN_UPDATE_STATUS(status,gsl_matrix_mul_elements(m,n));
	
	return status;
}

latan_errno mat_mulp(mat m, const mat n, const mat o)
{
	latan_errno status;
	
	status = LATAN_SUCCESS;
	
	LATAN_UPDATE_STATUS(status,mat_cp(m,n));
	LATAN_UPDATE_STATUS(status,mat_eqmulp(m,o));
	
	return status;
}

latan_errno mat_eqmuls(mat m, const double s)
{
	latan_errno status;
	
	status = LATAN_SUCCESS;
	
	LATAN_UPDATE_STATUS(status,gsl_matrix_scale(m,s));
	
	return status;
}

latan_errno mat_muls(mat m, const mat n, const double s)
{
	latan_errno status;
	
	status = LATAN_SUCCESS;
	
	LATAN_UPDATE_STATUS(status,mat_cp(m,n));
	LATAN_UPDATE_STATUS(status,mat_eqmuls(m,s));
	
	return status;
}

latan_errno mat_eqdivp(mat m, const mat n)
{
	latan_errno status;
	
	status = LATAN_SUCCESS;
	
	if (!(mat_is_samedim(m,n)))
	{
		LATAN_ERROR("operation between matrices with dimension mismatch",\
					LATAN_EBADLEN);
	}
	
	LATAN_UPDATE_STATUS(status,gsl_matrix_div_elements(m,n));
	
	return status;
}

latan_errno mat_divp(mat m, const mat n, const mat o)
{
	latan_errno status;
	
	status = LATAN_SUCCESS;
	
	LATAN_UPDATE_STATUS(status,mat_cp(m,n));
	LATAN_UPDATE_STATUS(status,mat_eqdivp(m,o));
	
	return status;
}

latan_errno mat_abs(mat m, const mat n)
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

latan_errno mat_sqrt(mat m, const mat n)
{
	size_t i,j;
	
	if (!mat_is_samedim(m,n))
	{
		LATAN_ERROR("operation between matrices with dimension mismatch",\
					LATAN_EBADLEN);
	}
	
	FOR_VAL(m,i,j)
	{
		mat_set(m,i,j,sqrt(mat_get(n,i,j)));
	}
	
	return LATAN_SUCCESS;
}
