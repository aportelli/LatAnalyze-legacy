#include <latan/mat.h>
#include <latan/includes.h>
#include <latan/rand.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>

/*								allocation									*/
/****************************************************************************/
void mat_create(mat* m, const size_t init_nrow, const size_t init_ncol)
{
	gsl_error_handler_t* error_handler;
	
	error_handler = gsl_set_error_handler(&latan_error);
	*m = gsl_matrix_alloc(init_nrow,init_ncol);
	gsl_set_error_handler(error_handler);
}

void mat_create_from_mat(mat* m, const mat n)
{
	mat_create(m,nrow(n),ncol(n));
	mat_cp(*m,n);
}

void mat_create_ar(mat** m, const size_t nmat, const size_t init_nrow,\
			   const size_t init_ncol)
{
	size_t i;
	
	MALLOC_NOERRET(*m,mat*,nmat);
	for (i=0;i<nmat;i++)
	{
		mat_create((*m)+i,init_nrow,init_ncol);
	}
}

void mat_destroy(mat* m)
{
	gsl_matrix_free(*m);
}

void mat_destroy_ar(mat** m, const size_t nmat)
{
	size_t i;
	
	for (i=0;i<nmat;i++)
	{
		mat_destroy((*m)+i);
	}
	FREE(*m);
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

int mat_cp_subm(mat m, const mat n, const size_t k1, const size_t l1, \
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

int mat_set_subm(mat m, const mat n, const size_t k1, const size_t l1, \
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

/* tests */
bool mat_issamedim(const mat m, const mat n)
{
	return ((nrow(m) == nrow(n))&&(ncol(m) == ncol(n)));
}

bool mat_issquare(const mat m)
{
	return (nrow(m) == ncol(m));
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

int mat_cp(mat m, const mat n)
{
	int status;
	
	status = LATAN_SUCCESS;
	
	if (!mat_issamedim(m,n))
	{
		LATAN_ERROR("matrix copy with dimension mismatch",LATAN_EBADLEN);
	}
	
	LATAN_UPDATE_STATUS(status,gsl_matrix_memcpy(m,n));
	
	return status;
}
	
int mat_eqadd(mat m, const mat n)
{
	int status;
	
	status = LATAN_SUCCESS;
	
	if (!(mat_issamedim(m,n)))
	{
		LATAN_ERROR("operation between matrices with dimension mismatch",\
					LATAN_EBADLEN);
	}
	
	LATAN_UPDATE_STATUS(status,gsl_matrix_add(m,n));
	
	return status;
}

int mat_add(mat m, const mat n, const mat o)
{
	int status;
	
	status = LATAN_SUCCESS;
	
	LATAN_UPDATE_STATUS(status,mat_cp(m,n));
	LATAN_UPDATE_STATUS(status,mat_eqadd(m,o));
	
	return status;
}


int mat_eqsub(mat m, const mat n)
{
	int status;
	
	status = LATAN_SUCCESS;
	
	if (!(mat_issamedim(m,n)))
	{
		LATAN_ERROR("operation between matrices with dimension mismatch",\
					LATAN_EBADLEN);
	}
	
	LATAN_UPDATE_STATUS(status,gsl_matrix_sub(m,n));
	
	return status;
}

int mat_sub(mat m, const mat n, const mat o)
{
	int status;
	
	status = LATAN_SUCCESS;
	
	LATAN_UPDATE_STATUS(status,mat_cp(m,n));
	LATAN_UPDATE_STATUS(status,mat_eqsub(m,o));
	
	return status;
}

int mat_eqmul_l(mat m, const mat n)
{
	int status;
	mat buf;
	
	status = LATAN_SUCCESS;
	
	mat_create(&buf,nrow(m),ncol(m));
	
	LATAN_UPDATE_STATUS(status,mat_mul(buf,n,m));
	LATAN_UPDATE_STATUS(status,mat_cp(m,buf));
	
	return status;
}

int mat_eqmul_r(mat m, const mat n)
{
	int status;
	mat buf;
	
	status = LATAN_SUCCESS;
	
	mat_create(&buf,nrow(m),ncol(m));
	
	LATAN_UPDATE_STATUS(status,mat_mul(buf,m,n));
	LATAN_UPDATE_STATUS(status,mat_cp(m,buf));
	
	mat_destroy(&buf);
	
	return status;
}

int mat_mul(mat m, const mat n, const mat o)
{
	int status;
	
	status = LATAN_SUCCESS;
	
	if (!((nrow(m) == nrow(n))&&(ncol(m) == ncol(o))&&(ncol(n) == nrow(o))))
	{
		LATAN_ERROR("operation between matrices with dimension mismatch",\
					LATAN_EBADLEN);
	}
	if ((m == n)||(m == o))
	{
		LATAN_ERROR("output argument of mat_mul can not be also an input argument (BLAS limitation), use mat_eqmul_* functions",\
					LATAN_EINVAL);
	}
	
	LATAN_UPDATE_STATUS(status,gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,\
											  1.0,n,o,0.0,m));
	
	return status;
}

int mat_eqmulp(mat m, const mat n)
{
	int status;
	
	status = LATAN_SUCCESS;
	
	if (!(mat_issamedim(m,n)))
	{
		LATAN_ERROR("operation between matrices with dimension mismatch",\
					LATAN_EBADLEN);
	}
	
	LATAN_UPDATE_STATUS(status,gsl_matrix_mul_elements(m,n));
	
	return status;
}

int mat_mulp(mat m, const mat n, const mat o)
{
	int status;
	
	status = LATAN_SUCCESS;
	
	LATAN_UPDATE_STATUS(status,mat_cp(m,n));
	LATAN_UPDATE_STATUS(status,mat_eqmulp(m,o));
	
	return status;
}

int mat_eqmuls(mat m, const double s)
{
	int status;
	
	status = LATAN_SUCCESS;
	
	LATAN_UPDATE_STATUS(status,gsl_matrix_scale(m,s));
	
	return status;
}

int mat_muls(mat m, const mat n, const double s)
{
	int status;
	
	status = LATAN_SUCCESS;
	
	LATAN_UPDATE_STATUS(status,mat_cp(m,n));
	LATAN_UPDATE_STATUS(status,mat_eqmuls(m,s));
	
	return status;
}

int mat_eqdivp(mat m, const mat n)
{
	int status;
	
	status = LATAN_SUCCESS;
	
	if (!(mat_issamedim(m,n)))
	{
		LATAN_ERROR("operation between matrices with dimension mismatch",\
					LATAN_EBADLEN);
	}
	
	LATAN_UPDATE_STATUS(status,gsl_matrix_div_elements(m,n));
	
	return status;
}

int mat_divp(mat m, const mat n, const mat o)
{
	int status;
	
	status = LATAN_SUCCESS;
	
	LATAN_UPDATE_STATUS(status,mat_cp(m,n));
	LATAN_UPDATE_STATUS(status,mat_eqdivp(m,o));
	
	return status;
}

int mat_abs(mat m, const mat n)
{
	size_t i,j;
	
	if (!mat_issamedim(m,n))
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

int mat_sqrt(mat m, const mat n)
{
	size_t i,j;
	
	if (!mat_issamedim(m,n))
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

int mat_mean(mat mean, const mat *m, const size_t size)
{
	int status;
	size_t i;
	
	status = LATAN_SUCCESS;
	mat_zero(mean);
	
	for (i=0;i<size;i++)
	{
		LATAN_UPDATE_STATUS(status,mat_eqadd(mean,m[i]));
	}
	LATAN_UPDATE_STATUS(status,mat_eqmuls(mean,1.0/((double)size)));
	
	return status;
}

int mat_meansig(mat mean, mat sig, const mat *m, const size_t size)
{
	int status;
	size_t i;
	const double dsize = ((double)(size));
	mat sq;
	
	if (!mat_issamedim(mean,sig))
	{
		LATAN_ERROR("operation between matrices with dimension mismatch",\
					LATAN_EBADLEN);
	}
	
	status = LATAN_SUCCESS;
	mat_zero(mean);
	mat_zero(sig);
	
	mat_create(&sq,nrow(mean),ncol(mean));
	
	for (i=0;i<size;i++)
	{
		LATAN_UPDATE_STATUS(status,mat_eqadd(mean,m[i]));
		LATAN_UPDATE_STATUS(status,mat_mulp(sq,m[i],m[i]));
		LATAN_UPDATE_STATUS(status,mat_eqadd(sig,sq));
	}
	LATAN_UPDATE_STATUS(status,mat_eqmuls(mean,1.0/dsize));
	LATAN_UPDATE_STATUS(status,mat_eqmuls(sig,1.0/(dsize-1.0)));
	LATAN_UPDATE_STATUS(status,mat_mulp(sq,mean,mean));
	LATAN_UPDATE_STATUS(status,mat_eqmuls(sq,(dsize)/(dsize-1.0)));
	LATAN_UPDATE_STATUS(status,mat_eqsub(sig,sq));
	LATAN_UPDATE_STATUS(status,mat_eqsqrt(sig));
	
	mat_destroy(&sq);
	
	return status;
}

double mat_elsum(mat m)
{
	size_t i,j;
	double sum;
	
	sum = 0.0;
	
	FOR_VAL(m,i,j)
	{
		sum += mat_get(m,i,j);
	}
	
	return sum;
}

double mat_elmean(mat m)
{
	double mean;
	
	mean = mat_elsum(m)/((double)(nrow(m)*ncol(m)));
	
	return mean;
}

/*									I/O										*/
/****************************************************************************/
void mat_dump(FILE *stream, const mat m)
{
	size_t i,j;
	
	for (i=0;i<nrow(m);i++)
	{
		for (j=0;j<ncol(m)-1;j++)
		{
			fprintf(stream,"%.10e ",mat_get(m,i,j));
		}
		fprintf(stream,"%.10e\n",mat_get(m,i,ncol(m)-1));
	}
}