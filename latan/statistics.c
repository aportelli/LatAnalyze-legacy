#include <latan/statistics.h>
#include <latan/includes.h>
#include <latan/rand.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_errno.h>

/*						elementary estimators								*/
/****************************************************************************/
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

int mat_mean(mat mean, const mat *m, const size_t size)
{
	int status;
	size_t i;
	const double dsize = (double)(size);
	
	status = LATAN_SUCCESS;
	mat_zero(mean);
	
	for (i=0;i<size;i++)
	{
		LATAN_UPDATE_STATUS(status,mat_eqadd(mean,m[i]));
	}
	LATAN_UPDATE_STATUS(status,mat_eqmuls(mean,1.0/dsize));
	
	return status;
}

int mat_cov(mat cov, const mat *m, const mat *n, const size_t size)
{
	int status;
	mat m_mean, n_mean;
	
	status = LATAN_SUCCESS;
	
	m_mean = mat_create_from_dim(cov);
	n_mean = mat_create_from_dim(cov);
	
	LATAN_UPDATE_STATUS(status,mat_mean(m_mean,m,size));
	LATAN_UPDATE_STATUS(status,mat_mean(n_mean,n,size));
	LATAN_UPDATE_STATUS(status,mat_cov_m(cov,m,n,size,m_mean,n_mean));
	
	mat_destroy(m_mean);
	mat_destroy(n_mean);
	
	return status;
}

int mat_cov_m(mat cov, const mat* m, const mat* n, const size_t size,\
			  const mat m_mean, const mat n_mean)
{
	int status;
	size_t i;
	size_t subdim;
	double dsubdim;
	mat* mtn;
	mat mean_tprod;
	
	status = LATAN_SUCCESS;
	subdim = ncol(m[0]);
	dsubdim = (double)(subdim);
	
	mtn = mat_create_ar_from_dim(size,cov);
	mean_tprod = mat_create_from_dim(cov);
	
	for (i=0;i<size;i++) 
	{
		LATAN_UPDATE_STATUS(status,mat_mul_nt(mtn[i],m[i],n[i]));
		LATAN_UPDATE_STATUS(status,mat_eqmuls(mtn[i],1.0/dsubdim));
	}
	LATAN_UPDATE_STATUS(status,mat_mul_nt(mean_tprod,m_mean,n_mean));
	LATAN_UPDATE_STATUS(status,mat_mean(cov,mtn,size));
	LATAN_UPDATE_STATUS(status,mat_eqsub(cov,mean_tprod));
	
	mat_destroy_ar(mtn,size);
	mat_destroy(mean_tprod);
	
	return status;
}

int mat_covp(mat cov, const mat *m, const mat *n, const size_t size)
{
	int status;
	mat m_mean, n_mean;
	
	status = LATAN_SUCCESS;
	
	m_mean = mat_create_from_dim(cov);
	n_mean = mat_create_from_dim(cov);
	
	LATAN_UPDATE_STATUS(status,mat_mean(m_mean,m,size));
	LATAN_UPDATE_STATUS(status,mat_mean(n_mean,n,size));
	LATAN_UPDATE_STATUS(status,mat_covp_m(cov,m,n,size,m_mean,n_mean));
	
	mat_destroy(m_mean);
	mat_destroy(n_mean);
	
	return status;
}

int mat_covp_m(mat cov, const mat* m, const mat* n, const size_t size,\
				const mat m_mean, const mat n_mean)
{
	int status;
	size_t i;
	mat* mn;
	mat mean_prod;
	
	status = LATAN_SUCCESS;
	
	mn = mat_create_ar_from_dim(size,cov);
	mean_prod = mat_create_from_dim(cov);
	
	for (i=0;i<size;i++)
	{
		LATAN_UPDATE_STATUS(status,mat_mulp(mn[i],m[i],n[i]));
	}
	LATAN_UPDATE_STATUS(status,mat_mulp(mean_prod,m_mean,n_mean));
	LATAN_UPDATE_STATUS(status,mat_mean(cov,mn,size));
	LATAN_UPDATE_STATUS(status,mat_eqsub(cov,mean_prod));
	
	mat_destroy_ar(mn,size);
	mat_destroy(mean_prod);
	
	return status;
}

/*								histogram									*/
/****************************************************************************/
int histogram(mat hist, const mat data, const double xmin, const double xmax,\
			  const size_t nint)
{
	size_t i;
	gsl_histogram* gsl_hist;
	gsl_error_handler_t* error_handler;
	
	
	if (nrow(hist) != nint)
	{
		LATAN_ERROR("histogram matrix row number do not match number of histogram intervals",\
					LATAN_EBADLEN);
	}
	
	error_handler = gsl_set_error_handler(&latan_error);
	gsl_hist = gsl_histogram_alloc(nint);
	gsl_set_error_handler(error_handler);
	mat_zero(hist);
	
	gsl_histogram_set_ranges_uniform(gsl_hist,xmin,xmax);
	for (i=0;i<nrow(data);i++)
	{
		gsl_histogram_increment(gsl_hist,mat_get(data,i,0));
	}
	mat_set_from_ar(hist,gsl_hist->bin);
	
	gsl_histogram_free(gsl_hist);
	
	return LATAN_SUCCESS;
}
