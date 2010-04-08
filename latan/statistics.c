#include <latan/statistics.h>
#include <latan/includes.h>
#include <latan/rand.h>

const boot_io boot_no_io = {false,"","",NULL};

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
	
	mat_create_from_dim(&m_mean,cov);
	mat_create_from_dim(&n_mean,cov);
	
	LATAN_UPDATE_STATUS(status,mat_mean(m_mean,m,size));
	LATAN_UPDATE_STATUS(status,mat_mean(n_mean,n,size));
	LATAN_UPDATE_STATUS(status,mat_cov_m(cov,m,n,size,m_mean,n_mean));
	
	mat_destroy(&m_mean);
	mat_destroy(&n_mean);
	
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
	
	mat_create_ar_from_dim(&mtn,size,cov);
	mat_create_from_dim(&mean_tprod,cov);
	
	for (i=0;i<size;i++) 
	{
		LATAN_UPDATE_STATUS(status,mat_mul_nt(mtn[i],m[i],n[i]));
		LATAN_UPDATE_STATUS(status,mat_eqmuls(mtn[i],1.0/dsubdim));
	}
	LATAN_UPDATE_STATUS(status,mat_mul_nt(mean_tprod,m_mean,n_mean));
	LATAN_UPDATE_STATUS(status,mat_mean(cov,mtn,size));
	LATAN_UPDATE_STATUS(status,mat_eqsub(cov,mean_tprod));
	
	mat_destroy_ar(&mtn,size);
	mat_destroy(&mean_tprod);
	
	return status;
}

int mat_covp(mat cov, const mat *m, const mat *n, const size_t size)
{
	int status;
	mat m_mean, n_mean;
	
	status = LATAN_SUCCESS;
	
	mat_create_from_dim(&m_mean,cov);
	mat_create_from_dim(&n_mean,cov);
	
	LATAN_UPDATE_STATUS(status,mat_mean(m_mean,m,size));
	LATAN_UPDATE_STATUS(status,mat_mean(n_mean,n,size));
	LATAN_UPDATE_STATUS(status,mat_covp_m(cov,m,n,size,m_mean,n_mean));
	
	mat_destroy(&m_mean);
	mat_destroy(&n_mean);
	
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
	
	mat_create_ar_from_dim(&mn,size,cov);
	mat_create_from_dim(&mean_prod,cov);
	
	for (i=0;i<size;i++)
	{
		LATAN_UPDATE_STATUS(status,mat_mulp(mn[i],m[i],n[i]));
	}
	LATAN_UPDATE_STATUS(status,mat_mulp(mean_prod,m_mean,n_mean));
	LATAN_UPDATE_STATUS(status,mat_mean(cov,mn,size));
	LATAN_UPDATE_STATUS(status,mat_eqsub(cov,mean_prod));
	
	mat_destroy_ar(&mn,size);
	mat_destroy(&mean_prod);
	
	return status;
}