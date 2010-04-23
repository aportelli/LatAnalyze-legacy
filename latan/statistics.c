#include <latan/statistics.h>
#include <latan/includes.h>
#include <latan/rand.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_errno.h>

static latan_errno resample_bootstrap(mat cent_val, mat* sample,			\
									  const size_t nboot, const mat* dat,	\
									  const size_t ndat, rs_func* f,		\
									  void* param);
static latan_errno resample_jackknife(mat cent_val, mat* sample,			\
									  const int jk_depth, const mat* dat,	\
									  const size_t ndat, rs_func* f,		\
									  void* param);
static size_t jackknife_nsample(const size_t ndat, const size_t jk_depth); 

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

latan_errno mat_mean(mat mean, const mat *m, const size_t size)
{
	latan_errno status;
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

latan_errno mat_cov(mat cov, const mat *m, const mat *n, const size_t size)
{
	latan_errno status;
	mat m_mean, n_mean;
	
	status = LATAN_SUCCESS;
	
	m_mean = mat_create_from_dim(m[0]);
	n_mean = mat_create_from_dim(n[0]);
	
	LATAN_UPDATE_STATUS(status,mat_mean(m_mean,m,size));
	LATAN_UPDATE_STATUS(status,mat_mean(n_mean,n,size));
	LATAN_UPDATE_STATUS(status,mat_cov_m(cov,m,n,size,m_mean,n_mean));
	
	mat_destroy(m_mean);
	mat_destroy(n_mean);
	
	return status;
}

latan_errno mat_cov_m(mat cov, const mat* m, const mat* n, const size_t size,\
					  const mat m_mean, const mat n_mean)
{
	latan_errno status;
	size_t i;
	size_t subdim;
	double dsubdim;
	mat* mctnc;
	mat* mc;
	mat* nc;
	
	status = LATAN_SUCCESS;
	subdim = ncol(m[0]);
	dsubdim = (double)(subdim);
	
	mctnc = mat_create_ar_from_dim(size,cov);
	mc = mat_create_ar_from_dim(size,m[0]);
	nc = mat_create_ar_from_dim(size,n[0]);
	
	for (i=0;i<size;i++) 
	{
		LATAN_UPDATE_STATUS(status,mat_sub(mc[i],m[i],m_mean));
		LATAN_UPDATE_STATUS(status,mat_sub(nc[i],n[i],n_mean));
		LATAN_UPDATE_STATUS(status,mat_mul_nt(mctnc[i],mc[i],nc[i]));
		LATAN_UPDATE_STATUS(status,mat_eqmuls(mctnc[i],1.0/dsubdim));
	}
	LATAN_UPDATE_STATUS(status,mat_mean(cov,mctnc,size));
	
	mat_destroy_ar(mctnc,size);
	mat_destroy_ar(mc,size);
	mat_destroy_ar(nc,size);
	
	return status;
}

latan_errno mat_covp(mat cov, const mat *m, const mat *n, const size_t size)
{
	latan_errno status;
	mat m_mean, n_mean;
	
	status = LATAN_SUCCESS;
	
	m_mean = mat_create_from_dim(m[0]);
	n_mean = mat_create_from_dim(n[0]);
	
	LATAN_UPDATE_STATUS(status,mat_mean(m_mean,m,size));
	LATAN_UPDATE_STATUS(status,mat_mean(n_mean,n,size));
	LATAN_UPDATE_STATUS(status,mat_covp_m(cov,m,n,size,m_mean,n_mean));
	
	mat_destroy(m_mean);
	mat_destroy(n_mean);
	
	return status;
}

latan_errno mat_covp_m(mat cov, const mat* m, const mat* n, const size_t size,\
					   const mat m_mean, const mat n_mean)
{
	latan_errno status;
	size_t i;
	mat* mcnc;
	mat* nc;
	
	status = LATAN_SUCCESS;
	
	mcnc = mat_create_ar_from_dim(size,cov);
	nc = mat_create_ar_from_dim(size,cov);
	
	for (i=0;i<size;i++)
	{
		LATAN_UPDATE_STATUS(status,mat_sub(mcnc[i],m[i],m_mean));
		LATAN_UPDATE_STATUS(status,mat_sub(nc[i],n[i],n_mean));
		LATAN_UPDATE_STATUS(status,mat_eqmulp(mcnc[i],nc[i]));
	}
	LATAN_UPDATE_STATUS(status,mat_mean(cov,mcnc,size));
	
	mat_destroy_ar(mcnc,size);
	mat_destroy_ar(nc,size);

	return status;
}

/*								histogram									*/
/****************************************************************************/
latan_errno histogram(mat hist, const mat data, const double xmin,\
					  const double xmax, const size_t nint)
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

/*					resampled samples manipulation							*/
/****************************************************************************/
/** jackknife sample number calculation **/
static size_t jackknife_nsample(const size_t ndat, const size_t jk_depth)
{
	unsigned int i;
	size_t nsample;
	const unsigned int ui_ndat = (unsigned int)(ndat);
	
	nsample = 0;

	for (i=1;i<=jk_depth;i++)
	{
		nsample += (size_t)(latan_binomial(ui_ndat,i));
	}
	
	return nsample;
}

/** allocation **/
rs_sample rs_sample_create_boot(const size_t init_nrow, const size_t nboot,\
								const stringbuf name)
{
	rs_sample s;
	
	MALLOC_ERRVAL(s,rs_sample,1,NULL);
	
	s->nsample = nboot;
	strcpy(s->name,name);
	s->resamp_method = BOOT;
	
	s->cent_val = mat_create(init_nrow,1);
	s->sample = mat_create_ar_from_dim(s->nsample,s->cent_val);
	
	return s;
}

rs_sample rs_sample_create_jack(const size_t init_nrow, const size_t ndat,\
								const size_t jk_depth, const stringbuf name)
{
	rs_sample s;
	
	MALLOC_ERRVAL(s,rs_sample,1,NULL);
	
	s->nsample = jackknife_nsample(ndat,jk_depth);
	strcpy(s->name,name);
	s->resamp_method = JACK;
	
	s->cent_val = mat_create(init_nrow,1);
	s->sample = mat_create_ar_from_dim(s->nsample,s->cent_val);
	
	return s;
}

void rs_sample_destroy(rs_sample s)
{
	mat_destroy(s->cent_val);
	mat_destroy_ar(s->sample,s->nsample);
	s->nsample = 0;
}

/** access **/
mat rs_sample_get_sample(const rs_sample s, const size_t i)
{
	return (s->sample)[i];
}

/** estimators **/
latan_errno rs_sample_cov(mat cov, const rs_sample s, const rs_sample t)
{
	latan_errno status;
	
	if (s->nsample != t->nsample)
	{
		LATAN_ERROR("operation between resampled samples with dimension mismatch",\
					LATAN_EBADLEN);
	}
	
	status = mat_cov_m(cov,s->sample,t->sample,s->nsample,\
					   s->cent_val,t->cent_val);
	
	return status;
}

latan_errno rs_sample_covp(mat cov, const rs_sample s, const rs_sample t)
{
	latan_errno status;
	
	if (s->nsample != t->nsample)
	{
		LATAN_ERROR("operation between resampled samples with dimension mismatch",\
					LATAN_EBADLEN);
	}
	
	status = mat_covp_m(cov,s->sample,t->sample,s->nsample,\
						s->cent_val,t->cent_val);
	
	return status;
}
/*						resampling functions								*/
/****************************************************************************/
static latan_errno resample_bootstrap(mat cent_val, mat* sample,			\
									  const size_t nboot, const mat* dat,	\
									  const size_t ndat, rs_func* f,		\
									  void* param)
{
	mat* fakedat;
	size_t i,j;
	unsigned int rj;
	latan_errno status;
	
	status = LATAN_SUCCESS;
	
	fakedat = mat_create_ar_from_dim(ndat,dat[0]);
	
	LATAN_UPDATE_STATUS(status,f(cent_val,dat,ndat,0,param));
	for (i=0;i<nboot;i++)
	{
		for (j=0;j<ndat;j++) 
		{
			rj = rand_ud((unsigned int)(ndat));
			fakedat[j] = dat[rj];
		}
		LATAN_UPDATE_STATUS(status,f(sample[i],fakedat,ndat,i+1,param));
	}
	
	return status;
}

latan_errno resample(rs_sample s, const mat* dat, size_t ndat,\
					 rs_func* f, void* param)
{
	latan_errno status;
	
	switch (s->resamp_method) 
	{
		case BOOT:
			randgen_get_state(s->gen_state);
			status = resample_bootstrap(s->cent_val,s->sample,s->nsample,dat,\
										ndat,f,param);
			break;
		default:
			LATAN_ERROR("resampling method flag invalid",LATAN_EINVAL);
			break;
	}
	
	return status;
}

/* some basic rs_func */
latan_errno rs_mean(mat res, const mat* dat, const size_t ndat,\
					const size_t sampno, void* param)
{
	latan_errno status;
	size_t st_nothing;
	
	param = NULL;
	st_nothing = sampno;
	
	status = mat_mean(res,dat,ndat);
	
	return status;
}