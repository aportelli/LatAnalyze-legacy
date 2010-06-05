#ifndef LATAN_STATISTICS_H_
#define LATAN_STATISTICS_H_

#include <latan/latan_globals.h>

#define BOOT 0
#define JACK 1
#define GENERIC 2

__BEGIN_DECLS

typedef struct
{
	stringbuf name;
	mat cent_val;
	mat* sample;
	size_t nsample;
	int resamp_method;
	randgen_state gen_state;
}* rs_sample;

typedef latan_errno rs_func(mat res, const mat* dat, const size_t ndat,\
							const size_t sampno, void* param);

/* elementary estimators */
double mat_elsum(mat m);
double mat_elmean(mat m);
latan_errno mat_mean(mat mean, const mat* m, const size_t size);
latan_errno mat_cov(mat cov, const mat* m, const mat* n, const size_t size);
latan_errno mat_cov_m(mat cov, const mat* m, const mat* n, const size_t size,\
					  const mat m_mean, const mat n_mean);
#define mat_var(var,m,size) mat_cov(var,m,m,size)
#define mat_var_m(var,m,size,mean) mat_cov(var,m,m,size,mean,mean)
latan_errno mat_covp(mat cov, const mat* m, const mat* n, const size_t size);
latan_errno mat_covp_m(mat cov, const mat* m, const mat* n, const size_t size,\
					   const mat m_mean, const mat n_mean);
#define mat_varp(var,m,size) mat_covp(var,m,m,size)
#define mat_varp_m(var,m,size,mean) mat_covp_m(var,m,m,size,mean,mean)

/* data binning */
latan_errno mat_ar_bin(mat* bindat, const mat* dat, const size_t ndat,\
					   const size_t binsize);

/* histogram */
latan_errno histogram(mat hist, const mat data, const double xmin,\
					  const double xmax, const size_t nint);

/* resampled samples manipulation */
/** allocation **/
rs_sample rs_sample_create_boot(const size_t init_nrow, const size_t nboot,\
								const stringbuf name);
rs_sample rs_sample_create_jack(const size_t init_nrow, const size_t ndat,\
								const size_t jk_depth, const stringbuf name);
rs_sample rs_sample_create(const size_t init_nrow, const size_t nsample,\
						   const stringbuf name);
void rs_sample_destroy(rs_sample s);

/** access **/
size_t rs_sample_get_nsample(const rs_sample s);
int rs_sample_get_method(const rs_sample s);
mat rs_sample_pt_cent_val(const rs_sample s);
mat rs_sample_pt_sample(const rs_sample s, const size_t i);

/** estimators **/
latan_errno rs_sample_cov(mat cov, const rs_sample s, const rs_sample t);
#define rs_sample_var(cov,s) rs_sample_cov(cov,s,s);
latan_errno rs_sample_covp(mat cov, const rs_sample s, const rs_sample t);
#define rs_sample_varp(cov,s) rs_sample_covp(cov,s,s);

/* resampling function */
latan_errno resample(rs_sample s, const mat* dat, const size_t ndat,\
					 const size_t nobs, rs_func* f, void* param);

/* useful rs_func */
latan_errno rs_mean(mat res, const mat* dat, const size_t ndat,\
					const size_t sampno, void* nothing);
latan_errno rs_finite_diff(mat res, const mat* dat, const size_t ndat,\
						   const size_t sampno, void* nothing);
latan_errno rs_effmass(mat res, const mat* dat, const size_t ndat,\
					   const size_t sampno, void* parity);
latan_errno rs_effmass_PCAC(mat res, const mat* dat, const size_t ndat,\
							const size_t sampno, void* nothing);
latan_errno rs_data_fit(mat res, const mat* dat, const size_t ndat,\
							const size_t sampno, void* d);

__END_DECLS

#endif