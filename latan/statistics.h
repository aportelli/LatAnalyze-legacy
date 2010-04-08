#ifndef LATAN_STATISTICS_H_
#define LATAN_STATISTICS_H_

#include <latan/mat.h>

typedef int estimator_t(mat res, mat* dat, size_t ndat, void* param);

typedef struct
{
	bool save_samples;
	stringbuf fname;
	stringbuf est_name;
	mat samples_buf;
} boot_io;

/* elementary estimators */
double mat_elsum(mat m);
double mat_elmean(mat m);
/*!
 @fn void mat_mean(mat mean, const mat *m, const int size)
 @brief Compute the "coefficient by coefficient" mean of an array of matrices.
 @see #mat_print
 
 @param mean matrix to store the mean
 @param m input array of matrices
 @param size number of matrices in the array
 @warning \b mean and all the elements of \b m must have the same dimensions, fatal error will be generated if it is not the case.
 */
int mat_mean(mat mean, const mat* m, const size_t size);
/*!
 @fn void mat_meansig(mat mean, mat sig, const mat *m, const int size)
 @brief Compute the "coefficient by coefficient" mean and standard deviation of an array of matrices.
 
 @param mean matrix to store the mean
 @param sig matrix to store the standard deviation
 @param m input array of matrices
 @param size number of matrices in the array
 @warning \b mean, \b sig and all the elements of \b m must have the same dimensions, fatal error will be generated if it is not the case.
 */
int mat_cov(mat cov, const mat* m, const mat* n, const size_t size);
int mat_cov_m(mat cov, const mat* m, const mat* n, const size_t size,\
			  const mat m_mean, const mat n_mean);
#define mat_var(var,m,size) mat_cov(var,m,m,size)
#define mat_var_m(var,m,size,mean) mat_cov(var,m,m,size,mean,mean)
int mat_covp(mat cov, const mat* m, const mat* n, const size_t size);
int mat_covp_m(mat cov, const mat* m, const mat* n, const size_t size,\
			   const mat m_mean, const mat n_mean);
#define mat_varp(var,m,size) mat_covar(var,m,m,size)
#define mat_varp_m(var,m,size,mean) mat_covar_m(var,m,m,size,mean,mean)

/* resampling functions */
int resamp_bootstrap(boot_res* res, const mat* dat, size_t ndat,\
					 estimator_t* estimator, void* param);
int resamp_jackknife(jack_res* sample, const mat* dat, size_t ndat,	\
					 estimator_t* estimator, void* param);

#endif