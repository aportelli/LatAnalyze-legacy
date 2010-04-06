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
int mat_mean(mat mean, const mat *m, const size_t size);
/*!
 @fn void mat_meansig(mat mean, mat sig, const mat *m, const int size)
 @brief Compute the "coefficient by coefficient" mean and standard deviation of an array of matrices.
 
 @param mean matrix to store the mean
 @param sig matrix to store the standard deviation
 @param m input array of matrices
 @param size number of matrices in the array
 @warning \b mean, \b sig and all the elements of \b m must have the same dimensions, fatal error will be generated if it is not the case.
 */
int mat_meansig(mat mean, mat sig, const mat *m, const size_t size);

/* resampling functions */
int resamp_bootstrap(mat res, mat sigres, mat* dat, size_t ndat,\
					 estimator_t* estimator, void* param,		\
					 size_t nboot, const boot_io* io_param);
int resamp_jackknife(mat res, mat sigres, mat* dat, size_t ndat,\
					 estimator_t* estimator, void* param);

extern const boot_io boot_no_io;

#endif