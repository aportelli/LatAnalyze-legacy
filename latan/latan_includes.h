#ifndef LATAN_INCLUDES_H_
#define LATAN_INCLUDES_H_

#include "../config.h"

/* replacement functions */
#ifndef HAVE_ACOSH
#ifdef acosh
#undef acosh
#endif
#define acosh gsl_acosh
#endif

/* system includes */
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <unistd.h>
#include <stdarg.h>
#include <limits.h>
#include <float.h>
#include <gsl/gsl_math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* setting names of GPU related functions */
#ifdef HAVE_LIBCUBLAS
#include <latan/latan_mat_cublas.h>
#define mat_gpu_free   mat_cublas_gpu_free
#define mat_on_cpu     mat_cublas_on_cpu
#define mat_on_gpu     mat_cublas_on_gpu
#define mat_sync       mat_cublas_sync
#define mat_gpu_mul_nn mat_cublas_gpu_mul_nn
#define mat_gpu_mul_nt mat_cublas_gpu_mul_nt
#define mat_gpu_mul_tn mat_cublas_gpu_mul_tn
#define mat_gpu_mul_tt mat_cublas_gpu_mul_tt
#else
#define mat_gpu_free(m)
#define mat_on_cpu(m)
#define mat_on_gpu(m)
#define mat_sync(m)
#define mat_gpu_mul_nn(m,n,o) LATAN_FAILURE
#define mat_gpu_mul_nt(m,n,o) LATAN_FAILURE
#define mat_gpu_mul_tn(m,n,o) LATAN_FAILURE
#define mat_gpu_mul_tt(m,n,o) LATAN_FAILURE
#endif

#endif