#ifndef LATAN_MAT_CUBLAS
#define LATAN_MAT_CUBLAS

#include <latan/latan_globals.h>

__BEGIN_DECLS

/* allocation */
latan_errno mat_cublas_gpu_alloc(mat *m);
latan_errno mat_cublas_gpu_free(mat *m);

/* CPU <-> GPU management */
latan_errno mat_cublas_cp_cpu_to_gpu(mat *m);
latan_errno mat_cublas_cp_gpu_to_cpu(mat *m);
latan_errno mat_cublas_on_cpu(mat *m);
latan_errno mat_cublas_on_gpu(mat *m);
latan_errno mat_cublas_sync(mat *m);

/* operations */
latan_errno mat_cublas_gpu_mul_nn(mat *m, mat *n, mat *o);
latan_errno mat_cublas_gpu_mul_nt(mat *m, mat *n, mat *o);
latan_errno mat_cublas_gpu_mul_tn(mat *m, mat *n, mat *o);
latan_errno mat_cublas_gpu_mul_tt(mat *m, mat *n, mat *o);

__END_DECLS

#endif
