/* latan_mat_cublas.h, part of LatAnalyze library
 *
 * Copyright (C) 2010, 2011 Antonin Portelli
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef LATAN_MAT_CUBLAS
#define LATAN_MAT_CUBLAS

#include <latan/latan_globals.h>

#ifdef HAVE_LIBCUBLAS
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

#endif
