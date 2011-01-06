/* latan_mat_cublas.c, part of LatAnalyze library
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

#include <latan/latan_mat_cublas.h>
#ifdef HAVE_LIBCUBLAS
#include <cublas.h>

#define INT(X) ((int)(X))

static latan_errno mat_cublas_error(cublasStatus stat);

/*                              internal                                    */
/****************************************************************************/
static latan_errno mat_cublas_error(cublasStatus stat)
{
    switch (stat)
    {
        case CUBLAS_STATUS_SUCCESS:
            break;
        case CUBLAS_STATUS_NOT_INITIALIZED:
            LATAN_ERROR("CUBLAS is not initialized",LATAN_ESYSTEM);
            break;
        case CUBLAS_STATUS_ALLOC_FAILED:
            LATAN_ERROR("(CUBLAS) allocation failed",LATAN_ENOMEM);
            break;
        case CUBLAS_STATUS_INVALID_VALUE:
            LATAN_ERROR("(CUBLAS) invalid value",LATAN_EINVAL);
            break;
        case CUBLAS_STATUS_ARCH_MISMATCH:
            LATAN_ERROR("(CUBLAS) operation not supported by architecture",\
                        LATAN_ESYSTEM);
            break;
        case CUBLAS_STATUS_MAPPING_ERROR:
            LATAN_ERROR("(CUBLAS) mapping error",LATAN_ESYSTEM);
            break;
        case CUBLAS_STATUS_EXECUTION_FAILED:
            LATAN_ERROR("(CUBLAS) execution failed",LATAN_ESYSTEM);
            break;
        case CUBLAS_STATUS_INTERNAL_ERROR:
            LATAN_ERROR("CUBLAS internal error",LATAN_FAILURE);
            break;
        default:
            LATAN_ERROR("invalid CUBLAS error code",LATAN_EINVAL);
            break;
    }
    
    return LATAN_SUCCESS;
}

/*                              allocation                                  */
/****************************************************************************/
latan_errno mat_cublas_gpu_alloc(mat *m)
{
    cublasStatus status;
    
    status = cublasAlloc(INT(nrow(m)*ncol(m)),sizeof(double),\
                         (void**)&(m->data_gpu));
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        LATAN_ERROR("GPU memory allocation failed",LATAN_ENOMEM);
    }
    m->mem_flag |= GPU_ALLOCATED;
    
    return LATAN_SUCCESS;
}

latan_errno mat_cublas_gpu_free(mat *m)
{
    cublasStatus status;
    
    status = cublasFree((void*)(m->data_gpu));
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        LATAN_ERROR("GPU memory desallocation failed",LATAN_ENOMEM);
    }
    m->mem_flag -= m->mem_flag&GPU_ALLOCATED;
    
    return LATAN_SUCCESS;
}

/*                          CPU<->GPU management                            */
/****************************************************************************/
latan_errno mat_cublas_cp_cpu_to_gpu(mat *m)
{
    cublasStatus status;
    size_t row;

    for (row=0;row<nrow(m);row++)
    {
        status = cublasSetVector(INT(ncol(m)),sizeof(double),         \
                                 m->data_cpu->data+INT(row*ncol(m)),1,\
                                 m->data_gpu+INT(row),INT(nrow(m)));
        if (status != CUBLAS_STATUS_SUCCESS)
        {
            LATAN_ERROR("CPU->GPU memory copy failed",LATAN_ESYSTEM);
        }
    }
    
    return LATAN_SUCCESS;
}

latan_errno mat_cublas_cp_gpu_to_cpu(mat *m)
{
    cublasStatus status;
    size_t row;

    for (row=0;row<nrow(m);row++)
    {
        status = cublasGetVector(INT(ncol(m)),sizeof(double),        \
                                 m->data_gpu+INT(row),INT(nrow(m)),  \
                                 m->data_cpu->data+INT(row*ncol(m)),1);
        if (status != CUBLAS_STATUS_SUCCESS)
        {
            LATAN_ERROR("GPU->CPU memory copy failed",LATAN_ESYSTEM);
        }
    }
    
    return LATAN_SUCCESS;
}

latan_errno mat_cublas_on_cpu(mat *m)
{
    latan_errno status;
    
    status = LATAN_SUCCESS;
    if (m->mem_flag & GPU_LAST)
    {
        status = mat_cublas_cp_gpu_to_cpu(m);
        MAT_SYNCED(m);
    }
    
    return status;
}

latan_errno mat_cublas_on_gpu(mat *m)
{
    latan_errno status;
    
    status = LATAN_SUCCESS;
    if (!(m->mem_flag & GPU_ALLOCATED))
    {
        LATAN_UPDATE_STATUS(status,mat_cublas_gpu_alloc(m));
    }
    if (m->mem_flag & CPU_LAST)
    {
        LATAN_UPDATE_STATUS(status,mat_cublas_cp_cpu_to_gpu(m));
        MAT_SYNCED(m);
    }
    
    return status;
}

latan_errno mat_cublas_sync(mat *m)
{
    latan_errno status;
    
    status = LATAN_SUCCESS;
    if (m->mem_flag & CPU_LAST)
    {
        status = mat_cublas_on_gpu(m);
    }
    else if (m->mem_flag & GPU_LAST)
    {
        status = mat_cublas_on_cpu(m);
    }
    
    return status;
}

latan_errno mat_cublas_gpu_mul_nn(mat *m, mat *n, mat *o)
{       
    cublasDgemm('n','n',INT(nrow(m)),INT(ncol(m)),INT(ncol(n)),1.0,n->data_gpu,\
                INT(nrow(n)),o->data_gpu,INT(nrow(o)),0.0,m->data_gpu,         \
                INT(nrow(m)));
    mat_cublas_error(cublasGetError());
    
    return EXIT_SUCCESS;
}

latan_errno mat_cublas_gpu_mul_nt(mat *m, mat *n, mat *o)
{
    cublasDgemm('n','t',INT(nrow(m)),INT(ncol(m)),INT(ncol(n)),1.0,n->data_gpu,\
                INT(nrow(n)),o->data_gpu,INT(nrow(o)),0.0,m->data_gpu,         \
                INT(nrow(m)));
    mat_cublas_error(cublasGetError());
    
    return EXIT_SUCCESS;
}

latan_errno mat_cublas_gpu_mul_tn(mat *m, mat *n, mat *o)
{
    cublasDgemm('t','n',INT(nrow(m)),INT(ncol(m)),INT(nrow(n)),1.0,n->data_gpu,\
                INT(nrow(n)),o->data_gpu,INT(nrow(o)),0.0,m->data_gpu,         \
                INT(nrow(m)));
    mat_cublas_error(cublasGetError());

    return EXIT_SUCCESS;
}

latan_errno mat_cublas_gpu_mul_tt(mat *m, mat *n, mat *o)
{
    cublasDgemm('t','t',INT(nrow(m)),INT(ncol(m)),INT(nrow(n)),1.0,n->data_gpu,\
                INT(nrow(n)),o->data_gpu,INT(nrow(o)),0.0,m->data_gpu,         \
                INT(nrow(m)));
    mat_cublas_error(cublasGetError());
    
    return EXIT_SUCCESS;
}

#endif