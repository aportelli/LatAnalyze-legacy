#include <latan/latan_mat_cublas.h>
#include <cublas.h>

#define INT(X) ((int)(X))

/*								allocation									*/
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

/*							CPU<->GPU management							*/
/****************************************************************************/
latan_errno mat_cublas_cp_cpu_to_gpu(mat *m)
{
	cublasStatus status;
	size_t row;
	
	for (row=0;row<nrow(m);row++)
	{
		status = cublasSetVector(INT(ncol(m)),sizeof(double),m->data_cpu,\
								 1,m->data_gpu,INT(ncol(m)));
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
		status = cublasGetVector(INT(ncol(m)),sizeof(double),m->data_gpu,\
								 INT(ncol(m)),m->data_cpu,1);
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
	mat_cublas_on_gpu(m);
	mat_cublas_on_gpu(n);
	mat_cublas_on_gpu(o);
	
	if ((nrow(m) != nrow(n))||(ncol(m) != ncol(o))||(ncol(n) != nrow(o)))
	{
		LATAN_ERROR("operation between matrices with dimension mismatch",\
					LATAN_EBADLEN);
	}
	
	cublasDgemm('n','n',INT(nrow(n)),INT(ncol(o)),INT(ncol(n)),1.0,n->data_gpu,\
				INT(nrow(n)),o->data_gpu,INT(nrow(o)),0.0,m->data_gpu,         \
				INT(nrow(m)));
	MAT_GPU_LAST(m);
	
	return EXIT_SUCCESS;
}

latan_errno mat_cublas_gpu_mul_nt(mat *m, mat *n, mat *o)
{
	mat_cublas_on_gpu(m);
	mat_cublas_on_gpu(n);
	mat_cublas_on_gpu(o);
	
	if ((nrow(m) != nrow(n))||(ncol(m) != nrow(o))||(ncol(n) != ncol(o)))
	{
		LATAN_ERROR("operation between matrices with dimension mismatch",\
					LATAN_EBADLEN);
	}
	
	cublasDgemm('n','t',INT(nrow(n)),INT(ncol(o)),INT(ncol(n)),1.0,n->data_gpu,\
				INT(nrow(n)),o->data_gpu,INT(nrow(o)),0.0,m->data_gpu,         \
				INT(nrow(m)));
	MAT_GPU_LAST(m);
	
	return EXIT_SUCCESS;
}

latan_errno mat_cublas_gpu_mul_tn(mat *m, mat *n, mat *o)
{
	mat_cublas_on_gpu(m);
	mat_cublas_on_gpu(n);
	mat_cublas_on_gpu(o);
	
	if ((nrow(m) != ncol(n))||(ncol(m) != ncol(o))||(nrow(n) != nrow(o)))
	{
		LATAN_ERROR("operation between matrices with dimension mismatch",\
					LATAN_EBADLEN);
	}
	
	cublasDgemm('t','n',INT(nrow(n)),INT(ncol(o)),INT(ncol(n)),1.0,n->data_gpu,\
				INT(nrow(n)),o->data_gpu,INT(nrow(o)),0.0,m->data_gpu,         \
				INT(nrow(m)));
	MAT_GPU_LAST(m);
	
	return EXIT_SUCCESS;
}

latan_errno mat_cublas_gpu_mul_tt(mat *m, mat *n, mat *o)
{
	mat_cublas_on_gpu(m);
	mat_cublas_on_gpu(n);
	mat_cublas_on_gpu(o);
	
	if ((nrow(m) != ncol(n))||(ncol(m) != nrow(o))||(nrow(n) != ncol(o)))
	{
		LATAN_ERROR("operation between matrices with dimension mismatch",\
					LATAN_EBADLEN);
	}
	
	cublasDgemm('t','t',INT(nrow(n)),INT(ncol(o)),INT(ncol(n)),1.0,n->data_gpu,\
				INT(nrow(n)),o->data_gpu,INT(nrow(o)),0.0,m->data_gpu,         \
				INT(nrow(m)));
	MAT_GPU_LAST(m);
	
	return EXIT_SUCCESS;
}