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
	
	status = LATAN_SUCCESS;
	mat_zero(mean);
	
	for (i=0;i<size;i++)
	{
		LATAN_UPDATE_STATUS(status,mat_eqadd(mean,m[i]));
	}
	LATAN_UPDATE_STATUS(status,mat_eqmuls(mean,1.0/((double)size)));
	
	return status;
}

int mat_meansig(mat mean, mat sig, const mat *m, const size_t size)
{
	int status;
	size_t i;
	const double dsize = ((double)(size));
	mat sq;
	
	if (!mat_issamedim(mean,sig))
	{
		LATAN_ERROR("operation between matrices with dimension mismatch",\
					LATAN_EBADLEN);
	}
	
	status = LATAN_SUCCESS;
	mat_zero(mean);
	mat_zero(sig);
	
	mat_create(&sq,nrow(mean),ncol(mean));
	
	for (i=0;i<size;i++)
	{
		LATAN_UPDATE_STATUS(status,mat_eqadd(mean,m[i]));
		LATAN_UPDATE_STATUS(status,mat_mulp(sq,m[i],m[i]));
		LATAN_UPDATE_STATUS(status,mat_eqadd(sig,sq));
	}
	LATAN_UPDATE_STATUS(status,mat_eqmuls(mean,1.0/dsize));
	LATAN_UPDATE_STATUS(status,mat_eqmuls(sig,1.0/(dsize-1.0)));
	LATAN_UPDATE_STATUS(status,mat_mulp(sq,mean,mean));
	LATAN_UPDATE_STATUS(status,mat_eqmuls(sq,(dsize)/(dsize-1.0)));
	LATAN_UPDATE_STATUS(status,mat_eqsub(sig,sq));
	LATAN_UPDATE_STATUS(status,mat_eqsqrt(sig));
	
	mat_destroy(&sq);
	
	return status;
}
