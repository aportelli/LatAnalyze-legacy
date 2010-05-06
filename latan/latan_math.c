#include <latan/latan_math.h>
#include <latan/latan_includes.h>

unsigned int binomial(const unsigned int n, const unsigned int p)
{
	unsigned int* b;
	unsigned int i,j;
	unsigned int res;
	
	MALLOC_NOERRET(b,unsigned int*,n+1);
	
	b[0]=1;
	for (i=1;i<=n;i++)
	{
		b[i] = 1;
		for (j=i-1;j>0;j--)
		{
			b[j] += b[j-1];
		}
	}
	res = b[p];
	
	FREE(b);
	
	return res;
}

latan_errno finite_diff(mat ddat, const mat dat)
{
	size_t i,j;
	double ddat_ij;
	
	if (nrow(ddat) != nrow(dat) - 2)
	{
		LATAN_ERROR("derivative matrix have wrong dimensions",LATAN_EBADLEN);
	}
	
	FOR_VAL(ddat,i,j)
	{
		ddat_ij = 0.5*(mat_get(dat,i+2,j) - mat_get(dat,i,j));
		mat_set(ddat,i,j,ddat_ij);
	}
	
	return LATAN_SUCCESS;
}