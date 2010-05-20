#include <stdio.h>
#include <latan/latan_mat.h>
#include <latan/latan_io.h>
#include <latan/latan_math.h>
#include <latan/latan_rand.h>
#include <latan/latan_statistics.h>

#define SEQ_LENGTH 1000000
#define NDIM 4

int main(void)
{
	mat* gvec;
	mat mean,var;
	size_t i,j;
	double sigma;
	
	mean	= mat_create(NDIM,1);
	var		= mat_create(NDIM,NDIM);
	gvec	= mat_ar_create_from_dim(SEQ_LENGTH,mean);
	randgen_init_from_time();
	
	for (i=0;i<SEQ_LENGTH;i++)
	{
		for (j=0;j<NDIM;j++)
		{
			sigma = DRATIO(j+1,5);
			mat_set(gvec[i],j,0,rand_n(0.0,sigma));
		}
	}
	mat_mean(mean,gvec,SEQ_LENGTH);
	mat_print(mean);
	printf("\n");
	mat_var(var,gvec,SEQ_LENGTH);
	mat_print(var);
}