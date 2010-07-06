#include <stdio.h>
#include <latan/latan_mat.h>
#include <latan/latan_io.h>
#include <latan/latan_math.h>
#include <latan/latan_rand.h>
#include <latan/latan_statistics.h>

#define SEQ_LENGTH 10
#define NDIM 4
#define NBOOT 5

int main(void)
{
	mat *gvec;
	mat mean,var;
	rs_sample s_mean;
	size_t i,j;
	double sigma;
	
	mean	= mat_create(NDIM,1);
	s_mean  = rs_sample_create_boot(NDIM,NBOOT);
	var		= mat_create(NDIM,NDIM);
	gvec	= mat_ar_create_from_dim(SEQ_LENGTH,mean);
	randgen_init_from_time();
	
	printf("-- generating %d gaussian %d-vectors...\n",SEQ_LENGTH,NDIM);
	for (j=0;j<NDIM;j++)
	{
		sigma = DRATIO(j+1,5);
		printf("dimension %d sigma\t: %f\n",j,sigma);
		for (i=0;i<SEQ_LENGTH;i++)
		{
			mat_set(gvec[i],j,0,rand_n(0.0,sigma));
		}
	}
	
	printf("-- computing mean...\n");
	mat_mean(mean,gvec,SEQ_LENGTH);
	printf("mean =\n");
	mat_print(mean);
	printf("\n");
	
	printf("-- computing naive variance...\n");
	mat_var(var,gvec,SEQ_LENGTH);
	printf("variance =\n");
	mat_print(var);
	printf("\n");
	
	printf("-- resampling mean...\n");
	resample(s_mean,gvec,SEQ_LENGTH,1,&rs_mean,NULL);
	
	printf("-- computing variance from resampled sample...\n");
	rs_sample_var(var,s_mean);
	mat_eqmuls(var,SEQ_LENGTH);
	printf("variance =\n");
	mat_print(var);
	printf("\n");
	
	mat_destroy(mean);
	rs_sample_destroy(s_mean);
	mat_destroy(var);
	mat_ar_destroy(gvec,SEQ_LENGTH);
}