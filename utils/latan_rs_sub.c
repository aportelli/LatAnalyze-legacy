#include <stdlib.h>
#include <stdio.h>
#include <latan/latan_statistics.h>
#include <latan/latan_io.h>

#define INF1_NAME argv[1]
#define INF2_NAME argv[2]
#define OUTF_NAME argv[3]

int main(int argc, char* argv[])
{
	rs_sample s1,s2,res;
	size_t s1_nrow,s2_nrow,s1_nsample,s2_nsample;
	size_t i;
	mat sig;
	bool do_save_res;
	
	/* argument parsing */
	switch (argc)
	{
		case 3:
			do_save_res = false;
			break;
		case 4:
			do_save_res = true;
			break;
		default:
			fprintf(stderr,"usage: %s <sample 1> <sample 2> [<output sample>]\n",\
					argv[0]);
			return EXIT_FAILURE;
			break;
	}
	
	/* getting sizes */
	s1_nrow    = rs_sample_load_nrow(INF1_NAME);
	s1_nsample = rs_sample_load_nsample(INF1_NAME);
	s2_nrow    = rs_sample_load_nrow(INF2_NAME);
	s2_nsample = rs_sample_load_nsample(INF2_NAME);
	
	/* error checking */
	if (s1_nsample != s2_nsample)
	{
		fprintf(stderr,"error: number of sample mismatch\n");
		return EXIT_FAILURE;
	}
	if (s1_nrow != s2_nrow)
	{
		fprintf(stderr,"error: sample dimensions mismatch\n");
	}
	
	/* allocation */
	s1 = rs_sample_create(s1_nrow,s1_nsample,"");
	s2 = rs_sample_create(s2_nrow,s2_nsample,"");
	res = rs_sample_create(s1_nrow,s1_nsample,"");
	sig = mat_create(s1_nrow,1);
	
	/* loading samples */
	printf("-- loading resampled sample from %s...\n",INF1_NAME);
	rs_sample_load(s1,INF1_NAME);
	printf("-- loading resampled sample from %s...\n",INF2_NAME);
	rs_sample_load(s2,INF2_NAME);
	
	/* substracting samples */
	printf("-- substracting samples...\n");
	mat_sub(rs_sample_pt_cent_val(res),rs_sample_pt_cent_val(s1),\
			rs_sample_pt_cent_val(s2));
	for (i=0;i<rs_sample_get_nsample(res);i++)
	{
		mat_sub(rs_sample_pt_sample(res,i),rs_sample_pt_sample(s1,i),\
				rs_sample_pt_sample(s2,i));
	}
	
	/* result output */
	rs_sample_varp(sig,res);
	mat_eqsqrt(sig);
	printf("central value:\n");
	mat_print(rs_sample_pt_cent_val(res));
	printf("standard deviation:\n");
	mat_print(sig);
	if (do_save_res)
	{
		rs_sample_save(res,OUTF_NAME);
	}
	
	/* desallocation */
	rs_sample_destroy(s1);
	rs_sample_destroy(s2);
	rs_sample_destroy(res);
	mat_destroy(sig);
	
	return EXIT_SUCCESS;
}