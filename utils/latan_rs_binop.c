#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <latan/latan_statistics.h>
#include <latan/latan_io.h>

#define INF1_NAME argv[1]
#define INF2_NAME argv[2]
#define OUTF_NAME argv[3]

#ifndef BINOP
#error BINOP macro must be defined to compile this program (use -DBINOP=<op> option)
#endif

typedef latan_errno mat_binop_f(mat *, const mat *, const mat *);

mat_binop_f *mat_binop = &BINOP;

int main(int argc, char *argv[])
{
    rs_sample *s1,*s2,*res;
    size_t s1_nrow,s2_nrow,s1_nsample,s2_nsample;
    size_t i;
    mat *sig;
    bool do_save_res;
    strbuf res_name;

    /* I/O init */
    io_init();
    
    /* argument parsing */
    switch (argc)
    {
        case 3:
            do_save_res = false;
            strbufcpy(res_name,"");
            break;
        case 4:
            do_save_res = true;
            strbufcpy(res_name,OUTF_NAME);
            break;
        default:
            fprintf(stderr,"usage: %s <sample 1> <sample 2> [<output sample>]\n",\
                    argv[0]);
            return EXIT_FAILURE;
            break;
    }
    
    /* getting sizes */
    rs_sample_load_nrow(&s1_nrow,INF1_NAME,"");
    rs_sample_load_nsample(&s1_nsample,INF1_NAME,"");
    rs_sample_load_nrow(&s2_nrow,INF2_NAME,"");
    rs_sample_load_nsample(&s2_nsample,INF2_NAME,"");
    
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
    s1 = rs_sample_create(s1_nrow,s1_nsample);
    s2 = rs_sample_create(s2_nrow,s2_nsample);
    res = rs_sample_create(s1_nrow,s1_nsample);
    sig = mat_create(s1_nrow,1);
    
    /* loading samples */
    printf("-- loading resampled sample from %s...\n",INF1_NAME);
    rs_sample_load(s1,INF1_NAME,"");
    printf("-- loading resampled sample from %s...\n",INF2_NAME);
    rs_sample_load(s2,INF2_NAME,"");
    
    /* multiplying samples */
    printf("-- executing operation on samples...\n");
    mat_binop(rs_sample_pt_cent_val(res),rs_sample_pt_cent_val(s1),\
              rs_sample_pt_cent_val(s2));
    for (i=0;i<rs_sample_get_nsample(res);i++)
    {
        mat_binop(rs_sample_pt_sample(res,i),rs_sample_pt_sample(s1,i),\
                  rs_sample_pt_sample(s2,i));
    }
    
    /* result output */
    rs_sample_varp(sig,res);
    mat_eqsqrt(sig);
    printf("central value:\n");
    mat_print(rs_sample_pt_cent_val(res),"%e");
    printf("standard deviation:\n");
    mat_print(sig,"%e");
    if (do_save_res)
    {
        rs_sample_set_name(res,res_name);
        rs_sample_save(OUTF_NAME,'w',res);
    }
    
    /* desallocation */
    rs_sample_destroy(s1);
    rs_sample_destroy(s2);
    rs_sample_destroy(res);
    mat_destroy(sig);

    /* I/O finish */
    io_finish();
    
    return EXIT_SUCCESS;
}