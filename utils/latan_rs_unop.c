#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <latan/latan_statistics.h>
#include <latan/latan_io.h>

#define INF1_NAME argv[1]
#define OUTF_NAME argv[2]

#ifndef UNOP
#error UNOP macro must be defined to compile this program (use -DUNOP=<op> option)
#endif

typedef latan_errno mat_unop_f(mat *, mat *);

mat_unop_f *mat_unop = &UNOP;

int main(int argc, char *argv[])
{
    rs_sample *s1,*res;
    size_t s1_nrow,s1_nsample;
    size_t i;
    mat *sig;
    bool do_save_res;
    strbuf res_name;
    
    /* I/O init */
    io_init();
    
    /* argument parsing */
    switch (argc)
    {
        case 2:
            do_save_res = false;
            strbufcpy(res_name,"");
            break;
        case 3:
            do_save_res = true;
            strbufcpy(res_name,OUTF_NAME);
            break;
        default:
            fprintf(stderr,"usage: %s <sample> [<output sample>]\n",\
                    argv[0]);
            return EXIT_FAILURE;
            break;
    }
    
    /* getting sizes */
    rs_sample_load_nrow(&s1_nrow,INF1_NAME,"");
    rs_sample_load_nsample(&s1_nsample,INF1_NAME,"");
    
    /* allocation */
    s1 = rs_sample_create(s1_nrow,s1_nsample);
    res = rs_sample_create(s1_nrow,s1_nsample);
    sig = mat_create(s1_nrow,1);
    
    /* loading samples */
    printf("-- loading resampled sample from %s...\n",INF1_NAME);
    rs_sample_load(s1,INF1_NAME,"");
    
    /* multiplying samples */
    printf("-- executing operation on sample...\n");
    mat_unop(rs_sample_pt_cent_val(res),rs_sample_pt_cent_val(s1));
    for (i=0;i<rs_sample_get_nsample(res);i++)
    {
        mat_unop(rs_sample_pt_sample(res,i),rs_sample_pt_sample(s1,i));
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
    rs_sample_destroy(res);
    mat_destroy(sig);

    /* I/O finish */
    io_finish();
    
    return EXIT_SUCCESS;
}