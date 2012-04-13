#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <latan/latan_statistics.h>
#include <latan/latan_io.h>

#ifndef BINOP
#error BINOP macro must be defined to compile this program (use -DBINOP=<op> option)
#endif

#define STRINGIFY(x) #x

int main(int argc, char *argv[])
{
    int i,j;
    rs_sample *s1,*s2,*res;
    size_t s1_dim[2],s2_dim[2],s1_nsample,s2_nsample;
    mat *sig;
    bool do_save_res, show_usage;
    strbuf res_name,inf_name[2],outf_name;
    io_fmt_no fmt;
    
    /* argument parsing */
    i           = 1;
    j           = 0;
    show_usage  = false;
    do_save_res = false;
    fmt         = io_get_fmt();
    
    if (argc <= 1)
    {
        show_usage = true;
    }
    else
    {
        while (i < argc)
        {
            if (strcmp(argv[i],"-o") == 0)
            {
                if (i == argc - 1)
                {
                    show_usage = true;
                    break;
                }
                else
                {
                    strbufcpy(outf_name,argv[i+1]);
                    do_save_res = true;
                    i += 2;
                }
            }
            else if (strcmp(argv[i],"-f") == 0)
            {
                if (i == argc - 1)
                {
                    show_usage = true;
                    break;
                }
                else
                {
                    if (strcmp(argv[i+1],"xml") == 0)
                    {
                        fmt = IO_XML;
                    }
                    else if (strcmp(argv[i+1],"ascii") == 0)
                    {
                        fmt = IO_ASCII;
                    }
                    else
                    {
                        fprintf(stderr,"error: format %s unknown\n",argv[i+1]);
                        return EXIT_FAILURE;
                    }
                    i += 2;
                }
            }
            else
            {
                if (j < 2)
                {
                    strbufcpy(inf_name[j],argv[i]);
                    j++;
                    i++;
                }
                else
                {
                    show_usage = true;
                    break;
                }
            }
        }
    }
    if (show_usage)
    {
        fprintf(stderr,"usage: %s <in sample 1> <in sample 2> [-o <out sample>] [-f {ascii|xml}]\n",\
                argv[0]);
        return EXIT_FAILURE;
    }
    
    /* I/O init */
    io_set_fmt(fmt);
    io_init();
    
    /* getting sizes */
    rs_sample_load(NULL,&s1_nsample,s1_dim,inf_name[0]);
    rs_sample_load(NULL,&s2_nsample,s2_dim,inf_name[1]);
    
    /* allocation */
    s1  = rs_sample_create(s1_dim[0],s1_dim[1],s1_nsample);
    s2  = rs_sample_create(s2_dim[0],s2_dim[1],s2_nsample);
    res = rs_sample_create(s1_dim[0],s1_dim[1],s1_nsample);
    sig = mat_create(s1_dim[0],s1_dim[1]);
    
    /* loading samples */
    printf("-- loading sample from %s...\n",inf_name[0]);
    rs_sample_load(s1,NULL,NULL,inf_name[0]);
    printf("-- loading sample from %s...\n",inf_name[1]);
    rs_sample_load(s2,NULL,NULL,inf_name[1]);
    
    /* multiplying samples */
    printf("-- executing operation on samples...\n");
    rs_sample_binop(res,s1,s2,&BINOP);
    
    /* result output */
    rs_sample_varp(sig,res);
    mat_eqsqrt(sig);
    printf("central value:\n");
    mat_print(rs_sample_pt_cent_val(res),"%e");
    printf("standard deviation:\n");
    mat_print(sig,"%e");
    if (do_save_res)
    {
        sprintf(res_name,"%s%c%s_%s_%s",outf_name,LATAN_PATH_SEP,"",\
                STRINGIFY(BINOP),"");
        rs_sample_save(res_name,'w',res);
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
