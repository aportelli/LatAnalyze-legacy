#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <latan/latan_statistics.h>
#include <latan/latan_io.h>

#ifndef BINOP
#error BINOP macro must be defined to compile this program (use -DBINOP=<op> option)
#endif

int main(int argc, char *argv[])
{
    int i,j;
    rs_sample *s1,*s2,*res;
    size_t s1_dim[2],s2_dim[2],s1_nsample,s2_nsample;
    mat *sig;
    bool do_save_res, show_usage;
    strbuf out_elname,out_fname,in_elname[2],in_fname[2],in_path[2],out_path;
    char *spt[2];
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
                    strbufcpy(out_path,argv[i+1]);
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
                    strbufcpy(in_path[j],argv[i]);
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
    rs_sample_load(NULL,&s1_nsample,s1_dim,in_path[0]);
    rs_sample_load(NULL,&s2_nsample,s2_dim,in_path[1]);
    
    /* allocation */
    s1  = rs_sample_create(s1_dim[0],s1_dim[1],s1_nsample);
    s2  = rs_sample_create(s2_dim[0],s2_dim[1],s2_nsample);
    res = rs_sample_create(s1_dim[0],s1_dim[1],s1_nsample);
    sig = mat_create(s1_dim[0],s1_dim[1]);
    
    /* loading samples */
    printf("-- loading sample from %s...\n",in_path[0]);
    rs_sample_load(s1,NULL,NULL,in_path[0]);
    printf("-- loading sample from %s...\n",in_path[1]);
    rs_sample_load(s2,NULL,NULL,in_path[1]);
    
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
        get_elname(out_fname,out_elname,out_path);
        if (strlen(out_elname) == 0)
        {
            get_elname(in_fname[0],in_elname[0],in_path[0]);
            get_elname(in_fname[1],in_elname[1],in_path[1]);
            spt[0] = (strlen(in_elname[0]) == 0) ? in_fname[0] : in_elname[0];
            spt[1] = (strlen(in_elname[1]) == 0) ? in_fname[1] : in_elname[1];
            sprintf(out_path,"%s%c%s_%s_%s",out_fname,LATAN_PATH_SEP,argv[0],\
                    spt[0],spt[1]);
        }
        rs_sample_save(out_path,'w',res);
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
