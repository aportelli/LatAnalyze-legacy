#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <latan/latan_statistics.h>
#include <latan/latan_io.h>

#ifndef UNOP
#error UNOP macro must be defined to compile this program (use -DUNOP=<op> option)
#endif

int main(int argc, char *argv[])
{
    rs_sample *s1,*res;
    size_t s1_dim[2],s1_nsample;
    int i,j;
    mat *sig;
    bool do_save_res, show_usage;
    char *spt;
    strbuf out_elname,out_fname,in_elname,in_fname,in_path,out_path;
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
                if (j == 0)
                {
                    strbufcpy(in_path,argv[i]);
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
        fprintf(stderr,"usage: %s <in sample> [-o <out sample>] [-f {ascii|xml}]\n",\
                argv[0]);
        return EXIT_FAILURE;
    }
    
    /* I/O init */
    io_set_fmt(fmt);
    io_init();
    
    /* getting sizes */
    rs_sample_load(NULL,&s1_nsample,s1_dim,in_path);
    
    /* allocation */
    s1  = rs_sample_create(s1_dim[0],s1_dim[1],s1_nsample);
    res = rs_sample_create(s1_dim[0],s1_dim[1],s1_nsample);
    sig = mat_create(s1_dim[0],s1_dim[1]);
    
    /* loading samples */
    printf("-- loading sample from %s...\n",in_path);
    rs_sample_load(s1,NULL,NULL,in_path);
    
    /* multiplying samples */
    printf("-- executing operation on sample...\n");
    rs_sample_unop(res,s1,&UNOP);
    
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
            get_elname(in_fname,in_elname,in_path);
            spt = (strlen(in_elname) == 0) ? in_fname : in_elname;
            sprintf(out_path,"%s%c%s_%s",out_fname,LATAN_PATH_SEP,argv[0],spt);
        }
        rs_sample_save(out_path,'w',res);
    }
    
    /* desallocation */
    rs_sample_destroy(s1);
    rs_sample_destroy(res);
    mat_destroy(sig);

    /* I/O finish */
    io_finish();
    
    return EXIT_SUCCESS;
}
