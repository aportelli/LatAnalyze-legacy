/* latan_io_ascii.c, part of LatAnalyze library
 *
 * Copyright (C) 2010, 2011 Antonin Portelli
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <latan/latan_io_ascii.h>
#include <latan/latan_includes.h>
#include <latan/latan_io.h>

/*                       file buffer management (internal)                  */
/****************************************************************************/
typedef struct
{
    FILE **f_buf;
    strbuf *fname;
    bool *file_is_loaded;
    int nfile;
} io_ascii_env;

static io_ascii_env env =
{
    NULL,\
    NULL,\
    NULL,\
    0    \
};

static latan_errno ascii_new_file_buf(const strbuf fname)
{
    latan_errno status;
    int nthread,thread,i;

#ifdef _OPENMP
    nthread = omp_get_num_threads();
    thread  = omp_get_thread_num();
#else
    nthread = 1;
    thread  = 0;
#endif

#ifdef _OPENMP
    #pragma omp critical
#endif
    {
        if (nthread > env.nfile)
        {
            REALLOC_NOERRET(env.f_buf,env.f_buf,FILE **,nthread);
            REALLOC_NOERRET(env.fname,env.fname,strbuf *,nthread);
            REALLOC_NOERRET(env.file_is_loaded,env.file_is_loaded,bool *,\
                            nthread);
            for (i=env.nfile;i<nthread;i++)
            {
                env.file_is_loaded[i] = false;
                strbufcpy(env.fname[i],"");
                env.f_buf[i]        = NULL;
            }
            env.nfile = nthread;
        }
    }
    if (env.file_is_loaded[thread])
    {
        fclose(env.f_buf[thread]);
    }
    FOPEN(env.f_buf[thread],fname,"w");
    strbufcpy(env.fname[thread],fname);
    env.file_is_loaded[thread] = true;

    return status;
}

static latan_errno ascii_open_file_buf(const strbuf fname, const char mode)
{
    latan_errno status;
    strbuf smode;
    int nthread,thread,i;

#ifdef _OPENMP
    nthread = omp_get_num_threads();
    thread  = omp_get_thread_num();
#else
    nthread = 1;
    thread  = 0;
#endif
    status   = LATAN_SUCCESS;
    smode[0] = mode;
    smode[1] = '\0';

#ifdef _OPENMP
    #pragma omp critical
#endif
    {
        if (nthread > env.nfile)
        {
            REALLOC_NOERRET(env.f_buf,env.f_buf,FILE **,nthread);
            REALLOC_NOERRET(env.fname,env.fname,strbuf *,nthread);
            REALLOC_NOERRET(env.file_is_loaded,env.file_is_loaded,bool *,\
                            nthread);
            for (i=env.nfile;i<nthread;i++)
            {
                env.file_is_loaded[i] = false;
                env.f_buf[i]        = NULL;
            }
            env.nfile = nthread;
        }
    }
    if (env.file_is_loaded[thread])
    {
        if (strcmp(env.fname[thread],fname) != 0)
        {
            fclose(env.f_buf[thread]);
            FOPEN(env.f_buf[thread],fname,smode);
            strbufcpy(env.fname[thread],fname);
        }
    }
    else
    {
        FOPEN(env.f_buf[thread],fname,smode);
        strbufcpy(env.fname[thread],fname);
        env.file_is_loaded[thread] = true;
    }

    return status;
}

/*                              I/O init/finish                             */
/****************************************************************************/
static bool io_is_init = false;

void io_init_ascii(void)
{
    if (!io_is_init)
    {
#ifdef _OPENMP
        if(omp_in_parallel())
        {
            LATAN_WARNING("I/O initialization called from a parallel region",\
                          LATAN_FAILURE);
        }
#endif
        io_is_init = true;
    }
}

void io_finish_ascii(void)
{
    int i;

    i = 0;

    if (io_is_init)
    {
#ifdef _OPENMP
        if(omp_in_parallel())
        {
            LATAN_WARNING("I/O finish called from a parallel region",\
                          LATAN_FAILURE);
        }
#endif
        for (i=0;i<env.nfile;i++)
        {
            if (env.file_is_loaded[i])
            {
                fclose(env.f_buf[i]);
                strbufcpy(env.fname[i],"");
                env.file_is_loaded[i] = false;
            }
        }
        io_is_init = false;
    }
}

/*                          propagator I/O                                  */
/****************************************************************************/
latan_errno prop_load_nt_ascii(size_t *nt, const channel_no channel,\
                               const quark_no q1, const quark_no q2,\
                               const ss_no source, const ss_no sink,\
                               strbuf fname)
{
    strbuf *field,channel_id,q1_id,q2_id,source_id,sink_id;
    bool is_in_prop;
    int nf,lc,buf,thread;

#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif
    field      = NULL;
    is_in_prop = false;

    ascii_open_file_buf(fname,'r');
    channel_id_get(channel_id,channel);
    quark_id_get(q1_id,q1);
    quark_id_get(q2_id,q2);
    ss_id_get(source_id,source);
    ss_id_get(sink_id,sink);
    BEGIN_FOR_LINE_TOK_F(field,env.f_buf[thread]," ",nf,lc)
    {
        if((nf >= 6)&&(strcmp(field[0],"#")==0)&&                          \
           (strcmp(field[1],source_id)==0)&&(strcmp(field[2],sink_id)==0)&&\
           (strcmp(field[3],channel_id)==0)&&                              \
           (strcmp(field[4],q1_id)==0)&&(strcmp(field[5],q2_id)==0))
        {
            is_in_prop = true;
        }
        else if (is_in_prop)
        {
            if (sscanf(field[0],"%d",&buf) <= 0)
            {
                strbuf errmsg;
                sprintf(errmsg,"unable to read time extent (%s:%d)",fname,lc);
                LATAN_ERROR(errmsg,LATAN_ELATSYN);
            }
            *nt = (size_t)buf;
            break;
        }
    }
    END_FOR_LINE_TOK_F(field);

    return LATAN_SUCCESS;
}

latan_errno prop_load_ascii(mat *prop, const channel_no channel, \
                            const quark_no q1, const quark_no q2,\
                            const ss_no source, const ss_no sink,\
                            strbuf fname)
{
    strbuf *field,channel_id,q1_id,q2_id,source_id,sink_id;
    bool is_in_prop,is_first_line_in_prop;
    int nf,lc,thread;
    int i;
    double buf;
    size_t t;

#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread  = 0;
#endif
    field                 = NULL;
    is_in_prop            = false;
    is_first_line_in_prop = false;
    i                     = 0;
    t                     = 0;

    ascii_open_file_buf(fname,'r');
    channel_id_get(channel_id,channel);
    quark_id_get(q1_id,q1);
    quark_id_get(q2_id,q2);
    ss_id_get(source_id,source);
    ss_id_get(sink_id,sink);
    BEGIN_FOR_LINE_TOK_F(field,env.f_buf[thread]," ",nf,lc)
    {
        if(is_in_prop&&(strcmp(field[0],"#") == 0))
        {
            is_in_prop = false;
            break;
        }
        else if((nf >= 6)&&(strcmp(field[0],"#")==0)&&                         \
               (strcmp(field[1],source_id)==0)&&(strcmp(field[2],sink_id)==0)&&\
               (strcmp(field[3],channel_id)==0)&&                              \
               (strcmp(field[4],q1_id)==0)&&(strcmp(field[5],q2_id)==0))
        {
            is_in_prop            = true;
            is_first_line_in_prop = true;
            continue;
        }
        else if (is_in_prop)
        {
            if (is_first_line_in_prop)
            {
                is_first_line_in_prop = false;
                continue;
            }
            else
            {
                for (i=0;i<nf;i++)
                {
                    if (sscanf(field[i],"%lf",&buf) > 0)
                    {
                        mat_set(prop,t,0,buf);
                        t++;
                        continue;
                    }
                    else
                    {
                        strbuf errmsg;
                        sprintf(errmsg,"error while parsing propagator (%s:%d)",\
                                fname,lc);
                        LATAN_ERROR(errmsg,LATAN_ELATSYN);
                    }
                }
            }
        }
    }
    END_FOR_LINE_TOK_F(field);

    return LATAN_SUCCESS;
}

latan_errno prop_save_ascii(strbuf fname, const char mode, mat *prop,\
                            const strbuf channel,                    \
                            const quark_no q1, const quark_no q2,    \
                            const ss_no source, const ss_no sink,    \
                            const strbuf name)
{
    strbuf smode,source_id,sink_id;
    int thread;
    
#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif
    name     = NULL;
    smode[0] = mode;
    smode[1] = '\0';

    if (mode == 'w')
    {
        ascii_new_file_buf(fname);
    }
    else
    {
        ascii_open_file_buf(fname,mode);
    }
    ss_id_get(source_id,source);
    ss_id_get(sink_id,sink);
    fprintf(env.f_buf[thread],"# %s %s %-8.8s %d %d\n",source_id,sink_id,\
            channel,q1,q2);
    fprintf(env.f_buf[thread],"%d\n",(int)nrow(prop));
    mat_dump(env.f_buf[thread],prop,"%.15e");
    fprintf(env.f_buf[thread],"\n");

    return LATAN_SUCCESS;
}

/*                      random generator state I/O                          */
/****************************************************************************/
latan_errno randgen_save_state_ascii(const strbuf fname, const char mode,   \
                                     const rg_state state, const strbuf name)
{
    int thread;
    int i;
    strbuf smode;
    
#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif
    smode[0] = mode;
    smode[1] = '\0';

    if (mode == 'w')
    {
        ascii_new_file_buf(fname);
    }
    else
    {
        LATAN_ERROR("only 'w' file mode is authorized for saving random generator state in ASCII format",\
                    LATAN_EINVAL);
    }
    fprintf(env.f_buf[thread],"# latan_randgen_state %s\n",name);
    for (i=0;i<RLXG_STATE_SIZE;i++)
    {
        fprintf(env.f_buf[thread],"%d\n",state[i]);
    }
    fprintf(env.f_buf[thread],"\n");
    
    return LATAN_SUCCESS;
}

latan_errno randgen_load_state_ascii(rg_state state, const strbuf fname,\
                                     const strbuf name)
{
    int thread,nf,lc;
    int i,j;
    strbuf *field;

#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif
    field    = NULL;
    name     = NULL;
    i        = 0;
    j        = 0;

    ascii_open_file_buf(fname,'r');
    BEGIN_FOR_LINE_TOK_F(field,env.f_buf[thread]," ",nf,lc)
    {
        if (j >= RLXG_STATE_SIZE)
        {
            break;
        }
        else if (lc == 1)
        {
            continue;
        }
        else
        {
            for (i=0;i<nf;i++)
            {
                if (sscanf(field[i],"%d",state+j) > 0)
                {
                    j++;
                }
                else
                {
                    strbuf errmsg;
                    sprintf(errmsg,"error while parsing random generator state (%s:%d)",\
                    fname,lc);
                    LATAN_ERROR(errmsg,LATAN_ELATSYN);
                }
            }
            continue;
        }
    }
    END_FOR_LINE_TOK_F(field);
    
    return LATAN_SUCCESS;
}
/*                          resampled sample I/O                            */
/****************************************************************************/
latan_errno rs_sample_save_ascii(const strbuf fname, const char mode,\
                                 const rs_sample *s)
{
    int thread;
    size_t i;
    size_t nsample;
    strbuf smode,name;

#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif
    smode[0] = mode;
    smode[1] = '\0';
    nsample  = rs_sample_get_nsample(s);

    if (mode == 'w')
    {
        ascii_new_file_buf(fname);
    }
    else
    {
        LATAN_ERROR("only 'w' file mode is authorized for saving sample in ASCII format",\
                    LATAN_EINVAL);
    }
    rs_sample_get_name(name,s);
    fprintf(env.f_buf[thread],"# latan_resampled_sample %s\n",name);
    fprintf(env.f_buf[thread],"%lu\n",(long unsigned int)nsample);
    mat_dump(env.f_buf[thread],rs_sample_pt_cent_val(s),"%.15e");
    for (i=0;i<nsample;i++)
    {
        mat_dump(env.f_buf[thread],rs_sample_pt_sample(s,i),"%.15e");
    }
    fprintf(env.f_buf[thread],"\n");
    
    return LATAN_SUCCESS;
}

latan_errno rs_sample_load_nrow_ascii(size_t *nr, const strbuf fname,\
                                      const strbuf name)
{
    int thread,nf,lc,sec,nsample;
    int i;
    strbuf *field;
    double buf;
    bool got_nsample;

#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif
    field       = NULL;
    name        = NULL;
    sec         = 0;
    i           = 0;
    got_nsample = false;
    nsample     = 0;

    ascii_open_file_buf(fname,'r');
    BEGIN_FOR_LINE_TOK_F(field,env.f_buf[thread]," ",nf,lc)
    {
        if (lc == 1)
        {
            continue;
        }
        else if (!got_nsample)
        {
            if (sscanf(field[0],"%d",&nsample) > 0)
            {
                got_nsample = true;
            }
            else
            {
                strbuf errmsg;
                sprintf(errmsg,"error while reading number of samples (%s:%d)",\
                        fname,lc);
                LATAN_ERROR(errmsg,LATAN_ELATSYN);
            }
            continue;
        }
        else if (got_nsample)
        {
            for (i=0;i<nf;i++)
            {
                if (sscanf(field[i],"%lf",&buf) > 0)
                {
                    sec++;
                }
                else
                {
                    strbuf errmsg;
                    sprintf(errmsg,"error while counting sample elements (%s:%d)",\
                            fname,lc);
                    LATAN_ERROR(errmsg,LATAN_ELATSYN);
                }
            }
            continue;
        }
    }
    END_FOR_LINE_TOK_F(field);
    if ( (sec%(nsample+1)) != 0)
    {
        LATAN_ERROR("error while getting number of sample rows",LATAN_ELATSYN);
    }
    else
    {
        *nr = (size_t)(sec/(nsample+1));
    }
    
    return LATAN_SUCCESS;
}

latan_errno rs_sample_load_nsample_ascii(size_t *nsample, const strbuf fname,\
                                         const strbuf name)
{
    int thread,nf,lc,buf;
    strbuf *field;

#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif
    field = NULL;
    name  = NULL;

    ascii_open_file_buf(fname,'r');
    BEGIN_FOR_LINE_TOK_F(field,env.f_buf[thread]," ",nf,lc)
    {
        if (lc == 1)
        {
            continue;
        }
        else if ((nf >= 1)&&(sscanf(field[0],"%d",&buf) > 0))
        {
            *nsample = (size_t)buf;
            break;
        }
        else
        {
            strbuf errmsg;
            sprintf(errmsg,"error while reading number of samples (%s:%d)",\
                    fname,lc);
            LATAN_ERROR(errmsg,LATAN_ELATSYN);
        }
    }
    END_FOR_LINE_TOK_F(field);

    return LATAN_SUCCESS;
}

latan_errno rs_sample_load_ascii(rs_sample *s, const strbuf fname,\
                                 const strbuf name)
{
    int thread,nf,lc,j,jmod,nsample,nr;
    int i;
    strbuf *field;
    double buf;
    bool got_nsample;

#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif
    field       = NULL;
    name        = NULL;
    i           = 0;
    j           = 0;
    jmod        = 0;
    got_nsample = false;
    nsample     = 0;
    nr          = (int)rs_sample_get_nrow(s);

    ascii_open_file_buf(fname,'r');
    BEGIN_FOR_LINE_TOK_F(field,env.f_buf[thread]," ",nf,lc)
    {
        if (lc == 1)
        {
            continue;
        }
        else if (!got_nsample)
        {
            if (sscanf(field[0],"%d",&nsample) > 0)
            {
                got_nsample = true;
            }
            else
            {
                strbuf errmsg;
                sprintf(errmsg,"error while reading number of samples (%s:%d)",\
                        fname,lc);
                LATAN_ERROR(errmsg,LATAN_ELATSYN);
            }
            continue;
        }
        else if (got_nsample)
        {
            for (i=0;i<nf;i++)
            {
                if (sscanf(field[i],"%lf",&buf) > 0)
                {
                    jmod = j%nr;
                    if (jmod == j)
                    {
                        mat_set(rs_sample_pt_cent_val(s),(size_t)jmod,0,buf);
                    }
                    else
                    {
                        mat_set(rs_sample_pt_sample(s,(size_t)(j/nr-1)),\
                                (size_t)jmod,0,buf);
                    }
                    j++;
                }
                else
                {
                    strbuf errmsg;
                    sprintf(errmsg,"error while parsing sample (%s:%d)",\
                            fname,lc);
                    LATAN_ERROR(errmsg,LATAN_ELATSYN);
                }
            }
            continue;
        }
    }
    END_FOR_LINE_TOK_F(field);

    return LATAN_SUCCESS;
}