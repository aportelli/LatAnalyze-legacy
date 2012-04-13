/* latan_io_ascii.c, part of LatAnalyze library
 *
 * Copyright (C) 2010, 2011, 2012 Antonin Portelli
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

#define _POSIX_C_SOURCE 199506L /* strtok_r is used here */

#ifndef LATAN_COMMENT
#define LATAN_COMMENT "#L"
#endif
#ifndef LATAN_MAT
#define LATAN_MAT "latan_mat"
#endif
#ifndef LATAN_RG_STATE
#define LATAN_RG_STATE "latan_rg_state"
#endif
#ifndef LATAN_RS_SAMPLE
#define LATAN_RS_SAMPLE "latan_rs_sample"
#endif

#include <latan/latan_io_ascii.h>
#include <latan/latan_includes.h>
#include <latan/latan_io.h>

static latan_errno ascii_open_file_buf(const strbuf fname, const char mode);

/*                       file buffer management (internal)                  */
/****************************************************************************/
typedef struct
{
    FILE *f;
    strbuf fname;
    char mode;
} ascii_file;

typedef struct
{
    ascii_file *ascii_buf;
    bool *file_is_loaded;
    int nfile;
} io_ascii_env;

static io_ascii_env env =
{
    NULL,\
    NULL,\
    0    \
};

#define FILE_BUF(thread)  env.ascii_buf[thread].f
#define FILE_NAME(thread) env.ascii_buf[thread].fname
#define FILE_MODE(thread) env.ascii_buf[thread].mode

static latan_errno ascii_open_file_buf(const strbuf fname, const char mode)
{
    latan_errno status;
    strbuf smode,errmsg;
    int nthread,thread,i;

#ifdef _OPENMP
    nthread = omp_get_num_threads();
    thread  = omp_get_thread_num();
#else
    nthread = 1;
    thread  = 0;
#endif
    status   = LATAN_SUCCESS;
    switch (mode)
    {
        case 'r' :
            strbufcpy(smode,"r");
            break;
        case 'w' :
            strbufcpy(smode,"w+");
            break;
        case 'a' :
            strbufcpy(smode,"a+");
            break;
        default:
            sprintf(errmsg,"ASCII file mode %c unknown",mode);
            LATAN_ERROR(errmsg,LATAN_EINVAL);
            break;
    }

#ifdef _OPENMP
    #pragma omp critical
#endif
    {
        if (nthread > env.nfile)
        {
            REALLOC_NOERRET(env.ascii_buf,env.ascii_buf,ascii_file *,nthread);
            REALLOC_NOERRET(env.file_is_loaded,env.file_is_loaded,bool *,\
                            nthread);
            for (i=env.nfile;i<nthread;i++)
            {
                env.file_is_loaded[i] = false;
                strbufcpy(FILE_NAME(i),"");
                FILE_BUF(i)           = NULL;
                FILE_MODE(i)          = '\0';
            }
            env.nfile = nthread;
        }
    }
    if (env.file_is_loaded[thread])
    {
        if ((strbufcmp(FILE_NAME(thread),fname) != 0)||(mode == 'w')||\
            (mode != FILE_MODE(thread)))
        {
            fclose(FILE_BUF(thread));
            FOPEN(FILE_BUF(thread),fname,smode);
            strbufcpy(FILE_NAME(thread),fname);
            FILE_MODE(thread) = mode;
        }
        else
        {
            rewind(FILE_BUF(thread));
        }
    }
    else
    {
        FOPEN(FILE_BUF(thread),fname,smode);
        strbufcpy(FILE_NAME(thread),fname);
        FILE_MODE(thread) = mode;
        env.file_is_loaded[thread] = true;
    }

    return status;
}

/*                              I/O init/finish                             */
/****************************************************************************/
void io_init_ascii(void)
{
#ifdef _OPENMP
    if(omp_in_parallel())
    {
        LATAN_WARNING("I/O initialization called from a parallel region",\
                      LATAN_FAILURE);
    }
#endif
}

void io_finish_ascii(void)
{
    int i;

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
            fclose(FILE_BUF(i));
            strbufcpy(FILE_NAME(i),"");
            env.file_is_loaded[i] = false;
        }
    }
    FREE(env.ascii_buf);
    FREE(env.file_is_loaded);
}

/*                             mat I/O                                      */
/****************************************************************************/
latan_errno mat_save_ascii(const strbuf fname, const char mode, const mat *m,\
                           const strbuf name)
{
    int thread;
    
#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif
    
    if ((mode == 'w')||(mode == 'a'))
    {
        ascii_open_file_buf(fname,mode);
    }
    else
    {
        LATAN_ERROR("unknown or read-only file mode",LATAN_EINVAL);
    }
    fprintf(FILE_BUF(thread),"%s %s %s\n",LATAN_COMMENT,LATAN_MAT,name);
    fprintf(FILE_BUF(thread),"%lu\n",(long unsigned int)ncol(m));
    mat_dump(FILE_BUF(thread),m,"% .15e");
    fprintf(FILE_BUF(thread),"\n");
    
    return LATAN_SUCCESS;
}

latan_errno mat_load_ascii(mat *m, size_t *dim, const strbuf fname,\
                           const strbuf name)
{
    int thread,nf,lc,nc,nr;
    int i,j;
    strbuf *field;
    double dbuf;
    bool got_ncol,is_inmat;
    
#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif
    field    = NULL;
    nc       = 0;
    j        = 0;
    got_ncol = false;
    is_inmat = false;
    
    ascii_open_file_buf(fname,'r');
    BEGIN_FOR_LINE_TOK_F(field,FILE_BUF(thread)," ",nf,lc)
    {
        if ((nf >= 3)&&!is_inmat)
        {
            is_inmat  = true;
            is_inmat  = is_inmat && (strbufcmp(field[0],LATAN_COMMENT) == 0);
            is_inmat  = is_inmat && (strbufcmp(field[1],LATAN_MAT) == 0);
            if (strlen(name) != 0)
            {
                is_inmat  = is_inmat && (strbufcmp(field[2],name) == 0);
            }
        }
        else if (nf > 0)
        {
            if (strbufcmp(field[0],LATAN_COMMENT) == 0)
            {
                break;
            }
            else if (field[0][0] == '#')
            {
                continue;
            }
            else if (!got_ncol)
            {
                if (sscanf(field[0],"%d",&nc) > 0)
                {
                    if (m)
                    {
                        if (ncol(m) != (size_t)(nc))
                        {
                            LATAN_ERROR("column number mismatch",LATAN_EBADLEN);
                        }
                    }
                    if (dim)
                    {
                        dim[1]   = (size_t)(nc);
                    }
                    got_ncol = true;
                }
                else
                {
                    strbuf errmsg;
                    sprintf(errmsg,"error while reading column number (%s:%d)",\
                            fname,lc);
                    LATAN_ERROR(errmsg,LATAN_ELATSYN);
                }
            }
            else
            {
                for (i=0;i<nf;i++)
                {
                    if (sscanf(field[i],"%lf",&dbuf) > 0)
                    {
                        if (m)
                        {
                            if (j >= (int)nel(m))
                            {
                                LATAN_ERROR("row number mismatch",\
                                            LATAN_EBADLEN);
                            }
                            mat_set(m,(size_t)(j/nc),(size_t)(j%nc),dbuf);
                        }
                        j++;
                    }
                    else
                    {
                        strbuf errmsg;
                        sprintf(errmsg,"impossible to read matrix element (%s:%d)",\
                                fname,lc);
                        LATAN_ERROR(errmsg,LATAN_ELATSYN);
                    }
                }
            }
        }
    }
    END_FOR_LINE_TOK_F(field);
    if (!is_inmat)
    {
        strbuf errmsg,buf;
        
        if (strlen(name) == 0)
        {
            strcpy(buf,"<no_name>");
        }
        else
        {
            sprintf(buf,"\"%s\"",name);
        }
        sprintf(errmsg,"matrix (name= %s) not found in file %s",buf,fname);
        LATAN_ERROR(errmsg,LATAN_EINVAL);
    }
    if (j%nc != 0)
    {
        strbuf errmsg;
        sprintf(errmsg,"matrix parsing unexpected end (%s:%d)",fname,lc);
        LATAN_ERROR(errmsg,LATAN_ELATSYN);
    }
    else
    {
        nr = j/nc;
        if (m)
        {
            if (nrow(m) != (size_t)(nr))
            {
                LATAN_ERROR("row number mismatch",LATAN_EBADLEN);
            }
        }
        if (dim)
        {
            dim[0]   = (size_t)(nr);
        }
    }
    
    return LATAN_SUCCESS;
}

/*                      random generator state I/O                          */
/****************************************************************************/
latan_errno randgen_save_state_ascii(const strbuf fname, const char mode,   \
                                     const rg_state state, const strbuf name)
{
    int thread;
    int i;
    
#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif

    if ((mode == 'w')||(mode == 'a'))
    {
        ascii_open_file_buf(fname,mode);
    }
    else
    {
        LATAN_ERROR("unknown or read-only file mode",LATAN_EINVAL);
    }
    fprintf(FILE_BUF(thread),"%s %s %s\n",LATAN_COMMENT,LATAN_RG_STATE,name);
    for (i=0;i<RLXG_STATE_SIZE;i++)
    {
        fprintf(FILE_BUF(thread),"%d\n",state[i]);
    }
    fprintf(FILE_BUF(thread),"\n");
    
    return LATAN_SUCCESS;
}

latan_errno randgen_load_state_ascii(rg_state state, const strbuf fname,\
                                     const strbuf name)
{
    int thread,nf,lc;
    int i,j;
    strbuf *field;
    bool is_inrgs;

#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif
    field    = NULL;
    j        = 0;
    is_inrgs = false;

    ascii_open_file_buf(fname,'r');
    BEGIN_FOR_LINE_TOK_F(field,FILE_BUF(thread)," ",nf,lc)
    {
        if ((nf >= 3)&&!is_inrgs)
        {
            is_inrgs = true;
            is_inrgs = is_inrgs&&(strbufcmp(field[0],LATAN_COMMENT) == 0);
            is_inrgs = is_inrgs&&(strbufcmp(field[1],LATAN_RG_STATE) == 0);
            if (strlen(name) != 0)
            {
                is_inrgs = is_inrgs&&(strbufcmp(field[2],name) == 0);
            }
        }
        else if (nf > 0)
        {
            if (strbufcmp(field[0],LATAN_COMMENT) == 0)
            {
                break;
            }
            else if (field[0][0] == '#')
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
                        if (j == RLXG_STATE_SIZE)
                        {
                            break;
                        }
                    }
                    else
                    {
                        strbuf errmsg;
                        sprintf(errmsg,"impossible to read random generator state element (%s:%d)",\
                                fname,lc);
                        LATAN_ERROR(errmsg,LATAN_ELATSYN);
                    }
                }
            }
        }
    }
    END_FOR_LINE_TOK_F(field);
    if (!is_inrgs)
    {
        strbuf errmsg,buf;
        
        if (strlen(name) == 0)
        {
            strcpy(buf,"<no_name>");
        }
        else
        {
            sprintf(buf,"\"%s\"",name);
        }
        sprintf(errmsg,"random generator state (name= %s) not found in file %s",\
                buf,fname);
        LATAN_ERROR(errmsg,LATAN_EINVAL);
    }
    if (j != RLXG_STATE_SIZE)
    {
        strbuf errmsg;
        sprintf(errmsg,"random generator state parsing unexpected end (%s:%d)",\
                fname,lc);
        LATAN_ERROR(errmsg,LATAN_ELATSYN);
    }
    
    return LATAN_SUCCESS;
}
/*                          resampled sample I/O                            */
/****************************************************************************/
latan_errno rs_sample_save_ascii(const strbuf fname, const char mode,\
                                 const rs_sample *s, const strbuf name)
{
    latan_errno status;
    int thread;
    size_t nsample;
    size_t i;
    strbuf sname;
    
#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif
    status   = LATAN_SUCCESS;
    nsample  = rs_sample_get_nsample(s);

    if ((mode == 'w')||(mode == 'a'))
    {
        ascii_open_file_buf(fname,mode);
    }
    else
    {
        LATAN_ERROR("unknown or read-only file mode",LATAN_EINVAL);
    }
    fprintf(FILE_BUF(thread),"%s %s %s\n",LATAN_COMMENT,LATAN_RS_SAMPLE,name);
    fprintf(FILE_BUF(thread),"%lu\n",(long unsigned int)nsample);
    sprintf(sname,"%s_C",name);
    USTAT(mat_save_ascii(fname,'a',rs_sample_pt_cent_val(s),sname));
    for (i=0;i<nsample;i++)
    {
        sprintf(sname,"%s_S_%lu",name,(long unsigned int)(i));
        USTAT(mat_save_ascii(fname,'a',rs_sample_pt_sample(s,i),sname));
    }
    
    return status;
}

latan_errno rs_sample_load_ascii(rs_sample *s, size_t *nsample, size_t *dim,\
                                 const strbuf fname, const strbuf name)
{
    latan_errno status;
    int thread,nf,ns,lc;
    size_t cvdim[2];
    size_t i;
    strbuf *field,read_name,sname;
    bool is_inrss;

#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif
    status   = LATAN_SUCCESS;
    field    = NULL;
    is_inrss = false;

    ascii_open_file_buf(fname,'r');
    BEGIN_FOR_LINE_TOK_F(field,FILE_BUF(thread)," ",nf,lc)
    {
        if ((nf >= 3)&&!is_inrss)
        {
            is_inrss = true;
            is_inrss = is_inrss&&(strbufcmp(field[0],LATAN_COMMENT) == 0);
            is_inrss = is_inrss&&(strbufcmp(field[1],LATAN_RS_SAMPLE) == 0);
            if (strlen(name))
            {
                is_inrss = is_inrss&&(strbufcmp(field[2],name) == 0);
                strbufcpy(read_name,name);
            }
            else
            {
                strbufcpy(read_name,field[2]);
            }
        }
        else if (nf > 0)
        {
            if (strbufcmp(field[0],LATAN_COMMENT) == 0)
            {
                strbuf errmsg;
                sprintf(errmsg,"error while reading number of samples (%s:%d)",\
                        fname,lc);
                LATAN_ERROR(errmsg,LATAN_ELATSYN);
            }
            else if (field[0][0] == '#')
            {
                continue;
            }
            else if (sscanf(field[0],"%d",&ns) > 0)
            {
                if (nsample)
                {
                    *nsample = (size_t)(ns);
                }
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
    }
    END_FOR_LINE_TOK_F(field);
    if (!is_inrss)
    {
        strbuf errmsg,buf;
        
        if (strlen(name) == 0)
        {
            strcpy(buf,"<no_name>");
        }
        else
        {
            sprintf(buf,"\"%s\"",name);
        }
        sprintf(errmsg,"resampled samples (name= %s) not found in file %s",\
                buf,fname);
        LATAN_ERROR(errmsg,LATAN_EINVAL);
    }
    sprintf(sname,"%s_C",read_name);
    USTAT(mat_load_ascii(rs_sample_pt_cent_val(s),cvdim,fname,sname));
    if (dim)
    {
        dim[0] = cvdim[0];
        dim[1] = cvdim[1];
    }
    for (i=0;i<((size_t)ns);i++)
    {
        sprintf(sname,"%s_S_%lu",read_name,(long unsigned int)(i));
        USTAT(mat_load_ascii(rs_sample_pt_sample(s,i),NULL,fname,sname));
    }

    return status;
}
