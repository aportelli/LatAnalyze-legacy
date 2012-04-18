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
#ifndef LATAN_BEGIN
#define LATAN_BEGIN "latan_begin"
#endif
#ifndef LATAN_END
#define LATAN_END "latan_end"
#endif
#ifndef LATAN_MAT
#define LATAN_MAT "mat"
#endif
#ifndef LATAN_RG_STATE
#define LATAN_RG_STATE "rg_state"
#endif
#ifndef LATAN_RS_SAMPLE
#define LATAN_RS_SAMPLE "rs_sample"
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
        else if (mode == 'r')
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

/*                   parsing kernel declarations                            */
/****************************************************************************/
typedef struct mat_ker_state_s
{
    int j,nr,nc;
    bool got_ncol;
} mat_ker_state;

typedef struct rs_sample_ker_state_s
{
    int ns,i;
    bool got_ns;
    strbuf read_name,sname;
    mat *pt;
    mat_ker_state sampks;
} rs_sample_ker_state;

static latan_errno mat_load_ascii_ker(mat *m, size_t *dim, const strbuf fname,\
                                      const strbuf name, strbuf *field,       \
                                      const size_t nf, const int lc,          \
                                      bool *is_inmat, bool *is_end,           \
                                      mat_ker_state *ks);
static latan_errno randgen_load_state_ascii_ker(rg_state state,                \
                                                const strbuf fname,            \
                                                const strbuf name,             \
                                                strbuf *field, const size_t nf,\
                                                const int lc, bool *is_inrgs,  \
                                                bool *is_end,                  \
                                                int *j);
static latan_errno rs_sample_load_ascii_ker(rs_sample *s, size_t *nsample,    \
                                            size_t *dim, const strbuf fname,  \
                                            const strbuf name, strbuf *field, \
                                            const size_t nf, const int lc,    \
                                            bool *is_inrss, bool *is_end,     \
                                            bool *is_insamp, bool *is_sampend,\
                                            rs_sample_ker_state *ks);

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
    fprintf(FILE_BUF(thread),"%s %s %s %s\n",LATAN_COMMENT,LATAN_BEGIN,\
            LATAN_MAT,name);
    fprintf(FILE_BUF(thread),"%lu\n",(long unsigned int)ncol(m));
    mat_dump(FILE_BUF(thread),m,"% .15e");
    fprintf(FILE_BUF(thread),"%s %s %s\n",LATAN_COMMENT,LATAN_END,LATAN_MAT);
    
    return LATAN_SUCCESS;
}

latan_errno mat_load_ascii(mat *m, size_t *dim, const strbuf fname,\
                           const strbuf name)
{
    latan_errno status;
    int thread,nf,lc;
    strbuf *field;
    bool is_inmat,is_end;
    mat_ker_state ks;
    
#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif
    status   = LATAN_SUCCESS;
    field    = NULL;
    is_inmat = false;
    is_end   = false;
    
    ascii_open_file_buf(fname,'r');
    BEGIN_FOR_LINE_TOK_F(field,FILE_BUF(thread)," ",nf,lc)
    {
        USTAT(mat_load_ascii_ker(m,dim,fname,name,field,nf,lc,&is_inmat,\
                                 &is_end,&ks));
        if (is_end)
        {
            break;
        }
    }
    END_FOR_LINE_TOK_F(field);
    if (!is_end)
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
        if (is_inmat)
        {
            sprintf(errmsg,"unexpected EOF in %s while reading matrix (name= %s)",\
                    fname,buf);
            LATAN_ERROR(errmsg,LATAN_ELATSYN);
        }
        else
        {
            sprintf(errmsg,"matrix (name= %s) not found in file %s",buf,fname);
            LATAN_ERROR(errmsg,LATAN_EINVAL);
        }
    }
    
    return status;
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
    fprintf(FILE_BUF(thread),"%s %s %s %s\n",LATAN_COMMENT,LATAN_BEGIN,\
            LATAN_RG_STATE,name);
    for (i=0;i<RLXG_STATE_SIZE;i++)
    {
        fprintf(FILE_BUF(thread),"%d\n",state[i]);
    }
    fprintf(FILE_BUF(thread),"%s %s %s\n",LATAN_COMMENT,LATAN_END,\
           LATAN_RG_STATE);
    
    return LATAN_SUCCESS;
}

latan_errno randgen_load_state_ascii(rg_state state, const strbuf fname,\
                                     const strbuf name)
{
    latan_errno status;
    int thread,nf,lc;
    strbuf *field;
    bool is_inrgs,is_end;
    int j;

#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif
    status   = LATAN_SUCCESS;
    field    = NULL;
    is_inrgs = false;
    is_end   = false;
    j        = 0;

    ascii_open_file_buf(fname,'r');
    BEGIN_FOR_LINE_TOK_F(field,FILE_BUF(thread)," ",nf,lc)
    {
        USTAT(randgen_load_state_ascii_ker(state,fname,name,field,nf,lc,\
                                           &is_inrgs,&is_end,&j));
        if (is_end)
        {
            break;
        }
    }
    END_FOR_LINE_TOK_F(field);
    if (!is_end)
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
        if (is_inrgs)
        {
            sprintf(errmsg,"unexpected EOF in %s while reading random generator state (name= %s)",\
                    fname,buf);
            LATAN_ERROR(errmsg,LATAN_ELATSYN);
        }
        else
        {
            sprintf(errmsg,"random generator state (name= %s) not found in file %s",\
                    buf,fname);
            LATAN_ERROR(errmsg,LATAN_EINVAL);
        }
    }
    
    return status;
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
    fprintf(FILE_BUF(thread),"%s %s %s %s\n",LATAN_COMMENT,LATAN_BEGIN,\
            LATAN_RS_SAMPLE,name);
    fprintf(FILE_BUF(thread),"%lu\n",(long unsigned int)nsample);
    sprintf(sname,"%s_C",name);
    USTAT(mat_save_ascii(fname,'a',rs_sample_pt_cent_val(s),sname));
    for (i=0;i<nsample;i++)
    {
        sprintf(sname,"%s_S_%lu",name,(long unsigned int)(i));
        USTAT(mat_save_ascii(fname,'a',rs_sample_pt_sample(s,i),sname));
    }
    fprintf(FILE_BUF(thread),"%s %s %s\n",LATAN_COMMENT,LATAN_END,\
            LATAN_RS_SAMPLE);
    
    return status;
}

latan_errno rs_sample_load_ascii(rs_sample *s, size_t *nsample, size_t *dim,\
                                 const strbuf fname, const strbuf name)
{
    latan_errno status;
    int thread,nf,lc;
    strbuf *field;
    bool is_inrss,is_end,is_insamp,is_sampend;
    rs_sample_ker_state ks;

#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif
    status     = LATAN_SUCCESS;
    field      = NULL;
    is_inrss   = false;
    is_end     = false;
    is_insamp  = false;
    is_sampend = false;

    ascii_open_file_buf(fname,'r');
    BEGIN_FOR_LINE_TOK_F(field,FILE_BUF(thread)," ",nf,lc)
    {
        USTAT(rs_sample_load_ascii_ker(s,nsample,dim,fname,name,field,nf,lc,\
                                       &is_inrss,&is_end,&is_insamp,        \
                                       &is_sampend,&ks));
        if (is_end)
        {
            break;
        }
    }
    END_FOR_LINE_TOK_F(field);
    if (!is_sampend)
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
        if (is_insamp)
        {
            sprintf(errmsg,"unexpected EOF in %s while reading sample %d (name= %s)",\
                    fname,ks.i,buf);
            LATAN_ERROR(errmsg,LATAN_ELATSYN);
        }
        else
        {
            sprintf(errmsg,"sample %d (name= %s) not found in file %s",ks.i,\
                    buf,fname);
            LATAN_ERROR(errmsg,LATAN_EINVAL);
        }
    }
    if (!is_end)
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
        if (is_inrss)
        {
            sprintf(errmsg,"unexpected EOF in %s while reading samples (name= %s)",\
                    fname,buf);
            LATAN_ERROR(errmsg,LATAN_ELATSYN);
        }
        else
        {
            sprintf(errmsg,"samples (name= %s) not found in file %s",buf,fname);
            LATAN_ERROR(errmsg,LATAN_EINVAL);
        }
    }
    
    return status;
}

/*                        parsing kernels (internal)                        */
/****************************************************************************/
static latan_errno mat_load_ascii_ker(mat *m, size_t *dim, const strbuf fname,\
                                      const strbuf name, strbuf *field,       \
                                      const size_t nf, const int lc,          \
                                      bool *is_inmat, bool *is_end,           \
                                      mat_ker_state *ks)
{
    size_t i;
    double dbuf;
    
    if ((nf >= 4)&&!*is_inmat)
    {
        *is_inmat = true;
        *is_inmat = *is_inmat && (strbufcmp(field[0],LATAN_COMMENT) == 0);
        *is_inmat = *is_inmat && (strbufcmp(field[1],LATAN_BEGIN) == 0);
        *is_inmat = *is_inmat && (strbufcmp(field[2],LATAN_MAT) == 0);
        if (strlen(name) != 0)
        {
            *is_inmat  = *is_inmat && (strbufcmp(field[3],name) == 0);
        }
        if (*is_inmat)
        {
            *is_end      = false;
            ks->got_ncol = false;
            ks->nc       = 0;
            ks->j        = 0;
        }
    }
    else if (*is_inmat)
    {
        if ((nf >= 3)&&(strbufcmp(field[0],LATAN_COMMENT) == 0))
        {
            *is_end = true;
            *is_end = *is_end && (strbufcmp(field[1],LATAN_END) == 0);
            *is_end = *is_end && (strbufcmp(field[2],LATAN_MAT) == 0);
            if (*is_end)
            {
                *is_inmat = false;
                if ((ks->j)%(ks->nc) != 0)
                {
                    strbuf errmsg;
                    sprintf(errmsg,"matrix parsing unexpected end (%s:%d)",\
                            fname,lc);
                    LATAN_ERROR(errmsg,LATAN_ELATSYN);
                }
                else
                {
                    ks->nr = (ks->nc != 0) ? (ks->j)/(ks->nc) : 0;
                    if (m)
                    {
                        if (nrow(m) != (size_t)(ks->nr))
                        {
                            strbuf errmsg;
                            sprintf(errmsg,"row number mismatch (%s:%d)",\
                                    fname,lc);
                            LATAN_ERROR(errmsg,LATAN_EBADLEN);
                        }
                    }
                    if (dim)
                    {
                        dim[0]   = (size_t)(ks->nr);
                        dim[1]   = (size_t)(ks->nc);
                    }
                }
            }
            else
            {
                strbuf errmsg;
                sprintf(errmsg,"matrix parsing unexpected end (%s:%d)",\
                        fname,lc);
                LATAN_ERROR(errmsg,LATAN_ELATSYN);
            }
        }
        else if (!(ks->got_ncol))
        {
            if (sscanf(field[0],"%d",&(ks->nc)) > 0)
            {
                if (m)
                {
                    if (ncol(m) != (size_t)(ks->nc))
                    {
                        strbuf errmsg;
                        sprintf(errmsg,"column number mismatch (%s:%d)",\
                                fname,lc);
                        LATAN_ERROR(errmsg,LATAN_EBADLEN);
                    }
                }
                ks->got_ncol = true;
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
                if (field[i][0] == '#')
                {
                    break;
                }
                else if (sscanf(field[i],"%lf",&dbuf) > 0)
                {
                    if (m)
                    {
                        if (ks->j >= (int)nel(m))
                        {
                            strbuf errmsg;
                            sprintf(errmsg,"row number mismatch (%s:%d)",\
                                    fname,lc);
                            LATAN_ERROR(errmsg,LATAN_EBADLEN);
                        }
                        mat_set(m,(size_t)((ks->j)/(ks->nc)), \
                                (size_t)((ks->j)%(ks->nc)),dbuf);
                    }
                    (ks->j)++;
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
    
    return LATAN_SUCCESS;
}

static latan_errno randgen_load_state_ascii_ker(rg_state state,                \
                                                const strbuf fname,            \
                                                const strbuf name,             \
                                                strbuf *field, const size_t nf,\
                                                const int lc, bool *is_inrgs,  \
                                                bool *is_end,                  \
                                                int *j)
{
    size_t i;
    
    if ((nf >= 4)&&!*is_inrgs)
    {
        *is_inrgs = true;
        *is_inrgs = *is_inrgs&&(strbufcmp(field[0],LATAN_COMMENT) == 0);
        *is_inrgs = *is_inrgs&&(strbufcmp(field[1],LATAN_BEGIN) == 0);
        *is_inrgs = *is_inrgs&&(strbufcmp(field[2],LATAN_RG_STATE) == 0);
        if (strlen(name) != 0)
        {
            *is_inrgs = *is_inrgs&&(strbufcmp(field[3],name) == 0);
        }
        if (*is_inrgs)
        {
            *is_end  = false;
            *j       = 0;
        }
    }
    else if (*is_inrgs)
    {
        if ((nf >= 3)&&(strbufcmp(field[0],LATAN_COMMENT) == 0))
        {
            *is_end = true;
            *is_end = *is_end && (strbufcmp(field[1],LATAN_END) == 0);
            *is_end = *is_end && (strbufcmp(field[2],LATAN_RG_STATE) == 0);
            if (*is_end)
            {
                *is_inrgs = false;
                if (*j != RLXG_STATE_SIZE)
                {
                    strbuf errmsg;
                    sprintf(errmsg,"random generator state parsing unexpected end (%s:%d)",\
                            fname,lc);
                    LATAN_ERROR(errmsg,LATAN_ELATSYN);
                }
            }
        }
        else
        {
            for (i=0;i<nf;i++)
            {
                if (sscanf(field[i],"%d",state+(*j)) > 0)
                {
                    (*j)++;
                    if (*j == RLXG_STATE_SIZE)
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
    
    return LATAN_SUCCESS;
}

static latan_errno rs_sample_load_ascii_ker(rs_sample *s, size_t *nsample,    \
                                            size_t *dim, const strbuf fname,  \
                                            const strbuf name, strbuf *field, \
                                            const size_t nf, const int lc,    \
                                            bool *is_inrss, bool *is_end,     \
                                            bool *is_insamp, bool *is_sampend,\
                                            rs_sample_ker_state *ks)
{
    latan_errno status;
    bool bbuf;
    
    status = LATAN_SUCCESS;
    
    if ((nf >= 4)&&!*is_inrss)
    {
        *is_inrss = true;
        *is_inrss = *is_inrss&&(strbufcmp(field[0],LATAN_COMMENT) == 0);
        *is_inrss = *is_inrss&&(strbufcmp(field[1],LATAN_BEGIN) == 0);
        *is_inrss = *is_inrss&&(strbufcmp(field[2],LATAN_RS_SAMPLE) == 0);
        if (strlen(name))
        {
            *is_inrss = *is_inrss&&(strbufcmp(field[3],name) == 0);
            strbufcpy(ks->read_name,name);
        }
        else
        {
            strbufcpy(ks->read_name,field[3]);
        }
        if (*is_inrss)
        {
            ks->ns      = 0;
            ks->i       = -1;
            *is_end     = false;
            *is_insamp  = false;
            *is_sampend = false;
            ks->got_ns  = false;
            bbuf        = false;
            sprintf(ks->sname,"%s_C",ks->read_name);
            if (s)
            {
                ks->pt = rs_sample_pt_cent_val(s);
            }
            else
            {
                ks->pt = NULL;
            }
        }
    }
    else if (*is_inrss)
    {
        if ((nf >= 3)&&(strbufcmp(field[0],LATAN_COMMENT) == 0))
        {
            *is_end = true;
            *is_end = *is_end && (strbufcmp(field[1],LATAN_END) == 0);
            *is_end = *is_end && (strbufcmp(field[2],LATAN_RS_SAMPLE) == 0);
            if (*is_end)
            {
                *is_inrss = false;
                if (s)
                {
                    if (ks->i != ks->ns)
                    {
                        LATAN_ERROR("sample number mismatch",LATAN_EBADLEN);
                    }
                }
            }
        }
        if (!(*is_end))
        {
            if (!(ks->got_ns))
            {
                if (sscanf(field[0],"%d",&(ks->ns)) > 0)
                {
                    if (nsample)
                    {
                        *nsample = (size_t)(ks->ns);
                    }
                    ks->got_ns = true;
                    if (s)
                    {
                        if (rs_sample_get_nsample(s) != (size_t)(ks->ns))
                        {
                            LATAN_ERROR("sample number mismatch",LATAN_EBADLEN);
                        }
                    }
                }
                else
                {
                    strbuf errmsg;
                    sprintf(errmsg,"error while reading number of samples (%s:%d)",\
                            fname,lc);
                    LATAN_ERROR(errmsg,LATAN_ELATSYN);
                }
            }
            else
            {
                USTAT(mat_load_ascii_ker(ks->pt,dim,fname,ks->sname,field,nf,\
                                         lc,&bbuf,is_sampend,&(ks->sampks)));
                if ((*is_sampend)&&(*is_insamp))
                {
                    if (s)
                    {
                        (ks->i)++;
                        ks->pt = rs_sample_pt_sample(s,(size_t)(ks->i));
                        sprintf(ks->sname,"%s_S_%lu",ks->read_name,\
                                (long unsigned int)(ks->i));
                    }
                    else
                    {
                        *is_end = true;
                    }
                }
                *is_insamp = bbuf;
            }
        }
    }
    
    return status;
}
