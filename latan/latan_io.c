/* latan_io.c, part of LatAnalyze library
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

#include <latan/latan_io.h>
#include <latan/latan_includes.h>
#include <latan/latan_io_ascii.h>
#ifdef HAVE_LIBXML2
#include <latan/latan_io_xml.h>
#endif
#include <latan/latan_math.h>

/*                            default I/O functions                         */
/****************************************************************************/
#define IO_FUNC(func_name,suf) func_name##_##suf
#define DEF_IO_FMT                 IO_ASCII
#define DEF_IO_INIT                IO_FUNC(io_init,ascii)
#define DEF_IO_FINISH              IO_FUNC(io_finish,ascii)
#define DEF_MAT_SAVE               IO_FUNC(mat_save,ascii)
#define DEF_MAT_LOAD               IO_FUNC(mat_load,ascii)
#define DEF_RANDGEN_SAVE_STATE     IO_FUNC(randgen_save_state,ascii)
#define DEF_RANDGEN_LOAD_STATE     IO_FUNC(randgen_load_state,ascii)
#define DEF_RS_SAMPLE_SAVE         IO_FUNC(rs_sample_save,ascii)
#define DEF_RS_SAMPLE_LOAD         IO_FUNC(rs_sample_load,ascii)

/*                               environment                                */
/****************************************************************************/
typedef struct io_env_s
{
    io_fmt_no fmt;
    bool      is_init;
} io_env;

static io_env env =\
{                  \
    DEF_IO_FMT,    \
    false          \
};

/*                            function pointers                             */
/****************************************************************************/
static void (*io_init_pt)(void)   = &DEF_IO_INIT;
static void (*io_finish_pt)(void) = &DEF_IO_FINISH;
static latan_errno (*mat_save_pt)(const strbuf fname, const char mode,\
                                  const mat *m, const strbuf name)    \
    = &DEF_MAT_SAVE;
static latan_errno (*mat_load_pt)(mat *m, size_t *dim, const strbuf fname,\
                                  const strbuf name)                      \
    = &DEF_MAT_LOAD;
static latan_errno (*randgen_save_state_pt)(const strbuf f_name, \
                                            const char mode,     \
                                            const rg_state state,\
                                            const strbuf name)   \
    = &DEF_RANDGEN_SAVE_STATE;
static latan_errno (*randgen_load_state_pt)(rg_state state,     \
                                            const strbuf f_name,\
                                            const strbuf name)  \
    = &DEF_RANDGEN_LOAD_STATE;
static latan_errno (*rs_sample_save_pt)(const strbuf fname, const char mode,  \
                                        const rs_sample *s, const strbuf name)\
    = &DEF_RS_SAMPLE_SAVE;
static latan_errno (*rs_sample_load_pt)(rs_sample *s, size_t *nsample,  \
                                        size_t *dim, const strbuf fname,\
                                        const strbuf name)              \
    = &DEF_RS_SAMPLE_LOAD;

/*                              I/O format                                  */
/****************************************************************************/
#define SET_IO_FUNC(func_name,suf) func_name##_pt = &IO_FUNC(func_name,suf)
#define SET_IO_FUNCS(suf)\
SET_IO_FUNC(io_init,suf);\
SET_IO_FUNC(io_finish,suf);\
SET_IO_FUNC(mat_save,suf);\
SET_IO_FUNC(mat_load,suf);\
SET_IO_FUNC(randgen_save_state,suf);\
SET_IO_FUNC(randgen_load_state,suf);\
SET_IO_FUNC(rs_sample_save,suf);\
SET_IO_FUNC(rs_sample_load,suf);

latan_errno io_set_fmt(const io_fmt_no fmt)
{
    bool reinit;
    
    reinit = (io_get_fmt() != fmt);
    
    if (reinit)
    {
        io_finish();
    }
    env.fmt = fmt;
    switch (env.fmt)
    {
        case IO_XML:
#ifdef HAVE_LIBXML2
            SET_IO_FUNCS(xml);
#else
            LATAN_ERROR("XML support was not compiled",LATAN_EINVAL);
#endif
            break;
        case IO_ASCII:
            SET_IO_FUNCS(ascii);
            break;
        default:
            LATAN_ERROR("I/O format flag unknown",LATAN_EINVAL);
            break;
    }
    if (reinit)
    {
        io_init();
    }

    return LATAN_SUCCESS;
}
#undef SET_IO_FUNCS
#undef SET_IO_FUNC

io_fmt_no io_get_fmt(void)
{
    return env.fmt;
}

/*                              I/O init/finish                             */
/****************************************************************************/
void io_init(void)
{
    if (!env.is_init)
    {
        io_init_pt();
        env.is_init = true;
    }
}

void io_finish(void)
{
    if (env.is_init)
    {
        io_finish_pt();
        env.is_init = false;
    }
}

/*                                general I/O                               */
/****************************************************************************/
int get_nfile(const strbuf manifestfname)
{
    strbuf *field;
    int nfile,lc,nf;
    
    nfile = 0;
    field = NULL;
    
    BEGIN_FOR_LINE_TOK(field,manifestfname," \t",nf,lc)
    if ((nf>0)&&(field[0][0] != '#'))
    {
        nfile++;
    }
    END_FOR_LINE_TOK(field)
    
    return nfile;
}

void get_firstfname(strbuf fname, const strbuf manifestfname)
{
    strbuf *field;
    int lc,nf;
    
    field   = NULL;
    
    BEGIN_FOR_LINE_TOK(field,manifestfname," \t",nf,lc)
    if ((nf>0)&&(field[0][0] != '#'))
    {
        strbufcpy(fname,field[0]);
        break;
    }
    END_FOR_LINE_TOK(field)
}

void get_elname(strbuf fname, strbuf elname, const strbuf latan_path)
{
    char *pt;
    strbuf buf;
    
    strbufcpy(buf,latan_path);
    pt = strrchr(buf,LATAN_PATH_SEP);
    if (pt)
    {
        if (elname)
        {
            strbufcpy(elname,pt+1);
        }
        *pt = '\0';
        if (fname)
        {
            strbufcpy(fname,buf);
        }
    }
    else
    {
        if (fname)
        {
            strbufcpy(fname,latan_path);
        }
        if (elname)
        {
            strbufcpy(elname,"");
        }
    }
}

/*                              mat I/O                                     */
/****************************************************************************/
#define FUNC_INIT(fname,elname)\
{\
    get_elname(fname,elname,latan_path);\
    io_init();\
}

void mat_dump(FILE* stream, const mat *m, const strbuf fmt)
{
    size_t i,j;
    
    for (i=0;i<nrow(m);i++)
    {
        for (j=0;j<ncol(m)-1;j++)
        {
            fprintf(stream,fmt,mat_get(m,i,j));
            fprintf(stream," ");
        }
        fprintf(stream,fmt,mat_get(m,i,ncol(m)-1));
        fprintf(stream,"\n");
    }
}

latan_errno mat_save(const strbuf latan_path, const char mode, const mat *m)
{
    latan_errno status;
    strbuf fname,elname;
    
    FUNC_INIT(fname,elname);
    if (strlen(elname) == 0)
    {
        LATAN_ERROR("no name specified",LATAN_EINVAL);
    }
    status = mat_save_pt(fname,mode,m,elname);
    
    return status;
}

latan_errno mat_save_subm(const strbuf latan_path, const char mode,      \
                          const mat *m, const size_t k1, const size_t l1,\
                          const size_t k2, const size_t l2)
{
    latan_errno status;
    mat *little_m;
    
    status = LATAN_SUCCESS;
    
    little_m = mat_create(k2-k1+1,l2-l1+1);
    
    USTAT(mat_get_subm(little_m,m,k1,l1,k2,l2));
    USTAT(mat_save(latan_path,mode,little_m));
    
    mat_destroy(little_m);
    
    return status;
}

latan_errno mat_load(mat *m, size_t *dim, const strbuf latan_path)
{
    latan_errno status;
    strbuf fname,elname;
    
    FUNC_INIT(fname,elname);
    status = mat_load_pt(m,dim,fname,elname);
    
    return status;
}

latan_errno mat_load_subm(mat *m, const strbuf latan_path, const size_t k1, \
                          const size_t l1, const size_t k2, const size_t l2)
{
    mat *big_m;
    size_t dim[2];
    latan_errno status;
    
    status = LATAN_SUCCESS;
    
    USTAT(mat_load(NULL,dim,latan_path));
    big_m = mat_create(dim[0],dim[1]);
    USTAT(mat_load(big_m,NULL,latan_path));
    USTAT(mat_get_subm(m,big_m,k1,l1,k2,l2));
    
    mat_destroy(big_m);
    
    return status;
}

latan_errno mat_ar_loadbin(mat **m, size_t *dim, const strbuf man_fname,\
                           const strbuf m_name, const size_t binsize)
{
    latan_errno status;
    strbuf *field,latan_path;
    mat **m_prebin;
    size_t nf,lc,nmat;
    size_t i;
    
    status = LATAN_SUCCESS;
    i      = 0;
    nmat   = get_nfile(man_fname);
    field  = NULL;

    m_prebin = (m) ? mat_ar_create_from_dim(nmat,m[0]) : NULL;
    
    BEGIN_FOR_LINE_TOK(field,man_fname," \t",nf,lc)
    {
        if (field[0][0] != '#')
        {
            sprintf(latan_path,"%s%c%s",field[0],LATAN_PATH_SEP,m_name);
            if (m)
            {
                USTAT(mat_load(m_prebin[i],dim,latan_path));
            }
            else
            {
                USTAT(mat_load(NULL,dim,latan_path));
                break;
            }
            i++;
        }
    }
    END_FOR_LINE_TOK(field)
    if (m)
    {
        USTAT(mat_ar_bin(m,m_prebin,nmat,binsize));
        mat_ar_destroy(m_prebin,nmat);
    }
    
    return status;
}

/*                      random generator state I/O                          */
/****************************************************************************/
latan_errno randgen_save_state(const strbuf latan_path, const char mode,\
                               const rg_state state)
{
    latan_errno status;
    strbuf fname,elname;
    
    FUNC_INIT(fname,elname);
    if (strlen(elname) == 0)
    {
        LATAN_ERROR("no name specified",LATAN_EINVAL);
    }
    status = randgen_save_state_pt(fname,mode,state,elname);
    
    return status;
}

latan_errno randgen_load_state(rg_state state, const strbuf latan_path)
{
    latan_errno status;
    strbuf fname,elname;
    
    FUNC_INIT(fname,elname);
    status = randgen_load_state_pt(state,fname,elname);
    
    return status;
}

/*                          resampled sample I/O                            */
/****************************************************************************/
latan_errno rs_sample_save(const strbuf latan_path, const char mode,\
                           const rs_sample *s)
{
    latan_errno status;
    strbuf fname,elname;
    
    FUNC_INIT(fname,elname);
    if (strlen(elname) == 0)
    {
        LATAN_ERROR("no name specified",LATAN_EINVAL);
    }
    status = rs_sample_save_pt(fname,mode,s,elname);
    
    return status;
}

latan_errno rs_sample_save_subsamp(const strbuf latan_path, const char mode,\
                                   const rs_sample *s, const size_t k1,     \
                                   const size_t l1, const size_t k2,        \
                                   const size_t l2)
{
    latan_errno status;
    rs_sample *little_s;
    
    status = LATAN_SUCCESS;
    
    little_s = rs_sample_create(k2-k1+1,l2-l1+1,rs_sample_get_nsample(s));

    USTAT(rs_sample_get_subsamp(little_s,s,k1,l1,k2,l2));
    USTAT(rs_sample_save(latan_path,mode,little_s));
    
    rs_sample_destroy(little_s);
    
    return status;
}

latan_errno rs_sample_load(rs_sample *s, size_t *nsample, size_t *dim,\
                           const strbuf latan_path)
{
    latan_errno status;
    strbuf fname,elname;
    
    FUNC_INIT(fname,elname);
    status = rs_sample_load_pt(s,nsample,dim,fname,elname);
    
    return status;
}

latan_errno rs_sample_load_subsamp(rs_sample *s, const strbuf latan_path,\
                                   const size_t k1, const size_t l1,     \
                                   const size_t k2, const size_t l2)
{
    latan_errno status;
    rs_sample *big_s;
    size_t dim[2],nsample;
    
    status = LATAN_SUCCESS;
    
    USTAT(rs_sample_load(NULL,&nsample,dim,latan_path));
    big_s = rs_sample_create(dim[0],dim[1],nsample);
    USTAT(rs_sample_load(big_s,NULL,NULL,latan_path));
    USTAT(rs_sample_get_subsamp(s,big_s,k1,l1,k2,l2));
    
    rs_sample_destroy(big_s);
    
    return status;
}
