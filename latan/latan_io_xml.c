/* latan_io_xml.c, part of LatAnalyze library
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


#include <latan/latan_includes.h>
#ifdef HAVE_LIBXML2
#include <latan/latan_io_xml.h>
#include <latan/latan_io.h>
#include <latan/latan_xml.h>

/*                       XML buffer management (internal)                   */
/****************************************************************************/
typedef struct
{
    xml_file **xml_buf;
    bool *file_is_loaded;
    int nfile;
} io_xml_env;

static io_xml_env env =
{
    NULL,\
    NULL,\
    0    \
};

#define FILE_BUF(thread) env.xml_buf[thread] /* type : xml_file * */

static latan_errno xml_open_file_buf(const strbuf fname, const char mode)
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
    status = LATAN_SUCCESS;

#ifdef _OPENMP
    #pragma omp critical
#endif
    {
        if (nthread > env.nfile)
        {
            REALLOC_NOERRET(env.xml_buf,env.xml_buf,xml_file **,nthread);
            REALLOC_NOERRET(env.file_is_loaded,env.file_is_loaded,bool *,\
                            nthread);
            for (i=env.nfile;i<nthread;i++)
            {
                env.file_is_loaded[i] = false;
                FILE_BUF(i)        = NULL;
            }
            env.nfile = nthread;
        }
    }
    if (env.file_is_loaded[thread])
    {
        if ((strbufcmp(FILE_BUF(thread)->fname,fname) != 0)||(mode == 'w')||\
            (mode != FILE_BUF(thread)->mode))
        {
            status           = xml_close_file(FILE_BUF(thread));
            FILE_BUF(thread) = xml_open_file(fname,mode);
        }
    }
    else
    {
        FILE_BUF(thread)           = xml_open_file(fname,mode);
        env.file_is_loaded[thread] = true;
    }

    return status;
}

/*                              I/O init/finish                             */
/****************************************************************************/
void io_init_xml(void)
{
#ifdef _OPENMP
    if(omp_in_parallel())
    {
        LATAN_WARNING("I/O initialization called from a parallel region",\
                      LATAN_FAILURE);
    }
#endif
    xmlInitParser();
}

void io_finish_xml(void)
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
            xml_close_file(FILE_BUF(i));
            env.file_is_loaded[i] = false;
        }
    }
    xmlCleanupParser();
}

/*                             mat I/O                                      */
/****************************************************************************/
latan_errno mat_save_xml(const strbuf fname, const char mode, const mat *m,\
                         const strbuf name)
{
    latan_errno status;
    int thread;
    
#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif
    
    if ((mode == 'w')||(mode == 'a'))
    {
        status = xml_open_file_buf(fname,mode);
    }
    else
    {
        LATAN_ERROR("unknown or read-only file mode",LATAN_EINVAL);
    }
    xml_insert_mat(FILE_BUF(thread)->root,m,name);
    
    return status;
}

latan_errno mat_load_xml(mat *m, size_t *dim, const strbuf fname,\
                         const strbuf name)
{
    xmlXPathObject *nodeset;
    strbuf xpath_expr,xpath_name;
    latan_errno status;
    int thread;
    
#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif
    status = LATAN_SUCCESS;
    
    USTAT(xml_open_file_buf(fname,'r'));
    if (strlen(name) == 0)
    {
        strbufcpy(xpath_name,"");
    }
    else
    {
        sprintf(xpath_name,"[@name='%s']",name);
    }
    sprintf(xpath_expr,"/%s:%s/%s:%s%s",LATAN_XMLNS_PREF,\
            xml_mark[i_main],LATAN_XMLNS_PREF,xml_mark[i_mat],xpath_name);
    nodeset = xml_get_nodeset(xpath_expr,FILE_BUF(thread));
    if (nodeset->nodesetval != NULL)
    {
        if (m)
        {
            USTAT(xml_get_mat(m,nodeset->nodesetval->nodeTab[0]));
        }
        if (dim) 
        {
            USTAT(xml_get_mat_size(dim,nodeset->nodesetval->nodeTab[0]));
        }
    }
    else
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
        LATAN_ERROR(errmsg,LATAN_ELATSYN);
    }
    
    xmlXPathFreeObject(nodeset);
    
    return status;
}

/*                      random generator state I/O                          */
/****************************************************************************/
latan_errno randgen_save_state_xml(const strbuf fname, const char mode,\
                                   const rg_state state, const strbuf name)
{
    latan_errno status;
    int thread;
    
#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif

    if ((mode == 'w')||(mode == 'a'))
    {
        status = xml_open_file_buf(fname,mode);
    }
    else
    {
        LATAN_ERROR("unknown or read-only file mode",LATAN_EINVAL);
    }
    xml_insert_rgstate(FILE_BUF(thread)->root,state,name);

    return status;
}

latan_errno randgen_load_state_xml(rg_state state, const strbuf fname,\
                                   const strbuf name)
{
    xmlXPathObject *nodeset;
    strbuf xpath_expr,xpath_name;
    latan_errno status;
    int thread;

#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif
    status = LATAN_SUCCESS;

    USTAT(xml_open_file_buf(fname,'r'));
    if (strlen(name) == 0)
    {
        strbufcpy(xpath_name,"");
    }
    else
    {
        sprintf(xpath_name,"[@name='%s']",name);
    }
    sprintf(xpath_expr,"/%s:%s/%s:%s%s",LATAN_XMLNS_PREF,\
            xml_mark[i_main],LATAN_XMLNS_PREF,xml_mark[i_rgstate],xpath_name);
    nodeset = xml_get_nodeset(xpath_expr,FILE_BUF(thread));
    if (nodeset->nodesetval != NULL)
    {
        USTAT(xml_get_rgstate(state,nodeset->nodesetval->nodeTab[0]));
    }
    else
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
        sprintf(errmsg,"rgstate (name= %s) not found in file %s",buf,fname);
        LATAN_ERROR(errmsg,LATAN_ELATSYN);
    }

    xmlXPathFreeObject(nodeset);

    return status;
}

/*                          resampled sample I/O                            */
/****************************************************************************/
latan_errno rs_sample_save_xml(const strbuf fname, const char mode,\
                               const rs_sample *s, const strbuf name)
{
    latan_errno status;
    int thread;
    
#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif

    if ((mode == 'w')||(mode == 'a'))
    {
        status = xml_open_file_buf(fname,mode);
    }
    else
    {
        LATAN_ERROR("unknown or read-only file mode",LATAN_EINVAL);
    }
    xml_insert_sample(FILE_BUF(thread)->root,s,name);

    return status;
}

latan_errno rs_sample_load_xml(rs_sample *s, size_t *nsample, size_t *dim,\
                               const strbuf fname, const strbuf name)
{
    xmlXPathObject *nodeset;
    strbuf xpath_expr,xpath_name;
    latan_errno status;
    int thread;
    
#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif
    status = LATAN_SUCCESS;

    USTAT(xml_open_file_buf(fname,'r'));
    if (strlen(name) == 0)
    {
        strbufcpy(xpath_name,"");
    }
    else
    {
        sprintf(xpath_name,"[@name='%s']",name);
    }
    sprintf(xpath_expr,"/%s:%s/%s:%s%s",LATAN_XMLNS_PREF,  \
            xml_mark[i_main],LATAN_XMLNS_PREF,xml_mark[i_sample],xpath_name);
    nodeset = xml_get_nodeset(xpath_expr,FILE_BUF(thread));
    if (nodeset->nodesetval != NULL)
    {
        if (s)
        {
            USTAT(xml_get_sample(s,nodeset->nodesetval->nodeTab[0]));
        }
        if (nsample) 
        {
            USTAT(xml_get_sample_nsample(nsample,                        \
                                         nodeset->nodesetval->nodeTab[0]));
        }
        if (dim)
        {
            USTAT(xml_get_sample_size(dim,nodeset->nodesetval->nodeTab[0]));
        }
    }
    else
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
        sprintf(errmsg,"sample (name= %s) not found in file %s",buf,fname);
        LATAN_ERROR(errmsg,LATAN_ELATSYN);
    }

    xmlXPathFreeObject(nodeset);

    return status;
}

#endif
