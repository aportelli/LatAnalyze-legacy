/* latan_io_xml.c, part of LatAnalyze library
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
        if ((strcmp(FILE_BUF(thread)->fname,fname) != 0)||(mode == 'w')||\
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
static bool io_is_init = false;

void io_init_xml(void)
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
        xmlInitParser();
        io_is_init = true;
    }
}

void io_finish_xml(void)
{
    int i;
    
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
                xml_close_file(FILE_BUF(i));
                env.file_is_loaded[i] = false;
            }
        }
        xmlCleanupParser();
        io_is_init = false;
    }
}

/*                             mat I/O                                      */
/****************************************************************************/
latan_errno mat_save_xml(const strbuf fname, const char mode, const mat *m,\
                         const strbuf name)
{
    int thread;
    
#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif
    
    xml_open_file_buf(fname,mode);
    xml_insert_mat(FILE_BUF(thread)->root,m,name);
    
    return LATAN_SUCCESS;
}

latan_errno mat_load_dim_xml(size_t dim[2], const strbuf fname,\
                             const strbuf name)
{
    xmlXPathObject *nodeset;
    strbuf xpath_expr;
    latan_errno status;
    int thread;
    
    status = LATAN_SUCCESS;
    
#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif
    
    USTAT(xml_open_file_buf(fname,'r'));
    sprintf(xpath_expr,"/%s:%s/%s:%s",LATAN_XMLNS_PREF,xml_mark[i_main],\
            LATAN_XMLNS_PREF,xml_mark[i_mat]);
    if (strlen(name) > 0)
    {
        sprintf(xpath_expr,"%s[@name='%s']",xpath_expr,name);
    }
    nodeset = xml_get_nodeset(xpath_expr,FILE_BUF(thread));
    if (nodeset->nodesetval != NULL)
    {
        USTAT(xml_get_mat_size(dim,nodeset->nodesetval->nodeTab[0]));
    }
    else
    {
        strbuf errmsg;
        if (strlen(name) > 0)
        {
            sprintf(errmsg,"(name=\"%s\") ",name);
        }
        else
        {
            strbufcpy(errmsg,"");
        }
        sprintf(errmsg,"matrix %snot found in file %s",errmsg,fname);
        LATAN_ERROR(errmsg,LATAN_ELATSYN);
    }
    
    xmlXPathFreeObject(nodeset);
    
    return status;
}

latan_errno mat_load_xml(mat *m, const strbuf fname, const strbuf name)
{
    xmlXPathObject *nodeset;
    strbuf xpath_expr;
    latan_errno status;
    int thread;
    
    status = LATAN_SUCCESS;
    
#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif
    
    USTAT(xml_open_file_buf(fname,'r'));
    sprintf(xpath_expr,"/%s:%s/%s:%s",LATAN_XMLNS_PREF,xml_mark[i_main],\
            LATAN_XMLNS_PREF,xml_mark[i_mat]);
    if (strlen(name) > 0)
    {
        sprintf(xpath_expr,"%s[@name='%s']",xpath_expr,name);
    }
    nodeset = xml_get_nodeset(xpath_expr,FILE_BUF(thread));
    if (nodeset->nodesetval != NULL)
    {
        USTAT(xml_get_mat(m,nodeset->nodesetval->nodeTab[0]));
    }
    else
    {
        strbuf errmsg;
        if (strlen(name) > 0)
        {
            sprintf(errmsg,"(name=\"%s\") ",name);
        }
        else
        {
            strbufcpy(errmsg,"");
        }
        sprintf(errmsg,"matrix %snot found in file %s",errmsg,fname);
        LATAN_ERROR(errmsg,LATAN_ELATSYN);
    }
    
    xmlXPathFreeObject(nodeset);
    
    return status;
}

/*                          propagator I/O                                  */
/****************************************************************************/
latan_errno prop_load_nt_xml(size_t *nt, const channel_no channel,\
                             const quark_no q1, const quark_no q2,\
                             const ss_no source, const ss_no sink,\
                             strbuf fname)
{
    xmlXPathObject *nodeset;
    strbuf channel_id,q1_id,q2_id,source_id,sink_id,xpath_expr;
    latan_errno status;
    int thread;

#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif

    xml_open_file_buf(fname,'r');
    channel_id_get(channel_id,channel);
    quark_id_get(q1_id,q1);
    quark_id_get(q2_id,q2);
    ss_id_get(source_id,source);
    ss_id_get(sink_id,sink);
    sprintf(xpath_expr,"/%s:%s/%s:%s[@channel='%s' and @mass1='%s' and @mass2='%s' and @source='%s' and @sink='%s']",\
            LATAN_XMLNS_PREF,xml_mark[i_main],\
            LATAN_XMLNS_PREF,xml_mark[i_prop],\
            channel_id,q1_id,q2_id,source_id,sink_id);
    nodeset = xml_get_nodeset(xpath_expr,FILE_BUF(thread));
    if (nodeset->nodesetval != NULL)
    {
        status = xml_get_prop_nt(nt,nodeset->nodesetval->nodeTab[0]);
    }
    else
    {
        strbuf errmsg;
        sprintf(errmsg,"propagator (ch=\"%s\" m1=\"%s\" m2=\"%s\" so=\"%s\" si=\"%s\") not found in file %s",\
                channel_id,q1_id,q2_id,source_id,sink_id,fname);
        LATAN_ERROR(errmsg,LATAN_ELATSYN);
    }

    xmlXPathFreeObject(nodeset);

    return status;
}

latan_errno prop_load_xml(mat *prop, const channel_no channel, \
                          const quark_no q1, const quark_no q2,\
                          const ss_no source, const ss_no sink,\
                          strbuf fname)
{
    xmlXPathObject *nodeset;
    strbuf channel_id,q1_id,q2_id,source_id,sink_id,xpath_expr;
    latan_errno status;
    int thread;

#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif

    xml_open_file_buf(fname,'r');
    channel_id_get(channel_id,channel);
    quark_id_get(q1_id,q1);
    quark_id_get(q2_id,q2);
    ss_id_get(source_id,source);
    ss_id_get(sink_id,sink);
    sprintf(xpath_expr,"/%s:%s/%s:%s[@channel='%s' and @mass1='%s' and @mass2='%s' and @source='%s' and @sink='%s']",\
            LATAN_XMLNS_PREF,xml_mark[i_main],\
            LATAN_XMLNS_PREF,xml_mark[i_prop],\
            channel_id,q1_id,q2_id,source_id,sink_id);
    nodeset = xml_get_nodeset(xpath_expr,FILE_BUF(thread));
    if (nodeset->nodesetval != NULL)
    {
        status = xml_get_prop(prop,nodeset->nodesetval->nodeTab[0]);
    }
    else
    {
        strbuf errmsg;
        sprintf(errmsg,"propagator (ch=\"%s\" m1=\"%s\" m2=\"%s\" so=\"%s\" si=\"%s\") not found in file %s",\
                channel_id,q1_id,q2_id,source_id,sink_id,fname);
        LATAN_ERROR(errmsg,LATAN_ELATSYN);
    }

    xmlXPathFreeObject(nodeset);

    return status;
}

latan_errno prop_save_xml(strbuf fname, const char mode, mat *prop, \
                          const strbuf channel,                     \
                          const quark_no q1, const quark_no q2,     \
                          const ss_no source, const ss_no sink,     \
                          const strbuf name)
{
    int thread;

#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif

    xml_open_file_buf(fname,mode);
    xml_insert_prop(FILE_BUF(thread)->root,prop,channel,q1,q2,source,sink,\
                    name);

    return LATAN_SUCCESS;
}

/*                      random generator state I/O                          */
/****************************************************************************/
latan_errno randgen_save_state_xml(const strbuf fname, const char mode,\
                                   const rg_state state, const strbuf name)
{
    int thread;

#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif

    xml_open_file_buf(fname,mode);
    xml_insert_rgstate(FILE_BUF(thread)->root,state,name);

    return LATAN_SUCCESS;
}

latan_errno randgen_load_state_xml(rg_state state, const strbuf fname,\
                                   const strbuf name)
{
    xmlXPathObject *nodeset;
    strbuf xpath_expr;
    latan_errno status;
    int thread;

    status = LATAN_SUCCESS;
    
#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif

    USTAT(xml_open_file_buf(fname,'r'));
    sprintf(xpath_expr,"/%s:%s/%s:%s",LATAN_XMLNS_PREF,xml_mark[i_main],\
            LATAN_XMLNS_PREF,xml_mark[i_rgstate]);
    if (strlen(name) > 0)
    {
        sprintf(xpath_expr,"%s[@name='%s']",xpath_expr,name);
    }
    nodeset = xml_get_nodeset(xpath_expr,FILE_BUF(thread));
    if (nodeset->nodesetval != NULL)
    {
        USTAT(xml_get_rgstate(state,nodeset->nodesetval->nodeTab[0]));
    }
    else
    {
        strbuf errmsg;
        if (strlen(name) > 0)
        {
            sprintf(errmsg,"(name=\"%s\") ",name);
        }
        else
        {
            strbufcpy(errmsg,"");
        }
        sprintf(errmsg,"rgstate %snot found in file %s",errmsg,fname);
        LATAN_ERROR(errmsg,LATAN_ELATSYN);
    }

    xmlXPathFreeObject(nodeset);

    return status;
}

/*                          resampled sample I/O                            */
/****************************************************************************/
latan_errno rs_sample_save_xml(const strbuf fname, const char mode,\
                               const rs_sample *s)
{
    strbuf name;
    int thread;

#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif

    xml_open_file_buf(fname,mode);
    rs_sample_get_name(name,s);
    if (strlen(name) == 0)
    {
        LATAN_ERROR("cannot save sample with an empty name",LATAN_EINVAL);
    }
    xml_insert_sample(FILE_BUF(thread)->root,s,name);

    return LATAN_SUCCESS;
}

latan_errno rs_sample_load_nrow_xml(size_t *nr, const strbuf fname,\
                                    const strbuf name)
{
    xmlXPathObject *nodeset;
    strbuf xpath_expr;
    latan_errno status;
    int thread;
    
    status = LATAN_SUCCESS;

#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif

    USTAT(xml_open_file_buf(fname,'r'));
    sprintf(xpath_expr,"/%s:%s/%s:%s",LATAN_XMLNS_PREF,xml_mark[i_main],\
            LATAN_XMLNS_PREF,xml_mark[i_sample]);
    if (strlen(name) > 0)
    {
        sprintf(xpath_expr,"%s[@name='%s']",xpath_expr,name);
    }
    nodeset = xml_get_nodeset(xpath_expr,FILE_BUF(thread));
    if (nodeset->nodesetval != NULL)
    {
        USTAT(xml_get_sample_nrow(nr,nodeset->nodesetval->nodeTab[0]));
    }
    else
    {
        strbuf errmsg;
        if (strlen(name) > 0)
        {
            sprintf(errmsg,"(name=\"%s\") ",name);
        }
        else
        {
            strbufcpy(errmsg,"");
        }
        sprintf(errmsg,"sample %snot found in file %s",errmsg,fname);
        LATAN_ERROR(errmsg,LATAN_ELATSYN);
    }

    xmlXPathFreeObject(nodeset);

    return status;
}

latan_errno rs_sample_load_nsample_xml(size_t *nsample, const strbuf fname,\
                                       const strbuf name)
{
    xmlXPathObject *nodeset;
    strbuf xpath_expr;
    latan_errno status;
    int thread;
    
    status = LATAN_SUCCESS;

#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif

    USTAT(xml_open_file_buf(fname,'r'));
    sprintf(xpath_expr,"/%s:%s/%s:%s",LATAN_XMLNS_PREF,xml_mark[i_main],\
            LATAN_XMLNS_PREF,xml_mark[i_sample]);
    if (strlen(name) > 0)
    {
        sprintf(xpath_expr,"%s[@name='%s']",xpath_expr,name);
    }
    nodeset = xml_get_nodeset(xpath_expr,FILE_BUF(thread));
    if (nodeset->nodesetval != NULL)
    {
        USTAT(xml_get_sample_nsample(nsample,nodeset->nodesetval->nodeTab[0]));
    }
    else
    {
        strbuf errmsg;
        if (strlen(name) > 0)
        {
            sprintf(errmsg,"(name=\"%s\") ",name);
        }
        else
        {
            strbufcpy(errmsg,"");
        }
        sprintf(errmsg,"sample %snot found in file %s",errmsg,fname);
        LATAN_ERROR(errmsg,LATAN_ELATSYN);
    }

    xmlXPathFreeObject(nodeset);

    return status;
}

latan_errno rs_sample_load_xml(rs_sample *s, const strbuf fname,\
                               const strbuf name)
{
    xmlXPathObject *nodeset;
    strbuf xpath_expr;
    latan_errno status;
    int thread;

    status = LATAN_SUCCESS;
    
#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif

    USTAT(xml_open_file_buf(fname,'r'));
    sprintf(xpath_expr,"/%s:%s/%s:%s",LATAN_XMLNS_PREF,xml_mark[i_main],\
            LATAN_XMLNS_PREF,xml_mark[i_sample]);
    if (strlen(name) > 0)
    {
        sprintf(xpath_expr,"%s[@name='%s']",xpath_expr,name);
    }
    nodeset = xml_get_nodeset(xpath_expr,FILE_BUF(thread));
    if (nodeset->nodesetval != NULL)
    {
        USTAT(xml_get_sample(s,nodeset->nodesetval->nodeTab[0]));
    }
    else
    {
        strbuf errmsg;
        if (strlen(name) > 0)
        {
            sprintf(errmsg,"(name=\"%s\") ",name);
        }
        else
        {
            strbufcpy(errmsg,"");
        }
        sprintf(errmsg,"sample %snot found in file %s",errmsg,fname);
        LATAN_ERROR(errmsg,LATAN_ELATSYN);
    }

    xmlXPathFreeObject(nodeset);

    return LATAN_SUCCESS;
}

#endif
