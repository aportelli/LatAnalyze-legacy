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

#include <latan/latan_io_xml.h>
#include <latan/latan_includes.h>
#include <latan/latan_xml.h>

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
    if (io_is_init)
    {
#ifdef _OPENMP
        if(omp_in_parallel())
        {
            LATAN_WARNING("I/O finish called from a parallel region",\
                          LATAN_FAILURE);
        }
#endif
        xmlCleanupParser();
        io_is_init = false;
    }
}

/*                          propagator I/O                                  */
/****************************************************************************/
latan_errno prop_load_nt_xml(size_t *nt, const channel_no channel,\
                             const quark_no q1, const quark_no q2,\
                             const ss_no source, const ss_no sink,\
                             strbuf fname)
{
    xml_file *f;
    xmlXPathObject *nodeset;
    strbuf channel_id,q1_id,q2_id,source_id,sink_id,xpath_expr;
    latan_errno status;

    f = xml_open_file(fname,'r');

    channel_id_get(channel_id,channel);
    quark_id_get(q1_id,q1);
    quark_id_get(q2_id,q2);
    ss_id_get(source_id,source);
    ss_id_get(sink_id,sink);
    sprintf(xpath_expr,"/%s:%s/%s:%s[@channel='%s' and @mass1='%s' and @mass2='%s' and @source='%s' and @sink='%s']",\
            LATAN_XMLNS_PREF,xml_mark[i_main],\
            LATAN_XMLNS_PREF,xml_mark[i_prop],\
            channel_id,q1_id,q2_id,source_id,sink_id);
    nodeset = xml_get_nodeset(xpath_expr,f);
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
    xml_close_file(f);

    return status;
}

latan_errno prop_load_xml(mat *prop, const channel_no channel, \
                          const quark_no q1, const quark_no q2,\
                          const ss_no source, const ss_no sink,\
                          strbuf fname)
{
    xml_file *f;
    xmlXPathObject *nodeset;
    strbuf channel_id,q1_id,q2_id,source_id,sink_id,xpath_expr;
    latan_errno status;

    f = xml_open_file(fname,'r');

    channel_id_get(channel_id,channel);
    quark_id_get(q1_id,q1);
    quark_id_get(q2_id,q2);
    ss_id_get(source_id,source);
    ss_id_get(sink_id,sink);
    sprintf(xpath_expr,"/%s:%s/%s:%s[@channel='%s' and @mass1='%s' and @mass2='%s' and @source='%s' and @sink='%s']",\
            LATAN_XMLNS_PREF,xml_mark[i_main],\
            LATAN_XMLNS_PREF,xml_mark[i_prop],\
            channel_id,q1_id,q2_id,source_id,sink_id);
    nodeset = xml_get_nodeset(xpath_expr,f);
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
    xml_close_file(f);

    return status;
}

/*                      random generator state I/O                          */
/****************************************************************************/
latan_errno randgen_save_state_xml(const strbuf fname, const char mode,\
                                   const rg_state state, const strbuf name)
{
    xml_file *f;

    if (mode == 'w')
    {
        f = xml_new_file(fname);
    }
    else
    {
        f = xml_open_file(fname,mode);
    }

    xml_insert_rgstate(f->root,state,name);

    xml_close_file(f);

    return LATAN_SUCCESS;
}

latan_errno randgen_load_state_xml(rg_state state, const strbuf fname,\
                                   const strbuf name)
{
    xml_file *f;
    xmlXPathObject *nodeset;
    strbuf xpath_expr;
    latan_errno status;

    f = xml_open_file(fname,'r');

    sprintf(xpath_expr,"/%s:%s/%s:%s",LATAN_XMLNS_PREF,xml_mark[i_main],\
            LATAN_XMLNS_PREF,xml_mark[i_rgstate]);
    if (strlen(name) > 0)
    {
        sprintf(xpath_expr,"%s[@name='%s']",xpath_expr,name);
    }
    nodeset = xml_get_nodeset(xpath_expr,f);
    if (nodeset->nodesetval != NULL)
    {
        status = xml_get_rgstate(state,nodeset->nodesetval->nodeTab[0]);
    }
    else
    {
        strbuf errmsg;
        if (strlen(name) > 0)
        {
            sprintf(errmsg,"(name=\"%s\")",name);
        }
        else
        {
            strbufcpy(errmsg,"");
        }
        sprintf(errmsg,"rgstate %snot found in file %s",errmsg,fname);
        LATAN_ERROR(errmsg,LATAN_ELATSYN);
    }

    xmlXPathFreeObject(nodeset);
    xml_close_file(f);

    return LATAN_SUCCESS;
}

/*                          resampled sample I/O                            */
/****************************************************************************/
latan_errno rs_sample_save_xml(const strbuf fname, const char mode,\
                               const rs_sample *s)
{
    xml_file *f;
    strbuf name;

    if (mode == 'w')
    {
        f = xml_new_file(fname);
    }
    else
    {
        f = xml_open_file(fname,mode);
    }

    rs_sample_get_name(name,s);
    if (strlen(name) == 0)
    {
        LATAN_ERROR("cannot save sample with an empty name",LATAN_EINVAL);
    }
    xml_insert_sample(f->root,s,name);

    xml_close_file(f);

    return LATAN_SUCCESS;
}

latan_errno rs_sample_load_nrow_xml(size_t *nr, const strbuf fname,\
                                    const strbuf name)
{
    xml_file *f;
    xmlXPathObject *nodeset;
    strbuf xpath_expr;
    latan_errno status;

    f = xml_open_file(fname,'r');

    sprintf(xpath_expr,"/%s:%s/%s:%s",LATAN_XMLNS_PREF,xml_mark[i_main],\
            LATAN_XMLNS_PREF,xml_mark[i_sample]);
    if (strlen(name) > 0)
    {
        sprintf(xpath_expr,"%s[@name='%s']",xpath_expr,name);
    }
    nodeset = xml_get_nodeset(xpath_expr,f);
    if (nodeset->nodesetval != NULL)
    {
        status = xml_get_sample_nrow(nr,nodeset->nodesetval->nodeTab[0]);
    }
    else
    {
        strbuf errmsg;
        if (strlen(name) > 0)
        {
            sprintf(errmsg,"(name=\"%s\")",name);
        }
        else
        {
            strbufcpy(errmsg,"");
        }
        sprintf(errmsg,"sample %snot found in file %s",errmsg,fname);
        LATAN_ERROR(errmsg,LATAN_ELATSYN);
    }

    xmlXPathFreeObject(nodeset);
    xml_close_file(f);

    return LATAN_SUCCESS;
}

latan_errno rs_sample_load_nsample_xml(size_t *nsample, const strbuf fname,\
                                       const strbuf name)
{
    xml_file *f;
    xmlXPathObject *nodeset;
    strbuf xpath_expr;
    latan_errno status;

    f = xml_open_file(fname,'r');

    sprintf(xpath_expr,"/%s:%s/%s:%s",LATAN_XMLNS_PREF,xml_mark[i_main],\
            LATAN_XMLNS_PREF,xml_mark[i_sample]);
    if (strlen(name) > 0)
    {
        sprintf(xpath_expr,"%s[@name='%s']",xpath_expr,name);
    }
    nodeset = xml_get_nodeset(xpath_expr,f);
    if (nodeset->nodesetval != NULL)
    {
        status = xml_get_sample_nsample(nsample,\
                                        nodeset->nodesetval->nodeTab[0]);
    }
    else
    {
        strbuf errmsg;
        if (strlen(name) > 0)
        {
            sprintf(errmsg,"(name=\"%s\")",name);
        }
        else
        {
            strbufcpy(errmsg,"");
        }
        sprintf(errmsg,"sample %snot found in file %s",errmsg,fname);
        LATAN_ERROR(errmsg,LATAN_ELATSYN);
    }

    xmlXPathFreeObject(nodeset);
    xml_close_file(f);

    return LATAN_SUCCESS;
}

latan_errno rs_sample_load_xml(rs_sample *s, const strbuf fname,\
                               const strbuf name)
{
    xml_file *f;
    xmlXPathObject *nodeset;
    strbuf xpath_expr;
    latan_errno status;

    f = xml_open_file(fname,'r');

    sprintf(xpath_expr,"/%s:%s/%s:%s",LATAN_XMLNS_PREF,xml_mark[i_main],\
            LATAN_XMLNS_PREF,xml_mark[i_sample]);
    if (strlen(name) > 0)
    {
        sprintf(xpath_expr,"%s[@name='%s']",xpath_expr,name);
    }
    nodeset = xml_get_nodeset(xpath_expr,f);
    if (nodeset->nodesetval != NULL)
    {
        status = xml_get_sample(s,nodeset->nodesetval->nodeTab[0]);
    }
    else
    {
        strbuf errmsg;
        if (strlen(name) > 0)
        {
            sprintf(errmsg,"(name=\"%s\")",name);
        }
        else
        {
            strbufcpy(errmsg,"");
        }
        sprintf(errmsg,"sample %snot found in file %s",errmsg,fname);
        LATAN_ERROR(errmsg,LATAN_ELATSYN);
    }

    xmlXPathFreeObject(nodeset);
    xml_close_file(f);

    return LATAN_SUCCESS;
}
