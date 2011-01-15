/* latan_xml.c, part of LatAnalyze library
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

#include <latan/latan_xml.h>
#include <latan/latan_includes.h>

/*                   LatAnalyze XML format specification                    */
/****************************************************************************/
const strbuf xml_mark[NXML_MARK] =
{
    "data",   \
    "int",    \
    "double", \
    "string", \
    "vect",   \
    "mat",    \
    "prop",   \
    "rgstate",\
    "sample"  \
};

/*                          data I/O functions                              */
/****************************************************************************/
/* input */
latan_errno xml_get_int(int *res, xmlNode *node)
{
    char *buf;

    IF_GOT_LATAN_MARK_ELSE_ERROR(node,i_int)
    {
        buf = (char *)xmlNodeListGetString(node->doc,node->children,1);
        sscanf(buf,"%d",res);
        xmlFree(buf);
    }

    return LATAN_SUCCESS;
}

latan_errno xml_get_double(double *res, xmlNode *node)
{
    char *buf;

    IF_GOT_LATAN_MARK_ELSE_ERROR(node,i_double)
    {
        buf = (char *)xmlNodeListGetString(node->doc,node->children,1);
        sscanf(buf,"%lf",res);
        xmlFree(buf);
    }

    return LATAN_SUCCESS;
}

latan_errno xml_get_string(strbuf res, xmlNode *node)
{
    char *buf;

    IF_GOT_LATAN_MARK_ELSE_ERROR(node,i_int)
    {
        buf = (char *)xmlNodeListGetString(node->doc,node->children,1);
        strcpy(res,buf);
        xmlFree(buf);
    }

    return LATAN_SUCCESS;
}

latan_errno xml_get_vect(mat *v, xmlNode *node)
{
    xmlNode *vcur;
    double buf;
    size_t row;
    latan_errno status;

    row    = 0;
    status = LATAN_SUCCESS;

    IF_GOT_LATAN_MARK_ELSE_ERROR(node,i_vect)
    {
        for (vcur=node->children;vcur!=NULL;vcur=vcur->next)
        {
            LATAN_UPDATE_STATUS(status,xml_get_double(&buf,vcur));
            mat_set(v,row,0,buf);
            row++;
        }
    }

    return status;
}

latan_errno xml_get_vect_size(size_t *row, xmlNode *node)
{
    xmlNode *vcur;

    *row = 0;

    IF_GOT_LATAN_MARK_ELSE_ERROR(node,i_vect)
    {
        for (vcur=node->children;vcur!=NULL;vcur=vcur->next)
        {
            IF_GOT_LATAN_MARK_ELSE_ERROR(vcur,i_double)
            {
                (*row)++;
            }
        }
    }

    return LATAN_SUCCESS;
}

latan_errno xml_get_mat(mat *m, xmlNode *node)
{
    xmlNode *ccur;
    size_t j;
    mat mcol;
    latan_errno status;

    j             = 0;
    status        = LATAN_SUCCESS;
    mcol.mem_flag = CPU_LAST;

    IF_GOT_LATAN_MARK_ELSE_ERROR(node,i_mat)
    {
        for (ccur=node->children;ccur!=NULL;ccur=ccur->next)
        {
            MAT_PT_SUBM(&mcol,m,0,j,nrow(m)-1,j);
            LATAN_UPDATE_STATUS(status,xml_get_vect(&mcol,ccur));
            j++;
        }
    }

    return LATAN_SUCCESS;
}

latan_errno xml_get_mat_size(size_t s[2], xmlNode *node)
{
    xmlNode *ccur;
    size_t buf;
    latan_errno status;

    s[0]   = 0;
    s[1]   = 0;
    status = LATAN_SUCCESS;

    IF_GOT_LATAN_MARK_ELSE_ERROR(node,i_mat)
    {
        for (ccur=node->children;ccur!=NULL;ccur=ccur->next)
        {
            LATAN_UPDATE_STATUS(status,xml_get_vect_size(&buf,ccur));
            if ((s[0] != 0)&&(buf != s[0]))
            {
                strbuf errmsg;
                sprintf(errmsg,"reading matrix with variable column dimension (%s:%u)",\
                        node->doc->URL,node->line);
                LATAN_ERROR(errmsg,LATAN_EBADLEN);
            }
            s[0] = buf;
            s[1]++;
        }
    }

    return LATAN_SUCCESS;
}

latan_errno xml_get_prop(mat *prop, xmlNode *node)
{
    latan_errno status;

    IF_GOT_LATAN_MARK_ELSE_ERROR(node,i_prop)
    {
        status = xml_get_vect(prop,node->children);
    }

    return status;
}

latan_errno xml_get_prop_nt(size_t *nt, xmlNode *node)
{
    latan_errno status;

    IF_GOT_LATAN_MARK_ELSE_ERROR(node,i_prop)
    {
        status = xml_get_vect_size(nt,node->children);
    }

    return status;
}

latan_errno xml_get_rgstate(rg_state state, xmlNode *node)
{
    xmlNode *vcur;
    int i;
    latan_errno status;

    i      = 0;
    status = LATAN_SUCCESS;

    IF_GOT_LATAN_MARK_ELSE_ERROR(node,i_rgstate)
    {
        for (vcur=node->children;vcur!=NULL;vcur=vcur->next)
        {
            LATAN_UPDATE_STATUS(status,xml_get_int(state+i,vcur));
            i++;
        }
    }

    return status;
}

latan_errno xml_get_sample(rs_sample *s, xmlNode *node)
{
    xmlNode *scur;
    size_t i;
    latan_errno status;
    mat *pt;

    i      = 0;
    status = LATAN_SUCCESS;

    IF_GOT_LATAN_MARK_ELSE_ERROR(node,i_sample)
    {
        scur = node->children;
        pt   = rs_sample_pt_cent_val(s);
        LATAN_UPDATE_STATUS(status,xml_get_vect(pt,scur));
        for (scur=scur->next;scur!=NULL;scur=scur->next)
        {
            pt = rs_sample_pt_sample(s,i);
            LATAN_UPDATE_STATUS(status,xml_get_vect(pt,scur));
            i++;
        }
    }

    return LATAN_SUCCESS;
}

latan_errno xml_get_sample_nsample(size_t *nsample, xmlNode *node)
{
    xmlNode *scur;

    *nsample = 0;

    IF_GOT_LATAN_MARK_ELSE_ERROR(node,i_sample)
    {
        for (scur=node->children;scur!=NULL;scur=scur->next)
        {
            IF_GOT_LATAN_MARK_ELSE_ERROR(scur,i_vect)
            {
                (*nsample)++;
            }
        }
    }
    (*nsample)--;

    return LATAN_SUCCESS;
}

latan_errno xml_get_sample_nrow(size_t *nr, xmlNode *node)
{
    xmlNode *scur;
    size_t buf;

    buf = 0;

    IF_GOT_LATAN_MARK_ELSE_ERROR(node,i_sample)
    {
        for (scur=node->children;scur!=NULL;scur=scur->next)
        {
            xml_get_vect_size(&buf,scur);
            if ((*nr != 0)&&(*nr != buf))
            {
                strbuf errmsg;
                sprintf(errmsg,"reading samples with variable length (%s:%u)",\
                        node->doc->URL,node->line);
                LATAN_ERROR(errmsg,LATAN_EBADLEN);
            }
            *nr = buf;
        }
    }

    return LATAN_SUCCESS;
}

/* output */
xmlNode * xml_insert_int(xmlNode *parent, const int i, const strbuf name)
{
    strbuf data;
    xmlNode *node_new;

    sprintf(data,"%d",i);
    node_new = xmlNewChild(parent,NULL,(const xmlChar *)xml_mark[i_int],\
                           (const xmlChar *)data);
    if (strlen(name) > 0)
    {
        xmlNewProp(node_new,(const xmlChar *)"name",(const xmlChar *)name);
    }

    return node_new;
}

xmlNode * xml_insert_double(xmlNode *parent, const double d, const strbuf name)
{
    strbuf data;
    xmlNode *node_new;

    sprintf(data,"%e",d);
    node_new = xmlNewChild(parent,NULL,(const xmlChar *)xml_mark[i_double],\
                           (const xmlChar *)data);
    if (strlen(name) > 0)
    {
        xmlNewProp(node_new,(const xmlChar *)"name",(const xmlChar *)name);
    }

    return node_new;
}

xmlNode * xml_insert_string(xmlNode *parent, const strbuf str,\
                            const strbuf name)
{
    xmlNode *node_new;

    node_new = xmlNewChild(parent,NULL,(const xmlChar *)xml_mark[i_string],\
                           (const xmlChar *)str);
    if (strlen(name) > 0)
    {
        xmlNewProp(node_new,(const xmlChar *)"name",(const xmlChar *)name);
    }

    return node_new;
}

xmlNode * xml_insert_vect(xmlNode *parent, mat *v, const strbuf name)
{
    xmlNode *node_new;
    size_t i;

    node_new = xmlNewChild(parent,NULL,(const xmlChar *)xml_mark[i_vect],\
                           (const xmlChar *)"");
    for (i=0;i<nrow(v);i++)
    {
        xml_insert_double(node_new,mat_get(v,i,0),"");
    }
    if (strlen(name) > 0)
    {
        xmlNewProp(node_new,(const xmlChar *)"name",(const xmlChar *)name);
    }

    return node_new;
}

xmlNode * xml_insert_mat(xmlNode *parent, mat *m, const strbuf name)
{
    xmlNode *node_new;
    mat mcol;
    strbuf buf;
    size_t j;

    mcol.mem_flag = CPU_LAST;
    
    node_new = xmlNewChild(parent,NULL,(const xmlChar *)xml_mark[i_mat],\
                           (const xmlChar *)"");
    for (j=0;j<ncol(m);j++)
    {
        MAT_PT_SUBM(&mcol,m,0,j,nrow(m)-1,j);
        sprintf(buf,"col%lu",(unsigned long)j);
        xml_insert_vect(node_new,&mcol,buf);
    }
    if (strlen(name) > 0)
    {
        xmlNewProp(node_new,(const xmlChar *)"name",(const xmlChar *)name);
    }

    return node_new;
}

xmlNode * xml_insert_prop(xmlNode *parent, mat *prop,            \
                          const channel_no ch, const quark_no q1,\
                          const quark_no q2, const ss_no source, \
                          const ss_no sink, const strbuf name)
{
    xmlNode *node_new;
    strbuf ch_id,q1_id,q2_id,source_id,sink_id;

    
    node_new = xmlNewChild(parent,NULL,(const xmlChar *)xml_mark[i_prop],\
                           (const xmlChar *)"");
    xml_insert_vect(node_new,prop,"");
    channel_id_get(ch_id,ch);
    quark_id_get(q1_id,q1);
    quark_id_get(q2_id,q2);
    ss_id_get(source_id,source);
    ss_id_get(sink_id,sink);
    if (strlen(name) > 0)
    {
        xmlNewProp(node_new,(const xmlChar *)"name",(const xmlChar *)name);
    }
    xmlNewProp(node_new,(const xmlChar *)"channel",(const xmlChar *)ch_id);
    xmlNewProp(node_new,(const xmlChar *)"mass1",(const xmlChar *)q1_id);
    xmlNewProp(node_new,(const xmlChar *)"mass2",(const xmlChar *)q2_id);
    xmlNewProp(node_new,(const xmlChar *)"source",(const xmlChar *)source_id);
    xmlNewProp(node_new,(const xmlChar *)"sink",(const xmlChar *)sink_id);

    return node_new;
}

xmlNode * xml_insert_rgstate(xmlNode *parent, const rg_state state,\
                             const strbuf name)
{
    xmlNode *node_new;
    int i;

    node_new = xmlNewChild(parent,NULL,(const xmlChar *)xml_mark[i_rgstate],\
                           (const xmlChar *)"");
    for (i=0;i<RLXG_STATE_SIZE;i++)
    {
        xml_insert_int(node_new,state[i],"");
    }
    if (strlen(name) > 0)
    {
        xmlNewProp(node_new,(const xmlChar *)"name",(const xmlChar *)name);
    }

    return node_new;
}

xmlNode * xml_insert_sample(xmlNode *parent, const rs_sample *s,\
                            const strbuf name)
{
    xmlNode *node_new;
    strbuf buf;
    size_t i,nsample;

    nsample = rs_sample_get_nsample(s);

    node_new = xmlNewChild(parent,NULL,(const xmlChar *)xml_mark[i_sample],\
                           (const xmlChar *)"");
    xml_insert_vect(node_new,rs_sample_pt_cent_val(s),"central");
    for (i=0;i<nsample;i++)
    {
        sprintf(buf,"sample%lu",(unsigned long)i);
        xml_insert_vect(node_new,rs_sample_pt_sample(s,i),buf);
    }
    if (strlen(name) > 0)
    {
        xmlNewProp(node_new,(const xmlChar *)"name",(const xmlChar *)name);
    }

    return node_new;
}

/*                          file writing function                           */
/****************************************************************************/
void xml_check_extension(strbuf f_name)
{
    char* ext = NULL;
    strbuf f_name_buf;

    ext = strrchr(f_name,'.');
    if (ext == NULL)
    {
        sprintf(f_name_buf,"%s.xml",f_name);
    }
    else if (strcmp(ext+1,"xml") != 0)
    {
        sprintf(f_name_buf,"%s.xml",f_name);
    }
    else
    {
        strcpy(f_name_buf,f_name);
    }
    
    strcpy(f_name,f_name_buf);
}

latan_errno xml_save(xml_workspace *ws, const strbuf fname)
{
    xmlDoc *doc;
    int blank_bak,size;
    char *buf;

    xmlReconciliateNs(ws->doc,ws->root);
    blank_bak = xmlKeepBlanksDefault(0);
    xmlDocDumpFormatMemoryEnc(ws->doc,(xmlChar **)(&buf),&size,LATAN_XML_ENC,0);
    doc = xmlReadMemory(buf,size,NULL,LATAN_XML_ENC,XML_PARSE_NOBLANKS);
    xmlSaveFormatFileEnc(fname,doc,LATAN_XML_ENC,1);
    xmlKeepBlanksDefault(blank_bak);

    xmlFreeDoc(doc);
    xmlFree(buf);

    return LATAN_SUCCESS;
}

/*                          XPath search function                           */
/****************************************************************************/
xmlXPathObject * xml_get_nodeset(strbuf xpath_expr, xml_workspace *ws)
{
    strbuf errmsg;
    xmlXPathObject* nodeset;

    nodeset = xmlXPathEvalExpression((const xmlChar *)xpath_expr,ws->ctxt);
    if (nodeset == NULL)
    {
        sprintf(errmsg,"impossible to evaluate XPath expression \"%s\" (%s)",\
                xpath_expr,ws->doc->URL);
        LATAN_ERROR_VAL(errmsg,LATAN_FAILURE,NULL);
    }

    return nodeset;
}
