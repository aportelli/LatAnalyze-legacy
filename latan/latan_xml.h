/* latan_xml.h, part of LatAnalyze library
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

#ifndef LATAN_XML_H_
#define	LATAN_XML_H_

#include <latan/latan_globals.h>
#include <latan/latan_hadron.h>
#include <latan/latan_statistics.h>
#include <libxml/parser.h>
#include <libxml/threads.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>

#ifndef LIBXML_XPATH_ENABLED
#error libxml2 was built without XPath interface support
#endif

/* LatAnalyze XML format specification */
#ifndef LATAN_XML_VER
#define LATAN_XML_VER "1.0"
#endif
#ifndef LATAN_XML_ENC
#define LATAN_XML_ENC "UTF-8"
#endif
#ifndef LATAN_XMLNS
#define LATAN_XMLNS "latanalyze/xmlfmt/1.0"
#endif
#ifndef LATAN_XMLNS_PREF
#define LATAN_XMLNS_PREF "latan"
#endif

#define NXML_MARK 10
enum
{
    i_main    = 0,
    i_int     = 1,
    i_double  = 2,
    i_string  = 3,
    i_vect    = 4,
    i_mat     = 5,
    i_prop    = 6,
    i_rgstate = 7,
    i_sample  = 8
};

extern const strbuf xml_mark[NXML_MARK];

/* libxml2 workspace structure */
typedef struct
{
    xmlDoc *doc;
    xmlNs *ns;
    xmlNode *root;
    xmlXPathContext* ctxt;
} xml_workspace;

/* macro to check XML extension */
#define CHECK_XML_EXTENSION(f_name)\
{\
    char* ext;\
    ext = strrchr(f_name,'.');\
    if (ext == NULL)\
    {\
        sprintf(f_name,"%s.xml",f_name);\
    }\
    else if (strcmp(ext+1,"xml") != 0)\
    {\
        sprintf(f_name,"%s.xml",f_name);\
    }\
}

/* elementary tree browsing macros */
#define SKIP_COMMENTS(node)\
{\
    while (strcmp((const char *)(node)->name,"comment") == 0)\
    {\
        node = (node)->next;\
    }\
}

#define GOTO_LAST(node)\
{\
    while ((node)->next != NULL)\
    {\
        node = (node)->next;\
    }\
}

/* namespace/mark test macros */
#define GOT_LATAN_NS(node)\
(((node)->ns) &&\
 (strcmp((const char *)((node)->ns->href),LATAN_XMLNS) == 0) &&\
 (strcmp((const char *)((node)->ns->prefix),LATAN_XMLNS_PREF) == 0))

#define IF_GOT_LATAN_MARK_ELSE_ERROR(node,ind)\
SKIP_COMMENTS(node);\
if (!GOT_LATAN_NS(node))\
{\
    strbuf _errmsg;\
    sprintf(_errmsg,"XML namespace mismatch, expecting xmlns:%s=\"%s\" (%s:%d)",\
            LATAN_XMLNS_PREF,LATAN_XMLNS,(node)->doc->URL,(node)->line);\
    LATAN_ERROR(_errmsg,LATAN_ELATSYN);\
}\
else if (strcmp((const char *)((node)->name),xml_mark[ind]) != 0)\
{\
    strbuf _errmsg;\
    sprintf(_errmsg,"XML mark %s:%s not found (%s:%d)",LATAN_XMLNS_PREF,\
            xml_mark[i_int],(node)->doc->URL,(node)->line);\
    LATAN_ERROR(_errmsg,LATAN_ELATSYN);\
}\
else

/* file I/O macros */
#define BEGIN_XML_PARSING(ws,fname)\
{\
    MALLOC(ws,xml_workspace *,1);\
    (ws)->doc  = NULL;\
    (ws)->root = NULL;\
    (ws)->ns   = NULL;\
    (ws)->ctxt = NULL;\
    LIBXML_TEST_VERSION\
    (ws)->doc = xmlReadFile(fname,NULL,XML_PARSE_NOBLANKS|XML_PARSE_NONET);\
    if ((ws)->doc == NULL)\
    {\
        strbuf _errmsg;\
        sprintf(_errmsg,"impossible to parse file %s",fname);\
        LATAN_ERROR(_errmsg,LATAN_EFAULT);\
    }\
    (ws)->ctxt = xmlXPathNewContext((ws)->doc);\
    if ((ws)->ctxt == NULL)\
    {\
        strbuf _errmsg;\
        sprintf(_errmsg,"impossible to create XPath context (%s)",fname);\
        LATAN_ERROR(_errmsg,LATAN_EFAULT);\
    }\
    xmlXPathRegisterNs((ws)->ctxt,(const xmlChar *)LATAN_XMLNS_PREF,\
                       (const xmlChar *)LATAN_XMLNS);\
    (ws)->root = xmlDocGetRootElement((ws)->doc);\
    (ws)->ns   = (ws)->root->ns;\
    IF_GOT_LATAN_MARK_ELSE_ERROR((ws)->root,i_main)

#define END_XML_PARSING(ws,fname)\
    xmlFreeDoc((ws)->doc);\
    xmlXPathFreeContext((ws)->ctxt);\
    (ws)->doc  = NULL;\
    (ws)->root = NULL;\
    (ws)->ns   = NULL;\
    (ws)->ctxt = NULL;\
    FREE(ws);\
}

#define BEGIN_XML_WRITING_NEW(ws,fname)\
{\
    strbuf _fname;\
    strbuf _buf,_ver,_name;\
    MALLOC(ws,xml_workspace *,1);\
    (ws)->doc  = NULL;\
    (ws)->root = NULL;\
    (ws)->ns   = NULL;\
    (ws)->ctxt = NULL;\
    LIBXML_TEST_VERSION\
    (ws)->doc  = xmlNewDoc((const xmlChar *)LATAN_XML_VER);\
    (ws)->ns   = xmlNewNs((ws)->root,(const xmlChar *)LATAN_XMLNS,\
                          (const xmlChar *)LATAN_XMLNS_PREF);\
    (ws)->root = xmlNewNode((ws)->ns,(const xmlChar *)xml_mark[i_main]);\
    xmlDocSetRootElement((ws)->doc,(ws)->root);\
    latan_get_name(_name);\
    latan_get_version(_ver);\
    sprintf(_buf," XML file generated by %s v%s ",_name,_ver);\
    xmlAddPrevSibling((ws)->root,xmlNewComment((const xmlChar *)_buf));\

#define END_XML_WRITING_NEW(ws,fname)\
    strcpy(_fname,fname);\
    CHECK_XML_EXTENSION(_fname);\
    xml_save(ws,_fname);\
    xmlFreeNs((ws)->ns);\
    xmlFreeDoc((ws)->doc);\
    (ws)->doc  = NULL;\
    (ws)->root = NULL;\
    (ws)->ns   = NULL;\
    (ws)->ctxt = NULL;\
    FREE(ws);\
}

#define BEGIN_XML_WRITING_APPEND(ws,fname)\
{\
    BEGIN_XML_PARSING(ws,fname)\
    {

#define END_XML_WRITING_APPEND(ws,fname)\
        xml_save(ws,fname);\
    }\
    END_XML_PARSING(ws,fname)\
}

#define XML_WRITE(ws,fname,mode,instruction)\
switch(mode)\
{\
    case 'w':\
        BEGIN_XML_WRITING_NEW(ws,fname)\
        {\
            instruction;\
        }\
        END_XML_WRITING_NEW(ws,fname)\
        break;\
    case 'a':\
        BEGIN_XML_WRITING_APPEND(ws,fname)\
        {\
            instruction;\
        }\
        END_XML_WRITING_APPEND(ws,fname)\
        break;\
    default:\
        LATAN_ERROR("unknown file mode",LATAN_EINVAL);\
        break;\
}

__BEGIN_DECLS

/* data I/O functions */
/** input **/
latan_errno xml_get_int(int *res, xmlNode *node);
latan_errno xml_get_double(double *res, xmlNode *node);
latan_errno xml_get_string(strbuf res, xmlNode *node);
latan_errno xml_get_vect(mat *v, xmlNode *node);
latan_errno xml_get_vect_size(size_t *row, xmlNode *node);
latan_errno xml_get_mat(mat *m, xmlNode *node);
latan_errno xml_get_mat_size(size_t s[2], xmlNode *node);
latan_errno xml_get_prop(mat *prop, xmlNode *node);
latan_errno xml_get_prop_nt(size_t *nt, xmlNode *node);
latan_errno xml_get_rgstate(rg_state state, xmlNode *node);
latan_errno xml_get_sample(rs_sample *s, xmlNode *node);
latan_errno xml_get_sample_nsample(size_t *nsample, xmlNode *node);
latan_errno xml_get_sample_nrow(size_t *nr, xmlNode *node);

/** output **/
xmlNode * xml_insert_int(xmlNode *parent, const int res, const strbuf name);
xmlNode * xml_insert_double(xmlNode *parent, const double d, const strbuf name);
xmlNode * xml_insert_string(xmlNode *parent, const strbuf res,\
                            const strbuf name);
xmlNode * xml_insert_vect(xmlNode *parent, mat *v, const strbuf name);
xmlNode * xml_insert_mat(xmlNode *parent, mat *m, const strbuf name);
xmlNode * xml_insert_prop(xmlNode *parent, mat *prop,            \
                          const channel_no ch, const quark_no q1,\
                          const quark_no q2, const ss_no source, \
                          const ss_no sink, const strbuf name);
xmlNode * xml_insert_rgstate(xmlNode *parent, const rg_state state,\
                             const strbuf name);
xmlNode * xml_insert_sample(xmlNode *parent, const rs_sample *s,\
                            const strbuf name);

/* file writing function */
latan_errno xml_save(xml_workspace *ws, const strbuf fname);

/* XPath search function */
xmlXPathObject * xml_get_nodeset(strbuf xpath_expr, xml_workspace *ws);

__END_DECLS

#endif