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

#define NXML_MARK 9
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
    strbuf fname;
    char mode;
} xml_file;

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
xmlNode * xml_insert_mat(xmlNode *parent, const mat *m, const strbuf name);
xmlNode * xml_insert_prop(xmlNode *parent, mat *prop,                   \
                          const strbuf ch, const quark_no q1,           \
                          const quark_no q2, const ss_no source,        \
                          const ss_no sink, const strbuf name);
xmlNode * xml_insert_rgstate(xmlNode *parent, const rg_state state,\
                             const strbuf name);
xmlNode * xml_insert_sample(xmlNode *parent, const rs_sample *s,\
                            const strbuf name);

/* file writing function */
void xml_check_extension(strbuf fname);
xml_file * xml_new_file(const strbuf fname);
xml_file * xml_open_file(const strbuf fname, const char mode);
latan_errno xml_save_file(xml_file *f);
void xml_file_destroy(xml_file *f);
latan_errno xml_close_file(xml_file *f);

/* XPath search function */
xmlXPathObject * xml_get_nodeset(strbuf xpath_expr, xml_file *f);

__END_DECLS

#endif
