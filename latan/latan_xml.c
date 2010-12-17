#include <latan/latan_xml.h>
#include <latan/latan_includes.h>

const strbuf xml_mark[NXML_MARK] =
{
    "data",   \
    "int",    \
    "double", \
    "string", \
    "vect",   \
    "mat",    \
    "row",    \
    "prop",   \
    "rgstate",\
    "sample"  \
};

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
    xmlNode *rcur,*ccur;
    double buf;
    size_t i,j;
    latan_errno status;

    i      = 0;
    j      = 0;
    status = LATAN_SUCCESS;
    
    IF_GOT_LATAN_MARK_ELSE_ERROR(node,i_mat)
    {
        for (rcur=node->children;rcur!=NULL;rcur=rcur->next)
        {
            IF_GOT_LATAN_MARK_ELSE_ERROR(rcur,i_row)
            {
                for (ccur=rcur->children;ccur!=NULL;ccur=ccur->next)
                {
                    LATAN_UPDATE_STATUS(status,xml_get_double(&buf,ccur));
                    mat_set(m,i,j,buf);
                    j++;
                }
            }
            i++;
        }
    }

    return status;
}

latan_errno xml_get_mat_size(size_t s[2], xmlNode *node)
{
    xmlNode *rcur,*ccur;

    s[0] = 0;
    s[1] = 0;

    IF_GOT_LATAN_MARK_ELSE_ERROR(node,i_mat)
    {
        for (rcur=node->children;rcur!=NULL;rcur=rcur->next)
        {
            IF_GOT_LATAN_MARK_ELSE_ERROR(rcur,i_row)
            {
                for (ccur=rcur->children;ccur!=NULL;ccur=ccur->next)
                {
                    IF_GOT_LATAN_MARK_ELSE_ERROR(ccur,i_double)
                    {
                        s[1]++;
                    }
                }
            }
            s[0]++;
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
