#include <latan/latan_io.h>
#include <latan/latan_includes.h>
#include <latan/latan_math.h>
#include <latan/latan_xml.h>

/*                              I/O init/finish                             */
/****************************************************************************/
static bool io_is_init = false;

void io_init(void)
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

void io_finish(void)
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


/*                                general I/O                               */
/****************************************************************************/
int get_nfile(const strbuf manifestfname)
{
    strbuf buf1, buf2;
    int nfile;
    FILE* manifest = NULL;
    
    nfile = 0;
    
    FOPEN_ERRVAL(manifest,manifestfname,"r",LATAN_FAILURE);
    while (!feof(manifest))
    {
        if ((fgets(buf1,STRING_LENGTH,manifest))&&(sscanf(buf1,"%s\n",buf2)>0))
        {
            nfile++;
        }
    }
    fclose(manifest);
    
    return nfile;
}

latan_errno get_firstfname(strbuf fname, const strbuf manifestfname)
{
    strbuf buf;
    FILE* manifest = NULL;
    
    FOPEN(manifest,manifestfname,"r");
    fgets(buf,STRING_LENGTH,manifest);
    fclose(manifest);
    sscanf(buf,"%s\n",fname);
    
    return LATAN_SUCCESS;
}

/*                              mat I/O                                     */
/****************************************************************************/
void mat_dump(FILE* stream, mat *m, const strbuf fmt)
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

latan_errno mat_save_plotdat(mat *x, mat *dat, const strbuf fname)
{
    FILE* f = NULL;
    size_t i;
    
    FOPEN(f,fname,"w");
    
    for (i=0;i<nrow(dat);i++)
    {
        fprintf(f,"%.10e\t%.10e\n",mat_get(x,i,0),mat_get(dat,i,0));
    }
    fclose(f);
    
    return LATAN_SUCCESS;
}

latan_errno mat_save_plotdat_yerr(mat *x, mat *dat, mat *yerr,\
                                  const strbuf fname)
{
    FILE* f = NULL;
    size_t i;
    
    FOPEN(f,fname,"w");
    for (i=0;i<nrow(dat);i++)
    {
        fprintf(f,"%.10e\t%.10e\t%.10e\n",mat_get(x,i,0),mat_get(dat,i,0),\
                mat_get(yerr,i,0));
    }
    fclose(f);
    
    return LATAN_SUCCESS;
}

latan_errno mat_save_plotdat_xyerr(mat *x, mat *dat, mat *xerr,\
                                   mat *yerr, const strbuf fname)
{
    FILE* f = NULL;
    size_t i;
    
    FOPEN(f,fname,"w");
    for (i=0;i<nrow(dat);i++)
    {
        fprintf(f,"%.10e\t%.10e\t%.10e\t%.10e\n",mat_get(x,i,0),mat_get(dat,i,0),\
                mat_get(xerr,i,0),mat_get(yerr,i,0));
    }
    fclose(f);
    
    return LATAN_SUCCESS;
}


/*                          propagator I/O                                  */
/****************************************************************************/
latan_errno prop_load(mat *prop, const channel_no channel, \
                      const quark_no q1, const quark_no q2,\
                      const ss_no source, const ss_no sink,\
                      strbuf fname)
{
    xml_workspace *ws;
    xmlXPathObject *nodeset;
    strbuf channel_id,q1_id,q2_id,source_id,sink_id,xpath_expr;
    latan_errno status;

    channel_id_get(channel_id,channel);
    quark_id_get(q1_id,q1);
    quark_id_get(q2_id,q2);
    ss_id_get(source_id,source);
    ss_id_get(sink_id,sink);
    BEGIN_XML_PARSING(ws,fname)
    {
        sprintf(xpath_expr,"/%s:%s/%s:%s[@channel='%s' and @mass1='%s' and @mass2='%s' and @source='%s' and @sink='%s']",\
                LATAN_XMLNS_PREF,xml_mark[i_main],\
                LATAN_XMLNS_PREF,xml_mark[i_prop],\
                channel_id,q1_id,q2_id,source_id,sink_id);
        nodeset = xml_get_nodeset(xpath_expr,ws);
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
    }
    END_XML_PARSING(ws,fname)
    
    return status;
}

latan_errno prop_load_nt(size_t *nt, const channel_no channel,\
                         const quark_no q1, const quark_no q2,\
                         const ss_no source, const ss_no sink,\
                         strbuf fname)
{
    xml_workspace *ws;
    xmlXPathObject *nodeset;
    strbuf channel_id,q1_id,q2_id,source_id,sink_id,xpath_expr;
    latan_errno status;

    channel_id_get(channel_id,channel);
    quark_id_get(q1_id,q1);
    quark_id_get(q2_id,q2);
    ss_id_get(source_id,source);
    ss_id_get(sink_id,sink);
    BEGIN_XML_PARSING(ws,fname)
    {
        sprintf(xpath_expr,"/%s:%s/%s:%s[@channel='%s' and @mass1='%s' and @mass2='%s' and @source='%s' and @sink='%s']",\
                LATAN_XMLNS_PREF,xml_mark[i_main],\
                LATAN_XMLNS_PREF,xml_mark[i_prop],\
                channel_id,q1_id,q2_id,source_id,sink_id);
        nodeset = xml_get_nodeset(xpath_expr,ws);
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
    }
    END_XML_PARSING(ws,fname)

    return status;
}

latan_errno hadron_prop_load_bin(mat **prop, const hadron *h,              \
                                 const ss_no source, const ss_no sink,     \
                                 const strbuf manfname,const size_t binsize)
{
    int i,p,s;
    size_t j;
    int ndat, chmix, stmix;
    size_t nt;
    double mean;
    mat **dat[MAXPROP][MAXQUARKST];
    mat **prop_prebin;
    strbuf buf,*fname;
    latan_errno status;
    
    nt      = nrow(prop[0]);
    ndat    = get_nfile(manfname);
    chmix   = ((h->chmix) == NOMIX) ? 1 : MAXPROP;
    stmix   = ((h->stmix) == NOMIX) ? 1 : MAXQUARKST;
    status  = LATAN_SUCCESS;
    
    if (ndat == -1)
    {
        LATAN_ERROR("error while reading manifest file",LATAN_ESYSTEM);
    }
    
    prop_prebin = mat_ar_create(ndat,nt,1);
    MALLOC(fname,strbuf *,ndat);
    
    for (p=0;p<chmix;p++)
    {
        for (s=0;s<stmix;s++)
        {
            dat[p][s] = mat_ar_create((size_t)(ndat),(size_t)(nt),1);
            i = 0;
            BEGIN_FOR_LINE(buf,manfname)
            {
                strcpy(fname[i],buf);
                i++;
            }
            END_FOR_LINE
#ifdef _OPENMP
            #pragma omp parallel for    
#endif
            for (i=0;i<ndat;i++)
            {
                LATAN_UPDATE_STATUS(status,prop_load(dat[p][s][i], \
                                    h->channel[p],h->quarkst[s][0],\
                                    h->quarkst[s][1],source,sink,fname[i]));
            }
#ifdef _OPENMP
            #pragma omp barrier
#endif
        }
    }
    for (i=0;i<ndat;i++)
    {
        mat_zero(prop_prebin[i]);
        for (p=0;p<chmix;p++)
        {
            for (s=0;s<stmix;s++)
            {
                LATAN_UPDATE_STATUS(status,mat_eqadd(prop_prebin[i],\
                                                     dat[p][s][i]));
            }
        }
        if ((h->chmix) == MEAN)
        {
            LATAN_UPDATE_STATUS(status,mat_eqmuls(prop_prebin[i],\
                                                  DRATIO(1,MAXPROP)));
        }
        if ((h->stmix) == MEAN)
        {
            LATAN_UPDATE_STATUS(status,mat_eqmuls(prop_prebin[i],\
                                                  DRATIO(1,MAXQUARKST)));
        }
        mat_eqabs(prop_prebin[i]);
    }
    if ((h->parity) == ODD)
    {
        for (i=0;i<ndat;i++)
        {
            for (j=1;j<nt/2;j++)
            {
                mean = 0.5*(mat_get(prop_prebin[i],(size_t)(j),0)   \
                            + mat_get(prop_prebin[i],(size_t)(nt-j),0));
                mat_set(prop_prebin[i],(size_t)(j),0,mean);
                mat_set(prop_prebin[i],(size_t)(nt-j),0,mean);
            }
        }
    }
    LATAN_UPDATE_STATUS(status,mat_ar_bin(prop,prop_prebin,ndat,binsize));
    
    mat_ar_destroy(prop_prebin,ndat);
    for (p=0;p<chmix;p++)
    {
        for (s=0;s<stmix;s++)
        {
            mat_ar_destroy(dat[p][s],(size_t)(ndat));
        }
    }
    FREE(fname);
    
    return status;
}

latan_errno hadron_prop_load_nt(size_t *nt, const hadron *h,\
                                const ss_no source,         \
                                const ss_no sink,           \
                                const strbuf manfname)
{
    latan_errno status;
    strbuf ffname;

    get_firstfname(ffname,manfname);
    status = prop_load_nt(nt,h->channel[0],h->quarkst[0][0],h->quarkst[0][1],\
                          source,sink,ffname);

    return status;
}

/*                      random generator state I/O                          */
/****************************************************************************/
latan_errno randgen_save_state(const strbuf fname, const char mode,\
                               const rg_state state, const strbuf name)
{
    xml_workspace *ws;

    XML_WRITE(ws,fname,mode,xml_insert_rgstate(ws->root,state,name));
    
    return LATAN_SUCCESS;
}

latan_errno randgen_load_state(rg_state state, const strbuf fname,\
                               const strbuf name)
{
    xml_workspace *ws;
    xmlXPathObject *nodeset;
    strbuf xpath_expr;
    latan_errno status;

    BEGIN_XML_PARSING(ws,fname)
    {
        sprintf(xpath_expr,"/%s:%s/%s:%s",LATAN_XMLNS_PREF,xml_mark[i_main],\
                LATAN_XMLNS_PREF,xml_mark[i_rgstate]);
        if (strlen(name) > 0)
        {
            sprintf(xpath_expr,"%s[@name='%s']",xpath_expr,name);
        }
        nodeset = xml_get_nodeset(xpath_expr,ws);
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
                strcpy(errmsg,"");
            }
            sprintf(errmsg,"rgstate %snot found in file %s",errmsg,fname);
            LATAN_ERROR(errmsg,LATAN_ELATSYN);
        }
    }
    END_XML_PARSING(ws,fname)

    xmlXPathFreeObject(nodeset);

    return LATAN_SUCCESS;
}

/*                          resampled sample I/O                            */
/****************************************************************************/
latan_errno rs_sample_save(const strbuf fname, const char mode,\
                           const rs_sample *s)
{
    xml_workspace *ws;
    strbuf name;

    rs_sample_get_name(name,s);
    if (strlen(name) == 0)
    {
        LATAN_ERROR("cannot save sample with an empty name",LATAN_EINVAL);
    }
    XML_WRITE(ws,fname,mode,xml_insert_sample(ws->root,s,name));
    
    return LATAN_SUCCESS;
}

latan_errno rs_sample_load_nrow(size_t *nr, const strbuf fname,\
                                const strbuf name)
{
    xml_workspace *ws;
    xmlXPathObject *nodeset;
    strbuf xpath_expr;
    latan_errno status;

    BEGIN_XML_PARSING(ws,fname)
    {
        sprintf(xpath_expr,"/%s:%s/%s:%s",LATAN_XMLNS_PREF,xml_mark[i_main],\
                LATAN_XMLNS_PREF,xml_mark[i_sample]);
        if (strlen(name) > 0)
        {
            sprintf(xpath_expr,"%s[@name='%s']",xpath_expr,name);
        }
        nodeset = xml_get_nodeset(xpath_expr,ws);
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
                strcpy(errmsg,"");
            }
            sprintf(errmsg,"sample %snot found in file %s",errmsg,fname);
            LATAN_ERROR(errmsg,LATAN_ELATSYN);
        }
    }
    END_XML_PARSING(ws,fname)

    xmlXPathFreeObject(nodeset);

    return LATAN_SUCCESS;
}

latan_errno rs_sample_load_nsample(size_t *nsample, const strbuf fname,\
                                   const strbuf name)
{
    xml_workspace *ws;
    xmlXPathObject *nodeset;
    strbuf xpath_expr;
    latan_errno status;

    BEGIN_XML_PARSING(ws,fname)
    {
        sprintf(xpath_expr,"/%s:%s/%s:%s",LATAN_XMLNS_PREF,xml_mark[i_main],\
                LATAN_XMLNS_PREF,xml_mark[i_sample]);
        if (strlen(name) > 0)
        {
            sprintf(xpath_expr,"%s[@name='%s']",xpath_expr,name);
        }
        nodeset = xml_get_nodeset(xpath_expr,ws);
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
                strcpy(errmsg,"");
            }
            sprintf(errmsg,"sample %snot found in file %s",errmsg,fname);
            LATAN_ERROR(errmsg,LATAN_ELATSYN);
        }
    }
    END_XML_PARSING(ws,fname)

    xmlXPathFreeObject(nodeset);

    return LATAN_SUCCESS;
}

latan_errno rs_sample_load(rs_sample *s, const strbuf fname, const strbuf name)
{
    xml_workspace *ws;
    xmlXPathObject *nodeset;
    strbuf xpath_expr;
    latan_errno status;

    BEGIN_XML_PARSING(ws,fname)
    {
        sprintf(xpath_expr,"/%s:%s/%s:%s",LATAN_XMLNS_PREF,xml_mark[i_main],\
                LATAN_XMLNS_PREF,xml_mark[i_sample]);
        if (strlen(name) > 0)
        {
            sprintf(xpath_expr,"%s[@name='%s']",xpath_expr,name);
        }
        nodeset = xml_get_nodeset(xpath_expr,ws);
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
                strcpy(errmsg,"");
            }
            sprintf(errmsg,"sample %snot found in file %s",errmsg,fname);
            LATAN_ERROR(errmsg,LATAN_ELATSYN);
        }
    }
    END_XML_PARSING(ws,fname)

    xmlXPathFreeObject(nodeset);

    return LATAN_SUCCESS;
}