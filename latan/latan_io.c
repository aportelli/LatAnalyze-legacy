#include <latan/latan_io.h>
#include <latan/latan_includes.h>
#include <latan/latan_math.h>
#include <latan/latan_xml.h>

/*                          general I/O                                     */
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
    strbuf channel_id,q1_id,q2_id,source_id,sink_id,errmsg,xpath_expr;
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
    strbuf channel_id,q1_id,q2_id,source_id,sink_id,errmsg,xpath_expr;
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
latan_errno randgen_save_state(const strbuf f_name,\
                               const rg_state state)
{
    strbuf full_f_name;
    FILE *f;
    size_t i;
    
    sprintf(full_f_name,"%s.rand",f_name);
    FOPEN(f,full_f_name,"w");
    for (i=0;i<RLXG_STATE_SIZE;i++)
    {
        fprintf(f,"%d ",state[i]);
    }
    fprintf(f,"\n");
    fclose(f);
    
    return LATAN_SUCCESS;
}

latan_errno randgen_load_state(rg_state state, const strbuf f_name)
{
    strbuf errmsg;
    FILE *f;
    size_t i;
    
    FOPEN(f,f_name,"r");
    for (i=0;i<RLXG_STATE_SIZE;i++)
    {
        if(fscanf(f,"%d ",state+i)<0)
        {
            sprintf(errmsg,"error while reading generator state component %lu in file %s",\
                    (unsigned long)i,f_name);
            LATAN_ERROR(errmsg,LATAN_ELATSYN);
        }
    }
    fclose(f);
    
    return LATAN_SUCCESS;
}

/*                          resampled sample I/O                            */
/****************************************************************************/
latan_errno rs_sample_save(const rs_sample *s, const strbuf f_name)
{
    strbuf s_name;
    FILE* f;
    size_t i,j;
    size_t s_nrow, s_nsample;
    
    s_nrow    = rs_sample_get_nrow(s);
    s_nsample = rs_sample_get_nsample(s);
    rs_sample_get_name(s_name,s);
    
    if (strlen(s->name) == 0)
    {
        LATAN_ERROR("cannot save sample with an empty name",LATAN_EINVAL);
    }

/*
    FOPEN(f,good_f_name,"w");
    fprintf(f,"%s %lu %lu\n",s_name,(unsigned long)s_nrow,\
            (unsigned long)s_nsample);
    for (i=0;i<s_nrow;i++)
    {
        fprintf(f,"%.10e ",mat_get(rs_sample_pt_cent_val(s),i,0));
        for (j=0;j<s_nsample-1;j++)
        {
            fprintf(f,"%.10e ",mat_get(rs_sample_pt_sample(s,j),i,0));
        }
        fprintf(f,"%.10e\n",mat_get(rs_sample_pt_sample(s,s_nsample-1),i,0));
    }
    fclose(f);
*/
    
    return LATAN_SUCCESS;
}

int rs_sample_load_nrow(const strbuf f_name)
{
    FILE* f;
    strbuf dumstr,errmsg;
    int s_nrow;
    
    FOPEN(f,f_name,"r");
    sprintf(errmsg,"error while reading sample dimension in file %s",f_name);
    if (fscanf(f,"%s ",dumstr)<0)
    {
        LATAN_ERROR(errmsg,LATAN_FAILURE);
    }
    if (fscanf(f,"%d ",&s_nrow)<0)
    {
        LATAN_ERROR(errmsg,LATAN_FAILURE);
    }
    fclose(f);
    
    return s_nrow;
}

int rs_sample_load_nsample(const strbuf f_name)
{
    FILE* f;
    strbuf dumstr,errmsg;
    int nsample;
    
    FOPEN(f,f_name,"r");
    sprintf(errmsg,"error while reading number of sample in file %s",f_name);
    if (fscanf(f,"%s ",dumstr)<0)
    {
        LATAN_ERROR(errmsg,LATAN_FAILURE);
    }
    if (fscanf(f,"%d ",&nsample)<0)
    {
        LATAN_ERROR(errmsg,LATAN_FAILURE);
    }
    if (fscanf(f,"%d ",&nsample)<0)
    {
        LATAN_ERROR(errmsg,LATAN_FAILURE);
    }
    fclose(f);
    
    return nsample;
}

latan_errno rs_sample_load(rs_sample *s, const strbuf f_name)
{
    FILE* f;
    size_t i,j;
    strbuf dumstr,errmsg;
    double dbuf;
    int ibuf;
    mat *pt;
    
    FOPEN(f,f_name,"r");
    if (fscanf(f,"%s ",dumstr)<0)
    {
        sprintf(errmsg,"error while reading sample name in file %s",f_name);
        LATAN_ERROR(errmsg,LATAN_ELATSYN);
    }
    strcpy(s->name,dumstr);
    if (fscanf(f,"%d ",&ibuf)<0)
    {
        sprintf(errmsg,"error while reading sample dimension in file %s",\
                f_name);
        LATAN_ERROR(errmsg,LATAN_ELATSYN);
    }
    if ((size_t)ibuf != rs_sample_get_nrow(s))
    {
        sprintf(errmsg,"sample dimension (%d) in file %s is invalid",ibuf,\
                f_name);
        LATAN_ERROR(errmsg,LATAN_EBADLEN);
    }
    if (fscanf(f,"%d ",&ibuf)<0)
    {
        sprintf(errmsg,"error while reading number of sample in file %s",\
                f_name);
        LATAN_ERROR(errmsg,LATAN_ELATSYN);
    }
    if ((size_t)ibuf != rs_sample_get_nsample(s))
    {
        sprintf(errmsg,"number of sample (%d) in file %s is invalid",ibuf,\
                f_name);
        LATAN_ERROR(errmsg,LATAN_EBADLEN);
    }
    for (i=0;i<rs_sample_get_nrow(s);i++)
    {
        for (j=0;j<rs_sample_get_nsample(s);j++)
        {
            if (fscanf(f,"%lf ",&dbuf)<0)
            {
                sprintf(errmsg,"error reading component %lu of sample %lu in file %s",\
                        (unsigned long)i,(unsigned long)j,f_name);
                LATAN_ERROR(errmsg,LATAN_ELATSYN);
            }
            if (j == 0)
            {
                pt = rs_sample_pt_cent_val(s);
            }
            else
            {
                pt = rs_sample_pt_sample(s,j-1);
            }
            mat_set(pt,i,0,dbuf);
        }
        j = rs_sample_get_nsample(s);
        if (fscanf(f,"%lf\n",&dbuf)<0)
        {
            sprintf(errmsg,"error reading component %lu of sample %lu in file %s",\
                    (unsigned long)i,(unsigned long)j,f_name);
            LATAN_ERROR(errmsg,LATAN_ELATSYN);
        }
        mat_set(rs_sample_pt_sample(s,j-1),i,0,dbuf);
    }
    fclose(f);
    
    return LATAN_SUCCESS;
}