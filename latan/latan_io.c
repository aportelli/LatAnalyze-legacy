/* latan_io.c, part of LatAnalyze library
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

#include <latan/latan_io.h>
#include <latan/latan_includes.h>
#include <latan/latan_math.h>

/*                            default I/O functions                         */
/****************************************************************************/
#define DEF_IO_INIT                io_init_xml
#define DEF_IO_FINISH              io_finish_xml
#define DEF_PROP_LOAD_NT           prop_load_nt_xml
#define DEF_PROP_LOAD              prop_load_xml
#define DEF_RANDGEN_SAVE_STATE     randgen_save_state_xml
#define DEF_RANDGEN_LOAD_STATE     randgen_load_state_xml
#define DEF_RS_SAMPLE_SAVE         rs_sample_save_xml
#define DEF_RS_SAMPLE_LOAD_NROW    rs_sample_load_nrow_xml
#define DEF_RS_SAMPLE_LOAD_NSAMPLE rs_sample_load_nsample_xml
#define DEF_RS_SAMPLE_LOAD         rs_sample_load_xml

/*                              I/O init/finish                             */
/****************************************************************************/
void (*io_init)(void)   = &DEF_IO_INIT;
void (*io_finish)(void) = &DEF_IO_FINISH;

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
latan_errno (*prop_load_nt)(size_t *nt, const channel_no channel,\
                            const quark_no q1, const quark_no q2,\
                            const ss_no source, const ss_no sink,\
                            strbuf fname)                        \
        = &DEF_PROP_LOAD_NT;
latan_errno (*prop_load)(mat *prop, const channel_no channel, \
                         const quark_no q1, const quark_no q2,\
                         const ss_no source, const ss_no sink,\
                         strbuf fname)                        \
        = &DEF_PROP_LOAD;

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
                STRBUFCPY(fname[i],buf);
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
latan_errno (*randgen_save_state)(const strbuf f_name, const char mode,   \
                                  const rg_state state, const strbuf name)\
        = &DEF_RANDGEN_SAVE_STATE;
latan_errno (*randgen_load_state)(rg_state state, const strbuf f_name,\
                                  const strbuf name)                  \
        = &DEF_RANDGEN_LOAD_STATE;

/*                          resampled sample I/O                            */
/****************************************************************************/
latan_errno (*rs_sample_save)(const strbuf fname, const char mode,\
                              const rs_sample *s)                 \
        = &DEF_RS_SAMPLE_SAVE;
latan_errno (*rs_sample_load_nrow)(size_t *nr, const strbuf fname,\
                                   const strbuf name)             \
                                                                  \
        = &DEF_RS_SAMPLE_LOAD_NROW;
latan_errno (*rs_sample_load_nsample)(size_t *nsample, const strbuf fname,\
                                      const strbuf name)                  \
        = &DEF_RS_SAMPLE_LOAD_NSAMPLE;
latan_errno (*rs_sample_load)(rs_sample *s, const strbuf fname,\
                              const strbuf name)               \
        = &DEF_RS_SAMPLE_LOAD;
