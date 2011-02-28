/* latan_statistics.c, part of LatAnalyze library
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

#include <latan/latan_statistics.h>
#include <latan/latan_includes.h>
#include <latan/latan_math.h>
#include <latan/latan_rand.h>
#include <latan/latan_io.h>
#include <latan/latan_mass.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_errno.h>

static latan_errno resample_bootstrap(mat *cent_val, mat **sample,          \
                                      const size_t nboot, mat **dat,    \
                                      const size_t ndat, const size_t nobs, \
                                      rs_func *f, void *param);

/* TODO : Jackknife resampling function
static latan_errno resample_jackknife(mat *cent_val, mat **sample,          \
                                      const size_t jk_depth, mat **dat,\
                                      const size_t ndat, const size_t nobs, \
                                      rs_func *f, void *param);
*/

/*                      elementary estimators                               */
/****************************************************************************/
double mat_elsum(mat *m)
{
    size_t i,j;
    double sum;
    
    sum = 0.0;
    
    FOR_VAL(m,i,j)
    {
        sum += mat_get(m,i,j);
    }
    
    return sum;
}

double mat_elmean(mat *m)
{
    double mean;
    
    mean = mat_elsum(m)/((double)(nrow(m)*ncol(m)));
    
    return mean;
}

latan_errno mat_mean(mat *mean, mat **m, const size_t size)
{
    latan_errno status;
    size_t i;
    const double dsize = (double)(size);
    
    status = LATAN_SUCCESS;
    mat_zero(mean);
    
    for (i=0;i<size;i++)
    {
        LATAN_UPDATE_STATUS(status,mat_eqadd(mean,m[i]));
    }
    LATAN_UPDATE_STATUS(status,mat_eqmuls(mean,1.0/dsize));
    
    return status;
}

latan_errno mat_cov(mat *cov, mat **m, mat **n, const size_t size)
{
    latan_errno status;
    mat *m_mean,*n_mean;
    
    status = LATAN_SUCCESS;
    
    m_mean = mat_create_from_dim(m[0]);
    n_mean = mat_create_from_dim(n[0]);
    
    LATAN_UPDATE_STATUS(status,mat_mean(m_mean,m,size));
    LATAN_UPDATE_STATUS(status,mat_mean(n_mean,n,size));
    LATAN_UPDATE_STATUS(status,mat_cov_m(cov,m,n,size,m_mean,n_mean));
    
    mat_destroy(m_mean);
    mat_destroy(n_mean);
    
    return status;
}

latan_errno mat_cov_m(mat *cov, mat **m, mat **n, const size_t size,\
                      mat *m_mean, mat *n_mean)
{
    latan_errno status;
    size_t i;
    size_t subdim;
    double dsubdim;
    mat **mctnc;
    mat **mc;
    mat **nc;
    
    status = LATAN_SUCCESS;
    subdim = ncol(m[0]);
    dsubdim = (double)(subdim);
    
    mctnc = mat_ar_create_from_dim(size,cov);
    mc    = mat_ar_create_from_dim(size,m[0]);
    nc    = mat_ar_create_from_dim(size,n[0]);
    
    for (i=0;i<size;i++) 
    {
        LATAN_UPDATE_STATUS(status,mat_sub(mc[i],m[i],m_mean));
        LATAN_UPDATE_STATUS(status,mat_sub(nc[i],n[i],n_mean));
        LATAN_UPDATE_STATUS(status,mat_mul(mctnc[i],mc[i],'n',nc[i],'t'));
        LATAN_UPDATE_STATUS(status,mat_eqmuls(mctnc[i],1.0/dsubdim));
    }
    LATAN_UPDATE_STATUS(status,mat_mean(cov,mctnc,size));
    
    mat_ar_destroy(mctnc,size);
    mat_ar_destroy(mc,size);
    mat_ar_destroy(nc,size);
    
    return status;
}

latan_errno mat_covp(mat *cov, mat **m, mat **n, const size_t size)
{
    latan_errno status;
    mat *m_mean,*n_mean;
    
    status = LATAN_SUCCESS;
    
    m_mean = mat_create_from_dim(m[0]);
    n_mean = mat_create_from_dim(n[0]);
    
    LATAN_UPDATE_STATUS(status,mat_mean(m_mean,m,size));
    LATAN_UPDATE_STATUS(status,mat_mean(n_mean,n,size));
    LATAN_UPDATE_STATUS(status,mat_covp_m(cov,m,n,size,m_mean,n_mean));
    
    mat_destroy(m_mean);
    mat_destroy(n_mean);
    
    return status;
}

latan_errno mat_covp_m(mat *cov, mat **m, mat **n, const size_t size,\
                       mat *m_mean, mat *n_mean)
{
    latan_errno status;
    size_t i;
    mat **mcnc;
    mat **nc;
    
    status = LATAN_SUCCESS;
    
    mcnc = mat_ar_create_from_dim(size,cov);
    nc = mat_ar_create_from_dim(size,cov);
    
    for (i=0;i<size;i++)
    {
        LATAN_UPDATE_STATUS(status,mat_sub(mcnc[i],m[i],m_mean));
        LATAN_UPDATE_STATUS(status,mat_sub(nc[i],n[i],n_mean));
        LATAN_UPDATE_STATUS(status,mat_eqmulp(mcnc[i],nc[i]));
    }
    LATAN_UPDATE_STATUS(status,mat_mean(cov,mcnc,size));
    
    mat_ar_destroy(mcnc,size);
    mat_ar_destroy(nc,size);

    return status;
}

/*                              data binning                                */
/****************************************************************************/
latan_errno mat_ar_bin(mat **bindat, mat **dat, const size_t ndat,\
                       const size_t binsize)
{
    size_t nbindat;
    size_t i;
    latan_errno status;
    
    if (ndat % binsize != 0)
    {
        LATAN_ERROR("bin size is not a divider of data number",LATAN_EBADLEN);
    }
    
    nbindat = ndat/binsize;
    status = LATAN_SUCCESS;
    
    for (i=0;i<nbindat;i++)
    {
        LATAN_UPDATE_STATUS(status,mat_mean(bindat[i],dat+i*binsize,binsize));
    }
    
    return status;
}

/*                              histogram                                   */
/****************************************************************************/
latan_errno histogram(mat *hist, mat *data, const double xmin,\
                      const double xmax, const size_t nint)
{
    size_t i;
    gsl_histogram *gsl_hist;
    gsl_error_handler_t *error_handler;
    
    
    if (nrow(hist) != nint)
    {
        LATAN_ERROR("histogram matrix row number do not match number of histogram intervals",\
                    LATAN_EBADLEN);
    }
    
    error_handler = gsl_set_error_handler(&latan_error);
    gsl_hist = gsl_histogram_alloc(nint);
    gsl_set_error_handler(error_handler);
    mat_zero(hist);
    
    gsl_histogram_set_ranges_uniform(gsl_hist,xmin,xmax);
    for (i=0;i<nrow(data);i++)
    {
        gsl_histogram_increment(gsl_hist,mat_get(data,i,0));
    }
    mat_set_from_ar(hist,gsl_hist->bin);
    
    gsl_histogram_free(gsl_hist);
    
    return LATAN_SUCCESS;
}

/*                  resampled samples manipulation                          */
/****************************************************************************/
/** jackknife sample number calculation **/
size_t jackknife_nsample(const size_t ndat, const size_t jk_depth)
{
    unsigned int i;
    size_t nsample;
    const unsigned int ui_ndat = (unsigned int)(ndat);
    
    nsample = 0;

    for (i=1;i<=jk_depth;i++)
    {
        nsample += (size_t)(binomial(ui_ndat,i));
    }
    
    return nsample;
}

/** allocation **/
rs_sample *rs_sample_create(const size_t init_nrow, const size_t nsample)
{
    rs_sample *s;

    MALLOC_ERRVAL(s,rs_sample *,1,NULL);

    s->nsample = nsample;

    s->cent_val = mat_create(init_nrow,1);
    s->sample   = mat_ar_create_from_dim(s->nsample,s->cent_val);

    return s;
}

void rs_sample_destroy(rs_sample *s)
{
    mat_destroy(s->cent_val);
    mat_ar_destroy(s->sample,s->nsample);
    FREE(s);
}

/** access **/
size_t rs_sample_get_nrow(const rs_sample *s)
{
    return nrow(s->cent_val);
}

size_t rs_sample_get_nsample(const rs_sample *s)
{
    return s->nsample;
}

void rs_sample_get_name(strbuf name, const rs_sample *s)
{
    strbufcpy(name,s->name);
}

mat *rs_sample_pt_cent_val(const rs_sample *s)
{
    return s->cent_val;
}

mat *rs_sample_pt_sample(const rs_sample *s, const size_t i)
{
    return (s->sample)[i];
}

void rs_sample_set_name(rs_sample *s, const strbuf name)
{
    strbufcpy(s->name,name);
}

/** estimators **/
latan_errno rs_sample_cov(mat *cov, const rs_sample *s, const rs_sample *t)
{
    latan_errno status;
    
    if (s->nsample != t->nsample)
    {
        LATAN_ERROR("operation between resampled samples with dimension mismatch",\
                    LATAN_EBADLEN);
    }
    
    status = mat_cov_m(cov,s->sample,t->sample,s->nsample,\
                       s->cent_val,t->cent_val);
    
    return status;
}

latan_errno rs_sample_covp(mat *cov, const rs_sample *s, const rs_sample *t)
{
    latan_errno status;
    
    if (s->nsample != t->nsample)
    {
        LATAN_ERROR("operation between resampled samples with dimension mismatch",\
                    LATAN_EBADLEN);
    }
    
    status = mat_covp_m(cov,s->sample,t->sample,s->nsample,\
                        s->cent_val,t->cent_val);
    
    return status;
}
/*                      resampling functions                                */
/****************************************************************************/
static latan_errno resample_bootstrap(mat *cent_val, mat **sample,         \
                                      const size_t nboot, mat **dat,       \
                                      const size_t ndat, const size_t nobs,\
                                      rs_func *f, void *param)
{
    mat **fakedat;
    size_t i,j,k;
    unsigned int rj;
    latan_errno status;
    
    status = LATAN_SUCCESS;
    
    MALLOC(fakedat,mat**,ndat*nobs);
    
    LATAN_UPDATE_STATUS(status,f(cent_val,dat,ndat,0,param));
    for (i=0;i<nboot;i++)
    {
        
        for (j=0;j<ndat;j++) 
        {
            rj = rand_ud((unsigned int)(ndat));
            for (k=0; k<nobs; k++)
            {
                fakedat[j+k*ndat] = dat[rj+k*ndat];
            }
        }
        LATAN_UPDATE_STATUS(status,f(sample[i],fakedat,ndat,i+1,param));
    }
    
    FREE(fakedat);
    
    return status;
}

latan_errno resample(rs_sample *s, mat **dat, const size_t ndat,               \
                     const size_t nobs, rs_func *f, unsigned int resamp_method,\
                     void *param)
{
    latan_errno status;

    /** bootstrap **/
    if (resamp_method == 0)
    {
        randgen_get_state(s->gen_state);
        status = resample_bootstrap(s->cent_val,s->sample,s->nsample,dat,ndat,\
                                    nobs,f,param);
    }
    /** jacknife **/
    else
    {
        if (resamp_method >= ndat)
        {
            LATAN_ERROR("jackknife resampling depth too large",LATAN_EINVAL);
        }
        LATAN_ERROR("jackknife resampling is not implemented yet",\
                    LATAN_FAILURE);
    }

    return status;
}

/* some basic rs_func */
latan_errno rs_mean(mat *res, mat **dat, const size_t ndat,\
                    const size_t sampno, void *nothing)
{
    latan_errno status;
    size_t st_nothing;
    
    nothing    = NULL;
    st_nothing = sampno;
    
    status = mat_mean(res,dat,ndat);
    
    return status;
}

latan_errno rs_finite_diff(mat *res, mat **dat, const size_t ndat,\
                           const size_t sampno, void *nothing)
{
    latan_errno status;
    mat *mean;
    size_t st_nothing;
    
    nothing    = NULL;
    st_nothing = sampno;
    status     = LATAN_SUCCESS;
    
    mean = mat_create_from_dim(dat[0]);
    
    LATAN_UPDATE_STATUS(status,mat_mean(mean,dat,ndat));
    LATAN_UPDATE_STATUS(status,finite_diff(res,mean));
    
    mat_destroy(mean);
    
    return status;
}

latan_errno rs_effmass(mat *res, mat **dat, const size_t ndat,\
                       const size_t sampno, void *parity)
{
    latan_errno status;
    mat *mean;
    int parityt;
    size_t st_nothing;
    
    st_nothing = sampno;
    status     = LATAN_SUCCESS;
    parityt    = *((int*)(parity));

    mean = mat_create_from_dim(dat[0]);
    
    LATAN_UPDATE_STATUS(status,mat_mean(mean,dat,ndat));
    LATAN_UPDATE_STATUS(status,effmass(res,mean,parityt));
    
    mat_destroy(mean);
    
    return status;
}

latan_errno rs_effmass_PCAC(mat *res, mat **dat, const size_t ndat,\
                            const size_t sampno, void *nothing)
{
    latan_errno status;
    mat **prop_AP;
    mat **prop_PP;
    mat *mprop_AP,*mprop_PP;
    size_t st_nothing;
    
    nothing    = NULL;
    st_nothing = sampno;
    status     = LATAN_SUCCESS;
    
    prop_AP  = dat;
    prop_PP  = dat + ndat;
    mprop_AP = mat_create_from_dim(prop_AP[0]);
    mprop_PP = mat_create_from_dim(prop_PP[0]);
    
    LATAN_UPDATE_STATUS(status,mat_mean(mprop_AP,prop_AP,ndat));
    LATAN_UPDATE_STATUS(status,mat_mean(mprop_PP,prop_PP,ndat));
    LATAN_UPDATE_STATUS(status,effmass_PCAC(res,mprop_AP,mprop_PP));
    
    mat_destroy(mprop_AP);
    mat_destroy(mprop_PP);
    
    return status;
}