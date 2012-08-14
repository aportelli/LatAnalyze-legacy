/* latan_statistics.c, part of LatAnalyze library
 *
 * Copyright (C) 2010, 2011, 2012 Antonin Portelli
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
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sort_double.h>

static latan_errno resample_bootstrap(mat *cent_val, mat **sample,         \
                                      const size_t nboot, mat **dat,       \
                                      const size_t ndat, rs_func *f,       \
                                      void *param);

/* TODO : Jackknife resampling function
static latan_errno resample_jackknife(mat *cent_val, mat **sample,         \
                                      const size_t jk_depth, mat **dat,    \
                                      const size_t ndat, const size_t nobs,\
                                      rs_func *f, void *param);
*/

/*                      elementary estimators                               */
/****************************************************************************/
double mat_elsum(const mat *m)
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

double mat_elmean(const mat *m)
{
    double mean;
    
    mean = mat_elsum(m)/((double)(nel(m)));
    
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
        USTAT(mat_eqadd(mean,m[i]));
    }
    USTAT(mat_eqmuls(mean,1.0/dsize));
    
    return status;
}

latan_errno mat_cov(mat *cov, mat **m, mat **n, const size_t size)
{
    latan_errno status;
    mat *m_mean,*n_mean;
    
    status = LATAN_SUCCESS;
    
    m_mean = mat_create_from_dim(m[0]);
    n_mean = mat_create_from_dim(n[0]);
    
    USTAT(mat_mean(m_mean,m,size));
    USTAT(mat_mean(n_mean,n,size));
    USTAT(mat_cov_m(cov,m,n,size,m_mean,n_mean));
    
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
        USTAT(mat_sub(mc[i],m[i],m_mean));
        USTAT(mat_sub(nc[i],n[i],n_mean));
        USTAT(mat_mul(mctnc[i],mc[i],'n',nc[i],'t'));
        USTAT(mat_eqmuls(mctnc[i],1.0/dsubdim));
    }
    USTAT(mat_mean(cov,mctnc,size));
    
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
    
    USTAT(mat_mean(m_mean,m,size));
    USTAT(mat_mean(n_mean,n,size));
    USTAT(mat_covp_m(cov,m,n,size,m_mean,n_mean));
    
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
        USTAT(mat_sub(mcnc[i],m[i],m_mean));
        USTAT(mat_sub(nc[i],n[i],n_mean));
        USTAT(mat_eqmulp(mcnc[i],nc[i]));
    }
    USTAT(mat_mean(cov,mcnc,size));
    
    mat_ar_destroy(mcnc,size);
    mat_ar_destroy(nc,size);

    return status;
}

/*                               percentiles                                */
/****************************************************************************/
double ar_percentile(const double *data, const double *w, const size_t *sind,\
                     const size_t ndata, const double p)
{
    double w_totsum, w_psum, w_i, p_i, p_im1, res;
    bool have_res;
    size_t i;
    
    if ((p < 0.0)||(p > 100.0))
    {
        LATAN_WARNING("percentile is outside the [0,100] range",LATAN_EINVAL);
    }
    
    /* compute percentile                             */
    /* (cf. http://en.wikipedia.org/wiki/Percentile ) */
    if (w)
    {
        w_totsum = 0.0;
        for (i=0;i<ndata;i++)
        {
            w_totsum += w[i];
        }
        w_psum = w[sind[0]];
    }
    else
    {
        w_totsum = (double)(ndata);
        w_psum   = 1.0;
    }
    p_i = (100.0/w_totsum)*w_psum*0.5;
    if (p < p_i)
    {
        res = data[sind[0]];
    }
    else
    {
        have_res = false;
        p_im1    = p_i;
        for (i=1;i<ndata;i++)
        {
            if (w)
            {
                w_i = w[sind[i]];
            }
            else
            {
                w_i = 1.0;
            }
            w_psum += w_i;
            p_i     = (100.0/w_totsum)*(w_psum-0.5*w_i);
            if ((p >= p_im1)&&(p < p_i))
            {
                res      = data[sind[i-1]]+(p-p_im1)/(p_i-p_im1)\
                           *(data[sind[i]]-data[sind[i-1]]);
                have_res = true;
                break;
            }
        }
        if (!have_res)
        {
            res = data[sind[ndata-1]];
        }
    }
    
    return res;
}

double mat_elpercentile_with_sind(const mat *m, const mat *w,       \
                                  const size_t *sind, const double p)
{
    mat *m_buf,*w_buf;
    const mat *m_pt, *w_pt;
    bool create_m_buf,create_w_buf;
    double *m_ar,*w_ar;
    double res;
    
    m_buf = NULL;
    w_buf = NULL;
    
    /* create buffers if matrices are submatrices */
    create_m_buf = (ncol(m) != m->data_cpu->tda);
    if (w)
    {
        create_w_buf = (ncol(w) != w->data_cpu->tda);
    }
    else
    {
        create_w_buf = false;
    }
    if (create_m_buf)
    {
        m_buf = mat_create_from_mat(m);
        m_pt  = m_buf;
        LATAN_WARNING("buffer created for data matrix",LATAN_EINVAL);
    }
    else
    {
        m_pt  = m;
    }
    if (create_w_buf)
    {
        w_buf = mat_create_from_mat(w);
        w_pt  = w_buf;
        LATAN_WARNING("buffer created for weight matrix",LATAN_EINVAL);
    }
    else
    {
        w_pt  = w;
    }
    
    /* compute percentile */
    m_ar = m_pt->data_cpu->data;
    if (w)
    {
        w_ar = w_pt->data_cpu->data;
    }
    else
    {
        w_ar = NULL;
    }
    res  = ar_percentile(m_ar,w_ar,sind,nel(m_pt),p);
    
    /* deallocation */
    if (create_m_buf)
    {
        mat_destroy(m_buf);
    }
    if (create_w_buf)
    {
        mat_destroy(w_buf);
    }
    
    return res;
}

double mat_elpercentile(const mat *m, const mat *w, const double p)
{
    size_t *sind;
    double res;
    
    MALLOC(sind,size_t *,nel(m));
    
    mat_get_sind(sind,m);
    res = mat_elpercentile_with_sind(m,w,sind,p);
    
    FREE(sind);
    
    return res;
}

/*                               chi^2 p-value                              */
/****************************************************************************/
double chi2_pvalue(const double chi2_val, const size_t ndof)
{
    return gsl_cdf_chisq_Q(chi2_val,(double)ndof);
}

/*                            confidence interval                           */
/****************************************************************************/
double conf_int(double ci[2], const mat *m, const mat *w, const double nsig)
{
    double cl,pl,pr;
    size_t *sind;
    
    MALLOC(sind,size_t *,nel(m));
    
    mat_get_sind(sind,m);
    cl    = gsl_sf_erf(nsig/sqrt(2.0));
    pl    = 50.0*(1.0-cl);
    pr    = 50.0*(1.0+cl);
    ci[0] = mat_elpercentile_with_sind(m,w,sind,pl);
    ci[1] = mat_elpercentile_with_sind(m,w,sind,pr);
    
    FREE(sind);
    
    return cl;
}

/*                              data binning                                */
/****************************************************************************/
latan_errno mat_ar_bin(mat **bindat, mat **dat, const size_t ndat,\
                       const size_t binsize)
{
    size_t i;
    size_t q,r;
    latan_errno status;
    
    status = LATAN_SUCCESS;
    q      = ndat/binsize;
    r      = ndat%binsize;
    
    for (i=0;i<q;i++)
    {
        USTAT(mat_mean(bindat[i],dat+i*binsize,binsize));
    }
    if (r != 0)
    {
        USTAT(mat_mean(bindat[q],dat+q*binsize,r));
    }
    
    return status;
}

/*                              histogram                                   */
/****************************************************************************/
latan_errno histogram(mat *hist, const mat *data, const mat* w,\
                      const double xmin, const double xmax, const size_t nbin)
{
    size_t i;
    gsl_histogram *gsl_hist;
    double w_i;
    
    if (nrow(hist) != nbin)
    {
        LATAN_ERROR("histogram matrix row number do not match number of histogram bins",\
                    LATAN_EBADLEN);
    }
    
    gsl_hist = gsl_histogram_alloc(nbin);

    mat_zero(hist);
    gsl_histogram_set_ranges_uniform(gsl_hist,xmin,xmax);
    for (i=0;i<nrow(data);i++)
    {
        if (w)
        {
            w_i = mat_get(w,i,0);
        }
        else
        {
            w_i = 1.0;
        }
        gsl_histogram_accumulate(gsl_hist,mat_get(data,i,0),w_i);
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
rs_sample *rs_sample_create(const size_t init_nrow, const size_t init_ncol,\
                            const size_t nsample)
{
    rs_sample *s;

    MALLOC_ERRVAL(s,rs_sample *,1,NULL);

    s->nsample = nsample;

    s->cent_val = mat_create(init_nrow,init_ncol);
    s->sample   = mat_ar_create_from_dim(s->nsample,s->cent_val);

    return s;
}

void rs_sample_destroy(rs_sample *s)
{
    if (s != NULL)
    {
        mat_destroy(s->cent_val);
        mat_ar_destroy(s->sample,s->nsample);
        FREE(s);
    }
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

mat *rs_sample_pt_cent_val(const rs_sample *s)
{
    return s->cent_val;
}

mat *rs_sample_pt_sample(const rs_sample *s, const size_t i)
{
    if (i >= s->nsample)
    {
        LATAN_ERROR_VAL("sample index out of range",LATAN_EBADLEN,NULL);
    }
    
    return (s->sample)[i];
}

latan_errno rs_sample_get_subsamp(rs_sample *s_a, const rs_sample *s_b,\
                                  const size_t k1, const size_t l1,    \
                                  const size_t k2, const size_t l2)
{
    latan_errno status;
    size_t nsample;
    size_t i;
    
    if (s_a->nsample != s_b->nsample)
    {
        LATAN_ERROR("operation between samples with different numbers of elements",\
                    LATAN_EINVAL);
    }
    
    status  = LATAN_SUCCESS;
    nsample = s_a->nsample;
    
    USTAT(mat_get_subm(s_a->cent_val,s_b->cent_val,k1,l1,k2,l2));
    for (i=0;i<nsample;i++)
    {
        USTAT(mat_get_subm(s_a->sample[i],s_b->sample[i],k1,l1,k2,l2));
    }
    
    return status;
}

latan_errno rs_sample_set_subsamp(rs_sample *s_a, const rs_sample *s_b,\
                                  const size_t k1, const size_t l1,    \
                                  const size_t k2, const size_t l2)
{
    latan_errno status;
    size_t nsample;
    size_t i;
    
    if (s_a->nsample != s_b->nsample)
    {
        LATAN_ERROR("operation between samples with different numbers of elements",\
                    LATAN_EINVAL);
    }
    
    status  = LATAN_SUCCESS;
    nsample = s_a->nsample;
    
    USTAT(mat_set_subm(s_a->cent_val,s_b->cent_val,k1,l1,k2,l2));
    for (i=0;i<nsample;i++)
    {
        USTAT(mat_set_subm(s_a->sample[i],s_b->sample[i],k1,l1,k2,l2));
    }
    
    return status;
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
    
    status = mat_cov(cov,s->sample,t->sample,s->nsample);
    
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
    
    status = mat_covp(cov,s->sample,t->sample,s->nsample);
    
    return status;
}
/*                      resampling functions                                */
/****************************************************************************/
static latan_errno resample_bootstrap(mat *cent_val, mat **sample,         \
                                      const size_t nboot, mat **dat,       \
                                      const size_t ndat, rs_func *f,       \
                                      void *param)
{
    mat **fakedat;
    size_t i,j;
    unsigned int rj;
    latan_errno status;
    
    status = LATAN_SUCCESS;
    
    MALLOC(fakedat,mat**,ndat);
    
    USTAT(f(cent_val,dat,ndat,param));
    for (i=0;i<nboot;i++)
    {
        
        for (j=0;j<ndat;j++) 
        {
            rj = rand_ud((unsigned int)(ndat));
            fakedat[j] = dat[rj];
        }
        USTAT(f(sample[i],fakedat,ndat,param));
    }
    
    FREE(fakedat);
    
    return status;
}

latan_errno resample(rs_sample *s, mat **dat, const size_t ndat, rs_func *f, \
                     unsigned int resamp_method, void *param)
{
    latan_errno status;

    /** bootstrap **/
    if (resamp_method == 0)
    {
        status = resample_bootstrap(s->cent_val,s->sample,s->nsample,dat,ndat,\
                                    f,param);
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

/*                              useful rs_func                              */
/****************************************************************************/
latan_errno rs_mean(mat *res, mat **dat, const size_t ndat,\
                    void *nothing __dumb)
{
    latan_errno status;

    status  = mat_mean(res,dat,ndat);
    
    return status;
}

/*                              operations                                  */
/****************************************************************************/
#define CENT_VAL(s)  rs_sample_pt_cent_val(s)
#define ELEMENT(s,i) rs_sample_pt_sample(s,i)

latan_errno rs_sample_unop(rs_sample *s_a, const rs_sample *s_b, mat_unop *f)
{
    size_t i;
    size_t nsample;
    latan_errno status;
    
    if (rs_sample_get_nsample(s_a) != rs_sample_get_nsample(s_b))
    {
        LATAN_ERROR("operation between samples with different numbers of elements",\
                    LATAN_EINVAL);
    }
    
    nsample = rs_sample_get_nsample(s_a);
    status  = LATAN_SUCCESS;
    
    USTAT(f(CENT_VAL(s_a),CENT_VAL(s_b)));
    for (i=0;i<nsample;i++)
    {
        USTAT(f(ELEMENT(s_a,i),ELEMENT(s_b,i)));
    }
    
    return status;
}

latan_errno rs_sample_unops(rs_sample *s_a, const double s, mat_unops *f)
{
    size_t i;
    size_t nsample;
    latan_errno status;
    
    nsample = rs_sample_get_nsample(s_a);
    status  = LATAN_SUCCESS;
    
    USTAT(f(CENT_VAL(s_a),s));
    for (i=0;i<nsample;i++)
    {
        USTAT(f(ELEMENT(s_a,i),s));
    }
    
    return status;
}

latan_errno rs_sample_binop(rs_sample *s_a, const rs_sample *s_b, \
                            const rs_sample *s_c, mat_binop *f)
{
    size_t i;
    size_t nsample;
    latan_errno status;
    
    if ((rs_sample_get_nsample(s_a) != rs_sample_get_nsample(s_b))\
      ||(rs_sample_get_nsample(s_b) != rs_sample_get_nsample(s_c)))
    {
        LATAN_ERROR("operation between samples with different numbers of elements",\
                    LATAN_EINVAL);
    }
    
    nsample = rs_sample_get_nsample(s_a);
    status  = LATAN_SUCCESS;
    
    USTAT(f(CENT_VAL(s_a),CENT_VAL(s_b),CENT_VAL(s_c)));
    for (i=0;i<nsample;i++)
    {
        USTAT(f(ELEMENT(s_a,i),ELEMENT(s_b,i),ELEMENT(s_c,i)));
    }
    
    return status;
}

latan_errno rs_sample_binops(rs_sample *s_a, const rs_sample *s_b,\
                             const double s, mat_binops *f)
{
    size_t i;
    size_t nsample;
    latan_errno status;
    
    if (rs_sample_get_nsample(s_a) != rs_sample_get_nsample(s_b))
    {
        LATAN_ERROR("operation between samples with different numbers of elements",\
                    LATAN_EINVAL);
    }
    
    nsample = rs_sample_get_nsample(s_a);
    status  = LATAN_SUCCESS;
    
    USTAT(f(CENT_VAL(s_a),CENT_VAL(s_b),s));
    for (i=0;i<nsample;i++)
    {
        USTAT(f(ELEMENT(s_a,i),ELEMENT(s_b,i),s));
    }
    
    return status;
}

#undef CENT_VAL
#undef ELEMENT
