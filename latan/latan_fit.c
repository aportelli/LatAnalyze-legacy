/* latan_fit.c, part of LatAnalyze library
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

#include <latan/latan_fit.h>
#include <latan/latan_includes.h>
#include <latan/latan_blas.h>
#include <latan/latan_io.h>
#include <latan/latan_math.h>
#include <latan/latan_minimizer.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

/*                       internal functions prototypes                      */
/****************************************************************************/
static size_t sym_rowmaj(const size_t i, const size_t j, const size_t dim);
static size_t get_Ysize(const fit_data *d);
static void init_chi2(fit_data *d, const int thread, const int nthread);


/*                          fit model structure                             */
/****************************************************************************/
/** some useful constant npar_func **/
#define DEFINE_CST_NPAR_FUNC(n)\
size_t npar_##n(void *nothing)\
{\
    nothing = NULL;\
    return n;\
}

DEFINE_CST_NPAR_FUNC(1)
DEFINE_CST_NPAR_FUNC(2)
DEFINE_CST_NPAR_FUNC(3)
DEFINE_CST_NPAR_FUNC(4)
DEFINE_CST_NPAR_FUNC(5)
DEFINE_CST_NPAR_FUNC(6)
DEFINE_CST_NPAR_FUNC(7)
DEFINE_CST_NPAR_FUNC(8)
DEFINE_CST_NPAR_FUNC(9)
DEFINE_CST_NPAR_FUNC(10)

/** access **/
void fit_model_get_name(strbuf name, const fit_model *model)
{
    strbufcpy(name,model->name);
}

size_t fit_model_get_npar(const fit_model *model,void *model_param)
{
    size_t npar;

    npar = model->npar(model_param);

    return npar;
}

double fit_model_eval(const fit_model *model, const mat *x, const mat *p,\
                      void *model_param)
{
    double res;

    res   =  model->func(x,p,model_param);
    
    return res;
}

void fit_model_plot2dstr(strbuf str, const fit_model *model, const size_t k,\
                         const mat *x, const mat *p, void *model_param)
{
    model->plot2dstr(str,k,x,p,model_param);
}

/*                          fit data structure                              */
/****************************************************************************/
static size_t sym_rowmaj(const size_t i, const size_t j, const size_t dim)
{
    size_t si,sj,ind;
    int k;

    si  = (i > j) ? j : i;
    sj  = (i > j) ? i : j;
    ind = sj - si;

    for (k=0;k<=((int)si)-1;k++)
    {
        ind += dim - ((size_t)k);
    }

    return ind;
}

/** allocation **/
fit_data *fit_data_create(const size_t ndata, const size_t ndim)
{
    fit_data *d;
    size_t i,k,k1,k2;
    
    MALLOC_ERRVAL(d,fit_data *,1,NULL);
    MALLOC_ERRVAL(d->to_fit,bool*,ndata,NULL);
    MALLOC_ERRVAL(d->have_x_covar,bool *,ndim,NULL);
    MALLOC_ERRVAL(d->have_xdata_covar,bool *,ndim,NULL);
    d->x            = mat_create(ndim,ndata);
    d->x_covar      = mat_ar_create(ndim*(ndim+1)/2,ndata,ndata);
    d->data         = mat_create(ndata,1);
    d->data_var     = mat_create(ndata,ndata);
    d->data_var_inv = mat_create(ndata,ndata);
    d->xdata_covar  = mat_ar_create(ndim,ndata,ndata);

    d->ndata           = ndata;
    d->ndim            = ndim;
    d->npar            = 0;
    d->x_var_inv       = NULL;
    d->is_x_correlated = false;
    for (k=0;k<ndim;k++)
    {
        d->have_x_covar[k] = false;
    }
    for (i=0;i<ndata;i++)
    {
        d->to_fit[i] = false;
    }
    d->is_data_correlated = false;
    for (k=0;k<ndim;k++)
    {
        d->have_xdata_covar[k] = false;
    }
    d->var_inv       = NULL;
    d->is_inverted   = false;
    d->model         = NULL;
    d->model_param   = NULL;
    d->chi2pdof      = -1.0;
    d->save_chi2pdof = true;
    d->buf           = NULL;
    d->nbuf          = 0;

    for (k1=0;k1<ndim;k1++)
    {
        for (k2=k1;k2<ndim;k2++)
        {
            if (k1 == k2)
            {
                mat_id(d->x_covar[sym_rowmaj(k1,k2,ndim)]);
            }
            else
            {
                mat_zero(d->x_covar[sym_rowmaj(k1,k2,ndim)]);
            }

        }
    }
    mat_id(d->data_var);
    mat_assume_sym(d->data_var,true);
    mat_assume_sym(d->data_var_inv,true);
    for (k=0;k<ndim;k++)
    {
        mat_zero(d->xdata_covar[k]);
    }

    return d;
}

void fit_data_destroy(fit_data *d)
{
    int i;

    mat_destroy(d->x);
    mat_ar_destroy(d->x_covar,d->ndim*(d->ndim+1)/2);
    if (d->x_var_inv != NULL)
    {
        mat_destroy(d->x_var_inv);
    }
    FREE(d->have_x_covar);
    FREE(d->have_xdata_covar);
    mat_destroy(d->data);
    mat_destroy(d->data_var);
    mat_destroy(d->data_var_inv);
    mat_ar_destroy(d->xdata_covar,d->ndim);
    if (d->var_inv != NULL)
    {
        mat_destroy(d->var_inv);
    }
    for(i=0;i<d->nbuf;i++)
    {
        mat_destroy(d->buf[i].x);
        mat_destroy(d->buf[i].X);
        mat_destroy(d->buf[i].CdX);
        if (d->buf[i].is_xpart_alloc)
        {
            mat_destroy(d->buf[i].Y);
            mat_destroy(d->buf[i].CxY);
            mat_destroy(d->buf[i].lX);
            mat_destroy(d->buf[i].ClX);
        }
    }
    FREE(d->buf);
    FREE(d->to_fit);
    FREE(d);
}

/** access **/
/*** sizes ***/
size_t fit_data_get_ndata(const fit_data *d)
{
    return d->ndata;
}

size_t fit_data_get_ndim(const fit_data *d)
{
    return d->ndim;
}

size_t fit_data_get_npar(const fit_data *d)
{
    return d->npar;
}

/*** chi2 value ***/
void fit_data_save_chi2pdof(fit_data *d, bool save)
{
    d->save_chi2pdof = save;
}

double fit_data_get_chi2pdof(const fit_data *d)
{
    return d->chi2pdof;
}

/*** fit points ***/
void fit_data_set_x(fit_data *d, const size_t i, const size_t k,\
                    const double x_ik)
{
    mat_set(d->x,k,i,x_ik);
}

latan_errno fit_data_set_x_vec(fit_data *d, const size_t i, const mat *x_i)
{
    latan_errno status;
    gsl_vector_view x_vview;

    x_vview = gsl_matrix_column(x_i->data_cpu,0);
    status  = gsl_matrix_set_row(d->x->data_cpu,i,&(x_vview.vector));

    return status;
}

double fit_data_get_x(const fit_data *d, const size_t i, const size_t k)
{
    return mat_get(d->x,k,i);
}

latan_errno fit_data_get_x_vec(mat *x_i, const fit_data *d, const size_t i)
{
    latan_errno status;
    gsl_vector_view x_vview;

    x_vview = gsl_matrix_column(x_i->data_cpu,0);
    status  = gsl_matrix_get_row(&(x_vview.vector),d->x->data_cpu,i);
    
    return status;
}

mat *fit_data_pt_x(const fit_data *d)
{
    return d->x;
}

latan_errno fit_data_set_x_covar(fit_data *d, const size_t k1,  \
                                 const size_t k2, const mat *var)
{
    latan_errno status;
    size_t ind;

    ind                 = sym_rowmaj(k1,k2,d->ndim);
    d->have_x_covar[k1] = true;
    d->have_x_covar[k2] = true;
    d->is_x_correlated  = d->is_x_correlated \
                           || (mat_is_square(var) || (k1 != k2));
    d->is_inverted      = false;
    
    if (mat_is_square(var))
    {
        status = mat_cp(d->x_covar[ind],var);
    }
    else
    {
        mat_zero(d->x_covar[ind]);
        status = mat_set_diag(d->x_covar[ind],var);
    }

    return status;
}

const mat * fit_data_pt_x_covar(const fit_data *d, const size_t k1,\
                                const size_t k2)
{
    size_t ind;
    
    ind = sym_rowmaj(k1,k2,d->ndim);
    
    return d->x_covar[ind];
}

bool fit_data_have_x_covar(const fit_data *d, const size_t k)
{
    return d->have_x_covar[k];
}

bool fit_data_have_x_var(const fit_data *d)
{
    size_t k;
    bool have_x_var;

    have_x_var = false;

    for (k=0;k<d->ndim;k++)
    {
        have_x_var = have_x_var || d->have_x_covar[k];
    }

    return have_x_var;
}

bool fit_data_is_x_correlated(const fit_data *d)
{
    return d->is_x_correlated;
}

void fit_data_fit_all_points(fit_data *d, const bool fit)
{
    size_t i;
    
    for (i=0;i<d->ndata;i++)
    {
        d->to_fit[i] = fit;
    }
}

void fit_data_fit_point(fit_data *d, const size_t i, const bool fit)
{
    if (i>=d->ndata)
    {
        LATAN_ERROR_VOID("index out of range",LATAN_EBADLEN);
    }
    
    d->to_fit[i] = fit;
}

void fit_data_fit_range(fit_data *d, const size_t start, const size_t end,\
                        const bool fit)
{
    size_t i;
    
    if ((start >= d->ndata)||(end >= d->ndata))
    {
        LATAN_ERROR_VOID("index out of range",LATAN_EBADLEN);
    }
    if (start > end)
    {
        LATAN_WARNING("trying to modify an empty fit range",LATAN_EBADLEN);
    }
    
    for (i=start;i<=end;i++)
    {
        d->to_fit[i] = fit;
    }
}

void fit_data_fit_region(fit_data *d, double **xb)
{
    size_t i,k;
    bool fit;
    double x_ik;
    
    for (i=0;i<d->ndata;i++)
    {
        fit = true;
        for (k=0;k<d->ndim;k++)
        {
            if (xb[k] != NULL)
            {
                x_ik = fit_data_get_x(d,i,k);
                fit  = fit && (xb[k][0] < x_ik) && (xb[k][1] > x_ik);
            }
        }
        fit_data_fit_point(d,i,fit);
    }
}

bool fit_data_is_fit_point(const fit_data *d, const size_t i)
{
    if (i>=d->ndata)
    {
        LATAN_ERROR_VAL("index out of range",LATAN_EBADLEN,false);
    }
    
    return d->to_fit[i];
}

size_t fit_data_fit_point_num(const fit_data *d)
{
    size_t nfitpt;
    size_t i;
    
    nfitpt = 0;
    
    for (i=0;i<d->ndata;i++)
    {
        if (d->to_fit[i])
        {
            nfitpt++;
        }
    }
    
    return nfitpt;
}

/*** data ***/
void fit_data_set_data(fit_data *d, const size_t i, const double data_i)
{
    mat_set(d->data,i,0,data_i);
}

double fit_data_get_data(const fit_data *d, const size_t i)
{
    return mat_get(d->data,i,0);
}

mat *fit_data_pt_data(const fit_data *d)
{
    return d->data;
}

latan_errno fit_data_set_data_var(fit_data *d, const mat *var)
{
    latan_errno status;
    
    d->is_data_correlated = mat_is_square(var);
    d->is_inverted        = false;
    
    if (mat_is_square(var))
    {
        status = mat_cp(d->data_var,var);
    }
    else
    {
        mat_zero(d->data_var);
        status = mat_set_diag(d->data_var,var);
    }
    
    return status;
}

const mat * fit_data_pt_data_var(const fit_data *d)
{
    return d->data_var;
}

bool fit_data_is_data_correlated(const fit_data *d)
{
    return d->is_data_correlated;
}

latan_errno fit_data_set_xdata_covar(fit_data *d, const size_t k,\
                                     const mat *covar)
{
    latan_errno status;

    d->have_xdata_covar[k] = true;
    d->is_inverted         = false;

    if (mat_is_square(covar))
    {
        status = mat_cp(d->xdata_covar[k],covar);
    }
    else
    {
        mat_zero(d->xdata_covar[k]);
        status = mat_set_diag(d->xdata_covar[k],covar);
    }

    return status;
}

bool fit_data_have_xdata_covar(const fit_data *d)
{
    size_t k;
    bool have_xdata_covar;

    have_xdata_covar = false;

    for (k=0;k<d->ndim;k++)
    {
        have_xdata_covar = have_xdata_covar || d->have_xdata_covar[k];
    }

    return have_xdata_covar;
}

/*** model ***/
latan_errno fit_data_set_model(fit_data *d, const fit_model *model,\
                               void *model_param)
{
    if (model->ndim != d->ndim)
    {
        LATAN_ERROR("fit model and fit data dimension mismatch",LATAN_EBADLEN);
    }
    
    d->model       = model;
    d->model_param = model_param;
    d->npar        = fit_model_get_npar(d->model,d->model_param);

    return LATAN_SUCCESS;
}

double fit_data_model_xeval(const fit_data *d, const mat *x, const mat *p)
{
    return fit_model_eval(d->model,x,p,d->model_param);
}

/* CRITICALLY CALLED FUNCTION
 * fit_data_model_eval is heavily called during one call of minimize function,
 * optimization is done using GSL matrix view
 */
double fit_data_model_eval(const fit_data *d, const size_t i, const mat *p)
{
    mat x_i;
    gsl_matrix_view x_view;
    double res;
    
    x_view       = gsl_matrix_submatrix(d->x->data_cpu,0,i,d->ndim,1);
    x_i.data_cpu = &(x_view.matrix);
    res          = fit_data_model_xeval(d,&x_i,p);
    
    return res;
}

void fit_data_plot2dstr(strbuf str, const fit_data *d, const size_t k,\
                        const mat *x, const mat *p)
{
    fit_model_plot2dstr(str,d->model,k,x,p,d->model_param);
}

/*** dof ***/
size_t fit_data_get_dof(const fit_data *d)
{
    size_t dof;
    
    dof = fit_data_fit_point_num(d) - fit_data_get_npar(d);
    
    return dof;
}

/*                          chi2 function                                   */
/****************************************************************************/
/*** compute size of Y vector for chi2 ***/
static size_t get_Ysize(const fit_data *d)
{
    size_t Ysize;
    size_t k;

    Ysize = 0;
    for (k=0;k<d->ndim;k++)
    {
        if (d->have_xdata_covar[k]||d->have_x_covar[k])
        {
            Ysize += d->ndata;
        }
    }

    return Ysize;
}

/*** (re)allocate chi2 buffers and invert variance matrices ***/
static void init_chi2(fit_data *d, const int thread, const int nthread)
{
    mat *txdata_covar_k,*tx_covar_k1k2;
    int t;
    size_t i,k,k1,k2;
    size_t ind,ndata,ndim,lXsize,Ysize,px_ind,px_ind1,px_ind2;
    double inv;
    bool have_xdata_covar;

#ifdef _OPENMP
    #pragma omp critical
#endif
    {
        txdata_covar_k   = NULL;
        tx_covar_k1k2    = NULL;
        have_xdata_covar = fit_data_have_xdata_covar(d);
        ndata            = d->ndata;
        ndim             = d->ndim;
        Ysize            = get_Ysize(d);
        lXsize           = ndata + Ysize;
        
        /* (re)allocating buffers if necessary */
        if (nthread > d->nbuf)
        {
            REALLOC_NOERRET(d->buf,d->buf,chi2_buf *,nthread);
            for (t=d->nbuf;t<nthread;t++)
            {
                d->buf[t].x   = mat_create(ndim,1);
                d->buf[t].X   = mat_create(ndata,1);
                d->buf[t].CdX = mat_create(ndata,1);
                if (Ysize > 0)
                {
                    d->buf[t].Y              = mat_create(Ysize,1);
                    d->buf[t].CxY            = mat_create(Ysize,1);
                    d->buf[t].lX             = mat_create(lXsize,1);
                    d->buf[t].ClX            = mat_create(lXsize,1);
                    d->buf[t].is_xpart_alloc = true;
                }
                else
                {
                    d->buf[t].Y              = NULL;
                    d->buf[t].CxY            = NULL;
                    d->buf[t].lX             = NULL;
                    d->buf[t].ClX            = NULL;
                    d->buf[t].is_xpart_alloc = false;
                }
            }
            d->nbuf = nthread;
        }
        else if (Ysize > 0)
        {
            if ((d->buf[thread].is_xpart_alloc\
                &&((nrow(d->buf[thread].Y)  != Ysize)\
                ||(nrow(d->buf[thread].CxY) != Ysize) \
                ||(nrow(d->buf[thread].lX)  != lXsize)\
                ||(nrow(d->buf[thread].ClX) != lXsize)))
                ||(!d->buf[thread].is_xpart_alloc))
            {
                if (d->buf[thread].is_xpart_alloc)
                {
                    mat_destroy(d->buf[thread].Y);
                    mat_destroy(d->buf[thread].CxY);
                    mat_destroy(d->buf[thread].lX);
                    mat_destroy(d->buf[thread].ClX);
                }
                d->buf[thread].Y              = mat_create(Ysize,1);
                d->buf[thread].CxY            = mat_create(Ysize,1);
                d->buf[thread].lX             = mat_create(lXsize,1);
                d->buf[thread].ClX            = mat_create(lXsize,1);
                d->buf[thread].is_xpart_alloc = true;
            }
            
        }
        /* (re)allocating global inverse variance matrix if necessary */
        if (have_xdata_covar)
        {
            if (d->var_inv == NULL)
            {
                d->var_inv = mat_create(lXsize,lXsize);
                mat_assume_sym(d->var_inv,true);
            }
            else if (nrow(d->var_inv) != lXsize)
            {
                mat_destroy(d->var_inv);
                d->var_inv = mat_create(lXsize,lXsize);
                mat_assume_sym(d->var_inv,true);
            }
        }
        /* inverting covariance matrices if necessary */
        if (!d->is_inverted)
        {
            /** inversing data variance matrix **/
            if (have_xdata_covar)
            {
                mat_set_subm(d->var_inv,d->data_var,0,0,ndata-1,ndata-1);
            }
            if (fit_data_is_data_correlated(d))
            {
                mat_inv(d->data_var_inv,d->data_var);
            }
            else
            {
                mat_zero(d->data_var_inv);
                for (i=0;i<ndata;i++)
                {
                    inv = 1.0/mat_get(d->data_var,i,i);
                    if (gsl_isinf(inv))
                    {
                        strbuf errmsg;

                        sprintf(errmsg,"errorless data %lu excluded from fit",\
                                (long unsigned)i);
                        LATAN_WARNING(errmsg,LATAN_EDOM);
                        inv = 0.0;
                    }
                    mat_set(d->data_var_inv,i,i,inv);
                }
            }
            latan_printf(DEBUG,"Cd^-1=\n");
            if (latan_get_verb() == DEBUG)
            {
                mat_print(d->data_var_inv,"% 6.1e");
            }
            /** inversing x variance matrix **/
            if (fit_data_have_x_var(d)||have_xdata_covar)
            {
                tx_covar_k1k2 = mat_create(ndata,ndata);

                /*** (re)allocating x inverse variance matrix if necessary ***/
                if (d->x_var_inv == NULL)
                {
                    d->x_var_inv = mat_create(Ysize,Ysize);
                    mat_assume_sym(d->x_var_inv,true);
                }
                else if (nrow(d->x_var_inv) != Ysize)
                {
                    mat_destroy(d->x_var_inv);
                    d->x_var_inv = mat_create(Ysize,Ysize);
                    mat_assume_sym(d->x_var_inv,true);
                }
                mat_zero(d->x_var_inv);
                /*** building x variance matrix by blocks ***/
                px_ind1 = 0;
                for (k1=0;k1<ndim;k1++)
                {
                    if (fit_data_have_x_covar(d,k1))
                    {
                        px_ind2 = px_ind1;
                        for (k2=k1;k2<ndim;k2++)
                        {
                            if (fit_data_have_x_covar(d,k2))
                            {
                                ind = sym_rowmaj(k1,k2,ndim);
                                mat_set_subm(d->x_var_inv,d->x_covar[ind],\
                                             px_ind1*ndata,px_ind2*ndata, \
                                             (px_ind1+1)*ndata-1,         \
                                             (px_ind2+1)*ndata-1);
                                if (k1 != k2)
                                {
                                    mat_transpose(tx_covar_k1k2,\
                                                  d->x_covar[ind]);
                                    mat_set_subm(d->x_var_inv,tx_covar_k1k2,\
                                                 px_ind2*ndata,px_ind1*ndata, \
                                                 (px_ind2+1)*ndata-1,         \
                                                 (px_ind1+1)*ndata-1);
                                }
                                px_ind2++;
                            }
                            
                        }
                        px_ind1++;
                    }
                }
                if (have_xdata_covar)
                {
                    mat_set_subm(d->var_inv,d->x_var_inv,ndata,ndata,\
                                 lXsize-1,lXsize-1);
                }
                /*** inversion ***/
                if (fit_data_is_x_correlated(d))
                {
                    mat_eqinv(d->x_var_inv);
                }
                else
                {
                    for (i=0;i<Ysize;i++)
                    {
                        inv = 1.0/mat_get(d->x_var_inv,i,i);
                        if (gsl_isinf(inv))
                        {
                            LATAN_ERROR_NORET("errorless point found in x variance matrix",\
                                              LATAN_EDOM);
                        }
                        mat_set(d->x_var_inv,i,i,inv);
                    }
                }
                latan_printf(DEBUG,"Cx^-1=\n");
                if (latan_get_verb() == DEBUG)
                {
                    mat_print(d->x_var_inv,"% 6.1e");
                }
                
                mat_destroy(tx_covar_k1k2);
            }
            /** inversing global variance matrix **/
            if (have_xdata_covar)
            {
                txdata_covar_k = mat_create(ndata,ndata);

                px_ind = 0;
                for (k=0;k<ndim;k++)
                {
                    if (d->have_xdata_covar[k])
                    {
                        mat_set_subm(d->var_inv,d->xdata_covar[k],0,\
                                     (px_ind+1)*ndata,ndata-1,      \
                                     (px_ind+2)*ndata-1);
                        mat_transpose(txdata_covar_k,d->xdata_covar[k]);
                        mat_set_subm(d->var_inv,txdata_covar_k,            \
                                     (px_ind+1)*ndata,0,(px_ind+2)*ndata-1,\
                                     ndata-1);
                        px_ind++;
                    }
                }
                mat_eqinv(d->var_inv);
                latan_printf(DEBUG,"C^-1=\n");
                if (latan_get_verb() == DEBUG)
                {
                    mat_print(d->var_inv,"% 6.1e");
                }

                mat_destroy(txdata_covar_k);
            }
            d->is_inverted = true;
        }
    }
}

/* chi^2 function :
 * ----------------
 * 
 * chi^2 is defined this way :
 *
 * chi^2 = t(lX)*C*lX  (t = transposition)
 *
 * where lX is the column vector defined by blocks :
 *
 *      (                    )    i    : index for the data
 *      ( X = f(q_i,p) - d_i )    k    : index for the dimension
 *      (                    )    d    : data
 * lX = ( .................. )    x    : points
 *      (                    )    f    : model
 *      ( Y =  q_ik - x_ik   )    p    : fit parameters
 *      (                    )    q_ik : x_ik if you don't consider x covariance
 *                                       on dimension k, an additional fit
 *                                       parameter else
 *
 * lX have size ndata + Ysize, where Ysize depends on how many covariance
 * relations you have between dimensions (it is given by the get_Ysize
 * function).
 * C is the inverse of the covariance matrix defined by blocks :
 *
 *      (                .                 )  i  : line data index
 *      (    <d_i*d_j>   .   <x_jl*d_i>    )  j  : column data index
 *      (                .                 )  k  : line dimension index
 *      ( ................................ )  l  : column dimension index
 *      (                .                 )
 *      (   <x_ik*d_j>   .  <x_ik*x_jl>    )  for the blocks involving x,
 *      (                .                 )  dimension major order indexing is
 *                                            used
 *
 * CRITICALLY CALLED FUNCTION
 * chi2 is heavily called during one call of minimize function,
 * optimization is done using vector/matrix views from GSL, matrix buffers
 * in fit_data structure, and BLAS operations
 * 
 */
double chi2(mat *p, void *vd)
{
    fit_data *d;
    int nthread,thread;
    size_t ndata,ndim,npar,px_ind;
    size_t i,k;
    mat *x,*X,*CdX,*Cd,*Y,*CxY,*Cx,*lX,*ClX,*C;
    double res,buf,eval;
    
    d       = (fit_data *)vd;
    ndata   = fit_data_get_ndata(d);
    ndim    = fit_data_get_ndim(d);
    npar    = fit_data_get_npar(d);
#ifdef _OPENMP
    nthread = omp_get_num_threads();
    thread  = omp_get_thread_num();
#else
    nthread = 1;
    thread  = 0;
#endif

    /* buffers and inverse variance matrices initialization */
    init_chi2(d,thread,nthread);
    x   = d->buf[thread].x;
    X   = d->buf[thread].X;
    CdX = d->buf[thread].CdX;
    Cd  = d->data_var_inv;
    Y   = d->buf[thread].Y;
    CxY = d->buf[thread].CxY;
    Cx  = d->x_var_inv;
    lX  = d->buf[thread].lX;
    ClX = d->buf[thread].ClX;
    C   = d->var_inv;

    /* setting X and Y in case of covariance in x */
    if (fit_data_have_x_var(d))
    {
        for (i=0;i<ndata;i++)
        {
            px_ind = 0;
            if (fit_data_is_fit_point(d,i))
            {
                for (k=0;k<ndim;k++)
                {
                    if (fit_data_have_x_covar(d,k))
                    {
                        buf = mat_get(p,npar+px_ind*ndata+i,0);
                        mat_set(Y,px_ind*ndata+i,0,buf-fit_data_get_x(d,i,k));
                        px_ind++;
                    }
                    else
                    {
                        buf = fit_data_get_x(d,i,k);
                    }
                    mat_set(x,k,0,buf);
                }
                eval = fit_model_eval(d->model,x,p,d->model_param);
                mat_set(X,i,0,eval - fit_data_get_data(d,i));
            }
            else
            {

                for (k=0;k<ndim;k++)
                {
                    if (fit_data_have_x_covar(d,k))
                    {
                        mat_set(Y,px_ind*ndata+i,0,0.0);
                        px_ind++;
                    }
                }
                mat_set(X,i,0,0.0);
            }
        }
    }
    /* setting X in case of no covariance in x */
    else
    {
        for (i=0;i<ndata;i++)
        {
            if (fit_data_is_fit_point(d,i))
            {
                mat_set(X,i,0,fit_data_model_eval(d,i,p)\
                        -fit_data_get_data(d,i));
            }
            else
            {
                mat_set(X,i,0,0.0);
            }
        }
    }
    /* setting lX and computing chi^2 in case of data/x covariance */
    if (fit_data_have_xdata_covar(d))
    {
        mat_set_subm(lX,X,0,0,d->ndata-1,0);
        mat_set_subm(lX,Y,d->ndata,0,nrow(lX)-1,0);
        mat_mul(ClX,C,'n',lX,'n');
        latan_blas_ddot(ClX,lX,&res);
    }
    /* computing chi^2 by blocks in case of no data/x covariance */
    else
    {
        mat_mul(CdX,Cd,'n',X,'n');
        latan_blas_ddot(CdX,X,&res);
        if (fit_data_have_x_var(d))
        {
            mat_mul(CxY,Cx,'n',Y,'n');
            latan_blas_ddot(CxY,Y,&buf);
            res += buf;
        }
    }
    
    return res;
}

/*                          fit functions                                   */
/****************************************************************************/
latan_errno data_fit(mat *p, fit_data *d)
{
    latan_errno status;
    strbuf cor_status;
    double chi2_min;
    
    strbufcpy(cor_status,"correlations :");
    if (fit_data_is_data_correlated(d))
    {
        strcat(cor_status," data/data");
    }
    if (fit_data_is_x_correlated(d))
    {
        strcat(cor_status," x/x");
    }
    if (fit_data_have_xdata_covar(d))
    {
        strcat(cor_status," data/x");
    }
    if (strcmp(cor_status,"correlations :") == 0)
    {
        strcat(cor_status," no");
    }
    latan_printf(VERB,"fitting (%s) %u data points with %s model...\n",
                 cor_status,(unsigned int)fit_data_fit_point_num(d),\
                 d->model->name);
    status = minimize(p,&chi2_min,&chi2,d);
    if (d->save_chi2pdof)
    {
        d->chi2pdof = DRATIO(chi2_min,fit_data_get_dof(d));
    }
    
    return status;
}

latan_errno rs_data_fit(rs_sample *p, const rs_sample *data, fit_data *d,\
                        const cor_flag flag)
{
    latan_errno status;
    mat *datavar;
    size_t ndata;
    size_t i;
    bool is_data_cor;
    int verb_backup;
    double chi2pdof_backup;
        
    verb_backup = latan_get_verb();
    status      = LATAN_SUCCESS;
    ndata       = fit_data_get_ndata(d);
    is_data_cor = ((flag & DATA_COR) != 0);
    
    datavar     = mat_create(ndata,is_data_cor ? ndata : 1);
    
    if (is_data_cor)
    {
        USTAT(rs_sample_var(datavar,data));
    }
    else
    {
        USTAT(rs_sample_varp(datavar,data));
    }
    USTAT(fit_data_set_data_var(d,datavar));
    USTAT(mat_cp(fit_data_pt_data(d),rs_sample_pt_cent_val(data)));
    USTAT(data_fit(rs_sample_pt_cent_val(p),d));
    chi2pdof_backup = fit_data_get_chi2pdof(d);
    latan_printf(DEBUG,"central value chi^2/dof = %e\n",chi2pdof_backup);
    if (verb_backup != DEBUG)
    {
        latan_set_verb(QUIET);
    }
    for (i=0;i<rs_sample_get_nsample(data);i++)
    {
        USTAT(mat_cp(rs_sample_pt_sample(p,i),rs_sample_pt_cent_val(p)));
        USTAT(mat_cp(fit_data_pt_data(d),rs_sample_pt_sample(data,i)));
        USTAT(data_fit(rs_sample_pt_sample(p,i),d));
        latan_printf(DEBUG,"sample %lu chi^2/dof = %e\n",(long unsigned)i,\
                     fit_data_get_chi2pdof(d));
    }
    d->chi2pdof = chi2pdof_backup;
    latan_set_verb(verb_backup);
    
    mat_destroy(datavar);
    
    return status;
}

latan_errno rs_x_data_fit(rs_sample *p, rs_sample * const *x,         \
                          const rs_sample *data, fit_data *d,         \
                          const cor_flag flag, const bool *use_x_var)
{
    latan_errno status;
    mat *datavar,*xcovar,*xdatacovar,*pbuf;
    size_t ndata,ndim,npar,Ysize,px_ind;
    size_t i,k,k1,k2;
    bool is_data_cor,is_x_cor,is_xdata_cor;
    int verb_backup;
    double chi2pdof_backup;
    
    verb_backup  = latan_get_verb();
    status       = LATAN_SUCCESS;
    ndata        = fit_data_get_ndata(d);
    ndim         = fit_data_get_ndim(d);
    npar         = fit_data_get_npar(d);
    is_data_cor  = ((flag & DATA_COR) != 0);
    is_x_cor     = ((flag & X_COR) != 0);
    is_xdata_cor = ((flag & XDATA_COR) != 0);
    
    datavar      = mat_create(ndata,is_data_cor ? ndata : 1);
    xcovar       = mat_create(ndata,(is_x_cor&&is_data_cor) ? ndata : 1);
    xdatacovar   = mat_create(ndata,is_data_cor ? ndata : 1);

    /* compute needed variances/covariances from samples */
    if (is_data_cor)
    {
        USTAT(rs_sample_var(datavar,data));
    }
    else
    {
        USTAT(rs_sample_varp(datavar,data));
    }
    USTAT(fit_data_set_data_var(d,datavar));
    if (is_x_cor)
    {
        if (is_data_cor)
        {
            for (k1=0;k1<ndim;k1++)
            {
                for (k2=k1;k2<ndim;k2++)
                {
                    if (use_x_var[k1]&&use_x_var[k2])
                    {
                        USTAT(rs_sample_cov(xcovar,x[k1],x[k2]));
                        USTAT(fit_data_set_x_covar(d,k1,k2,xcovar));
                    }
                }
            }
        }
        else
        {
            for (k1=0;k1<ndim;k1++)
            {
                for (k2=k1;k2<ndim;k2++)
                {
                    if (use_x_var[k1]&&use_x_var[k2])
                    {
                        USTAT(rs_sample_covp(xcovar,x[k1],x[k2]));
                        USTAT(fit_data_set_x_covar(d,k1,k2,xcovar));
                    }
                }
            }
        }
    }
    else
    {
        for (k=0;k<ndim;k++)
        {
            if (use_x_var[k])
            {
                USTAT(rs_sample_varp(xcovar,x[k]));
                USTAT(fit_data_set_x_covar(d,k,k,xcovar));
            }
        }
    }
    if (is_xdata_cor)
    {
        if (is_data_cor)
        {
            for (k=0;k<ndim;k++)
            {
                if (use_x_var[k])
                {
                    USTAT(rs_sample_cov(xdatacovar,data,x[k]));
                    USTAT(fit_data_set_xdata_covar(d,k,xdatacovar));
                }
            }
        }
        else
        {
            for (k=0;k<ndim;k++)
            {
                if (use_x_var[k])
                {
                    USTAT(rs_sample_covp(xdatacovar,data,x[k]));
                    USTAT(fit_data_set_xdata_covar(d,k,xdatacovar));
                }
            }
        }
    }
    Ysize = get_Ysize(d);
    pbuf  = mat_create(npar+Ysize,1);
    /* central value fit */
    /** setting data and initial parameters **/
    USTAT(mat_cp(fit_data_pt_data(d),rs_sample_pt_cent_val(data)));
    USTAT(mat_set_subm(pbuf,rs_sample_pt_cent_val(p),0,0,npar-1,0));
    px_ind = 0;
    for (k=0;k<ndim;k++)
    {
        USTAT(fit_data_set_x_vec(d,k,rs_sample_pt_cent_val(x[k])));
        if (use_x_var[k])
        {
            USTAT(mat_set_subm(pbuf,rs_sample_pt_cent_val(x[k]),             \
                               npar+px_ind*ndata,0,npar+(px_ind+1)*ndata-1,0));
            px_ind++;
        }
    }
    /** looking at starting chi^2/dof **/
    latan_printf(VERB,"starting chi^2/dof = %f\n",         \
                 chi2(pbuf,d)/((double)fit_data_get_dof(d)));
    /**  fit **/
    USTAT(data_fit(pbuf,d));
    USTAT(mat_get_subm(rs_sample_pt_cent_val(p),pbuf,0,0,npar-1,0));
    chi2pdof_backup = fit_data_get_chi2pdof(d);
    latan_printf(VERB,"fit: central value chi^2/dof = %e\n",chi2pdof_backup);
    
    /* sample fits */
    for (i=0;i<rs_sample_get_nsample(data);i++)
    {
        /** setting data and initial parameters **/
        USTAT(mat_cp(fit_data_pt_data(d),rs_sample_pt_sample(data,i)));
        USTAT(mat_set_subm(pbuf,rs_sample_pt_cent_val(p),0,0,npar-1,0));
        px_ind = 0;
        for (k=0;k<ndim;k++)
        {
            USTAT(fit_data_set_x_vec(d,k,rs_sample_pt_sample(x[k],i)));
            if (use_x_var[k])
            {
                USTAT(mat_set_subm(pbuf,rs_sample_pt_sample(x[k],i),\
                                   npar+px_ind*ndata,0,             \
                                   npar+(px_ind+1)*ndata-1,0));
                px_ind++;
            }
        }
        /** fit **/
        if (verb_backup != DEBUG)
        {
            USTAT(latan_set_verb(QUIET));
        }
        USTAT(data_fit(pbuf,d));
        USTAT(latan_set_verb(verb_backup));
        latan_printf(VERB,"fit: sample %d/%d chi^2/dof = %e\n",(int)i+1,\
                        (int)rs_sample_get_nsample(data),               \
                        fit_data_get_chi2pdof(d));
        USTAT(mat_get_subm(rs_sample_pt_sample(p,i),pbuf,0,0,npar-1,0));
    }
    d->chi2pdof = chi2pdof_backup;
    USTAT(mat_cp(fit_data_pt_data(d),rs_sample_pt_cent_val(data)));
    for (k=0;k<ndim;k++)
    {
        USTAT(fit_data_set_x_vec(d,k,rs_sample_pt_cent_val(x[k])));
    }
    
    mat_destroy(datavar);
    mat_destroy(xcovar);
    mat_destroy(xdatacovar);
    mat_destroy(pbuf);
    
    return status;
}

void fit_residual(mat *res, const mat *p, const fit_data *d)
{
    size_t i,ndata;
    double res_i;

    ndata = fit_data_get_ndata(d);
    
    for (i=0;i<ndata;i++)
    {
        res_i = fit_data_get_data(d,i) - fit_data_model_eval(d,i,p);
        mat_set(res,i,0,res_i);
    }
}

void fit_partresidual(mat *res, const mat *p, const fit_data *d,\
                      const mat *x_ex, const size_t k)
{
    size_t i,ndata,ndim;
    double res_i;
    mat *x_buf;
    
    ndim  = fit_data_get_ndim(d);
    
    if (ndim > 1)
    {
        ndata = fit_data_get_ndata(d);
        x_buf = mat_create(ndim,1);
        for (i=0;i<ndata;i++)
        {
            mat_cp(x_buf,x_ex);
            mat_set(x_buf,k,0,fit_data_get_x(d,i,k));
            res_i = fit_data_get_data(d,i) - fit_data_model_eval(d,i,p)\
                    + fit_data_model_xeval(d,x_buf,p);
            mat_set(res,i,0,res_i);
        }
        mat_destroy(x_buf);
    }
    else
    {
        mat_cp(res,fit_data_pt_data(d));
    }
}
