/* latan_fit.c, part of LatAnalyze library
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
static size_t rowmaj(const size_t i, const size_t j, const size_t dim1,\
                     const size_t dim2);
static size_t sym_rowmaj(const size_t i, const size_t j, const size_t dim);
static double zero(const mat *p, void *vd);
static size_t get_Xsize(const fit_data *d);
static size_t get_Ysize(const fit_data *d);
static void set_var_subm(mat *var, const mat *subm, const size_t k1,\
                         const size_t k2, const fit_data *d);
static void pseudoinvert_var(mat *var, const bool is_corr);
static void init_chi2(fit_data *d, const int thread, const int nthread);
static void set_X_Y(mat* X, mat *Y, mat *x_buf, const mat *p, \
                    const fit_data *d);
static double chi2_base(const mat *p, void *vd);

/*                          fit model structure                             */
/****************************************************************************/
/** some useful constant npar_func **/
#define DEFINE_CST_NPAR_FUNC(n)\
size_t npar_##n(void *nothing __dumb)\
{\
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

double fit_model_eval(const fit_model *model, const size_t k, const mat *x,\
                      const mat *p, void *model_param)
{
    double res;

    res   =  model->func[k](x,p,model_param);
    
    return res;
}

/*                          fit data structure                              */
/****************************************************************************/
static size_t rowmaj(const size_t i, const size_t j, const size_t dim1,\
                     const size_t dim2)
{
    size_t ind,x[2],dim[2];
    
    x[0]   = i;
    x[1]   = j;
    dim[0] = dim1;
    dim[1] = dim2;
    ind    = coord_to_rowmaj(x,dim,2);
    
    return ind;
}

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

static double zero(const mat *p __dumb, void *vd __dumb)
{    
    return 0.0;
}

/** allocation **/
fit_data *fit_data_create(const size_t ndata, const size_t nxdim,\
                          const size_t nydim)
{
    fit_data *d;
    size_t i,k,k1,k2;
    
    MALLOC_ERRVAL(d,fit_data *,1,NULL);
    MALLOC_ERRVAL(d->to_fit,bool*,ndata,NULL);
    MALLOC_ERRVAL(d->have_x_covar,bool *,nxdim,NULL);
    MALLOC_ERRVAL(d->have_xy_covar,bool *,nxdim,NULL);
    d->x          = mat_create(nxdim,ndata);
    d->x_covar    = mat_ar_create(nxdim*(nxdim+1)/2,ndata,ndata);
    d->y          = mat_create(nydim,ndata);
    d->y_covar    = mat_ar_create(nydim*(nydim+1)/2,ndata,ndata);
    d->xy_covar   = mat_ar_create(nxdim*nydim,ndata,ndata);
    d->cor_filter = mat_create(ndata,ndata);

    d->ndata           = ndata;
    d->nxdim           = nxdim;
    d->nydim           = nydim;
    d->model           = NULL;
    d->model_param     = NULL;
    d->npar            = 0;
    d->ndpar           = 0;
    d->x_var_inv       = NULL;
    d->y_var_inv       = NULL;
    d->is_x_correlated = false;
    for (k=0;k<nxdim;k++)
    {
        d->have_x_covar[k] = false;
    }
    for (i=0;i<ndata;i++)
    {
        d->to_fit[i] = false;
    }
    d->is_y_correlated = false;
    for (k=0;k<nxdim;k++)
    {
        d->have_xy_covar[k] = false;
    }
    d->var_inv       = NULL;
    d->is_inverted   = false;
    d->chi2_ext      = &zero;
    d->chi2_val      = -1.0;
    d->chi2_comp     = NULL;
    d->save_chi2pdof = true;
    d->buf           = NULL;
    d->nbuf          = 0;
    d->s             = 0;

    for (k1=0;k1<nxdim;k1++)
    {
        for (k2=k1;k2<nxdim;k2++)
        {
            if (k1 == k2)
            {
                mat_id(d->x_covar[sym_rowmaj(k1,k2,nxdim)]);
            }
            else
            {
                mat_zero(d->x_covar[sym_rowmaj(k1,k2,nxdim)]);
            }

        }
    }
    for (k1=0;k1<nydim;k1++)
    {
        for (k2=k1;k2<nydim;k2++)
        {
            if (k1 == k2)
            {
                mat_id(d->y_covar[sym_rowmaj(k1,k2,nydim)]);
            }
            else
            {
                mat_zero(d->y_covar[sym_rowmaj(k1,k2,nydim)]);
            }
        }
    }
    for (k1=0;k1<nydim;k1++)
    for (k2=0;k2<nxdim;k2++)
    {
        mat_zero(d->xy_covar[rowmaj(k1,k2,nydim,nxdim)]);
    }
    mat_cst(d->cor_filter,1.0);

    return d;
}

void fit_data_destroy(fit_data *d)
{
    int i;

    mat_destroy(d->x);
    mat_ar_destroy(d->x_covar,d->nxdim*(d->nxdim+1)/2);
    FREE(d->have_x_covar);
    FREE(d->have_xy_covar);
    mat_destroy(d->y);
    mat_ar_destroy(d->y_covar,d->nydim*(d->nydim+1)/2);
    mat_ar_destroy(d->xy_covar,d->nydim*d->nxdim);
    mat_destroy(d->cor_filter);
    if (d->x_var_inv != NULL)
    {
        mat_destroy(d->x_var_inv);
    }
    if (d->y_var_inv != NULL)
    {
        mat_destroy(d->y_var_inv);
    }
    if (d->var_inv != NULL)
    {
        mat_destroy(d->var_inv);
    }
    if (d->chi2_comp != NULL)
    {
        mat_destroy(d->chi2_comp);
    }
    for(i=0;i<d->nbuf;i++)
    {
        mat_destroy(d->buf[i].x_f);
        mat_destroy(d->buf[i].Y);
        mat_destroy(d->buf[i].CyY);
        if (d->buf[i].is_xpart_alloc)
        {
            mat_destroy(d->buf[i].X);
            mat_destroy(d->buf[i].CxX);
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

size_t fit_data_get_nydim(const fit_data *d)
{
    return d->nydim;
}

size_t fit_data_get_nxdim(const fit_data *d)
{
    return d->nxdim;
}

size_t fit_data_get_npar(const fit_data *d)
{
    return d->npar;
}

size_t fit_data_get_ndumbpar(const fit_data *d)
{
    return d->ndpar;
}

void fit_data_set_ndumbpar(fit_data *d, const size_t ndpar)
{
    d->ndpar = ndpar;
}

/*** chi^2 ***/
void fit_data_save_chi2(fit_data *d, bool save)
{
    d->save_chi2pdof = save;
}

double fit_data_get_chi2(const fit_data *d)
{
    return d->chi2_val;
}

double fit_data_get_chi2pdof(const fit_data *d)
{
    return d->chi2_val/((double)fit_data_get_dof(d));
}

void fit_data_set_chi2_ext(fit_data *d, min_func *f)
{
    d->chi2_ext = f;
}

latan_errno fit_data_get_chi2_comp(mat *comp, const fit_data *d)
{
    return mat_cp(comp,d->chi2_comp);
}

double fit_data_get_pvalue(const fit_data *d)
{
    return chi2_pvalue(fit_data_get_chi2(d),fit_data_get_dof(d));
}

/*** data ***/
void fit_data_set_y(fit_data *d, const size_t i, const size_t k,\
                    const double y_ik)
{
    mat_set(d->y,k,i,y_ik);
}

latan_errno fit_data_set_y_k(fit_data *d, const size_t k, const mat *y_k)
{
    latan_errno status;
    gsl_vector_view y_vview;
    
    y_vview = gsl_matrix_column(y_k->data_cpu,0);
    status  = gsl_matrix_set_row(d->y->data_cpu,k,&(y_vview.vector));
    
    return status;
}

double fit_data_get_y(const fit_data *d, const size_t i, const size_t k)
{
    return mat_get(d->y,k,i);
}

latan_errno fit_data_get_y_k(mat *y_k, const fit_data *d, const size_t k)
{
    latan_errno status;
    gsl_vector_view y_vview;
    
    y_vview = gsl_matrix_column(y_k->data_cpu,0);
    status  = gsl_matrix_get_row(&(y_vview.vector),d->y->data_cpu,k);
    
    return status;
}

mat *fit_data_pt_y(const fit_data *d)
{
    return d->y;
}

latan_errno fit_data_set_y_covar(fit_data *d, const size_t k1,  \
                                 const size_t k2, const mat *var)
{
    latan_errno status;
    size_t ind;
    
    ind                 = sym_rowmaj(k1,k2,d->nydim);
    d->is_y_correlated  = d->is_y_correlated \
                          || (mat_is_square(var) || (k1 != k2));
    d->is_inverted      = false;
    
    if (mat_is_square(var))
    {
        status = mat_cp(d->y_covar[ind],var);
    }
    else
    {
        mat_zero(d->y_covar[ind]);
        status = mat_set_diag(d->y_covar[ind],var);
    }
    
    return status;
}

const mat * fit_data_pt_y_covar(const fit_data *d, const size_t k1,\
                                const size_t k2)
{
    size_t ind;
    
    ind = sym_rowmaj(k1,k2,d->nydim);
    
    return d->y_covar[ind];
}

bool fit_data_is_y_correlated(const fit_data *d)
{
    return d->is_y_correlated;
}

latan_errno fit_data_set_xy_covar(fit_data *d, const size_t ky, \
                                  const size_t kx, const mat *covar)
{
    latan_errno status;
    size_t ind;
    
    d->have_xy_covar[kx] = true;
    d->is_inverted       = false;
    ind                  = rowmaj(ky,kx,d->nydim,d->nxdim);
    
    if (mat_is_square(covar))
    {
        status = mat_cp(d->xy_covar[ind],covar);
    }
    else
    {
        mat_zero(d->xy_covar[ind]);
        status = mat_set_diag(d->xy_covar[ind],covar);
    }
    
    return status;
}

bool fit_data_have_xy_covar(const fit_data *d)
{
    size_t k;
    bool have_xy_covar;
    
    have_xy_covar = false;
    
    for (k=0;k<d->nxdim;k++)
    {
        have_xy_covar = have_xy_covar || d->have_xy_covar[k];
    }
    
    return have_xy_covar;
}

/*** fit points ***/
void fit_data_set_x(fit_data *d, const size_t i, const size_t k,\
                    const double x_ik)
{
    mat_set(d->x,k,i,x_ik);
}

latan_errno fit_data_set_x_k(fit_data *d, const size_t k, const mat *x_k)
{
    latan_errno status;
    gsl_vector_view x_vview;

    x_vview = gsl_matrix_column(x_k->data_cpu,0);
    status  = gsl_matrix_set_row(d->x->data_cpu,k,&(x_vview.vector));

    return status;
}

double fit_data_get_x(const fit_data *d, const size_t i, const size_t k)
{
    return mat_get(d->x,k,i);
}

latan_errno fit_data_get_x_k(mat *x_k, const fit_data *d, const size_t k)
{
    latan_errno status;
    gsl_vector_view x_vview;

    x_vview = gsl_matrix_column(x_k->data_cpu,0);
    status  = gsl_matrix_get_row(&(x_vview.vector),d->x->data_cpu,k);
    
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

    ind                 = sym_rowmaj(k1,k2,d->nxdim);
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
    
    ind = sym_rowmaj(k1,k2,d->nxdim);
    
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

    for (k=0;k<d->nxdim;k++)
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
        for (k=0;k<d->nxdim;k++)
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

latan_errno fit_data_set_data_cor(fit_data *d, const size_t i, const size_t j,\
                                  bool is_cor)
{
    if (i == j) 
    {
        LATAN_ERROR("cannot uncorrelate a point with itself",LATAN_EINVAL);
    }
    
    d->is_inverted      = false;   
    if (is_cor)
    {
        mat_set(d->cor_filter,i,j,1.0);
        mat_set(d->cor_filter,j,i,1.0);
    }
    else
    {
        mat_set(d->cor_filter,i,j,0.0);
        mat_set(d->cor_filter,j,i,0.0);
    }
    
    return LATAN_SUCCESS;
}

bool fit_data_is_data_cor(const fit_data *d, const size_t i, const size_t j)
{
    return (mat_get(d->cor_filter,i,j) >= 0.5);
}

/*** model ***/
latan_errno fit_data_set_model(fit_data *d, fit_model *model,\
                               void *model_param)
{
    if ((model->nxdim != d->nxdim)||(model->nydim != d->nydim))
    {
        LATAN_ERROR("fit model and fit data dimension mismatch",LATAN_EBADLEN);
    }
    
    d->model       = model;
    d->model_param = model_param;
    d->npar        = fit_model_get_npar(d->model,d->model_param);

    return LATAN_SUCCESS;
}

const void * fit_data_pt_model_param(const fit_data *d)
{
    return d->model_param;
}

double fit_data_model_xeval(const fit_data *d, const size_t k, const mat *x,\
                            const mat *p)
{
    return fit_model_eval(d->model,k,x,p,d->model_param);
}

/* CRITICALLY CALLED FUNCTION
 * fit_data_model_eval is heavily called during one call of minimize function,
 * optimization is done using GSL matrix view
 */
double fit_data_model_eval(const fit_data *d, const size_t k, const size_t i,\
                           const mat *p)
{
    mat x_i;
    gsl_matrix_view x_view;
    double res;
    
    x_view       = gsl_matrix_submatrix(d->x->data_cpu,0,i,d->nxdim,1);
    x_i.data_cpu = &(x_view.matrix);
    res          = fit_data_model_xeval(d,k,&x_i,p);
    
    return res;
}

/*** dof ***/
size_t fit_data_get_dof(const fit_data *d)
{
    size_t dof;
    
    dof = fit_data_fit_point_num(d)*fit_data_get_nydim(d)\
          - fit_data_get_npar(d) + fit_data_get_ndumbpar(d);
    
    return dof;
}

/*** sample counter ***/
size_t fit_data_get_sample_counter(const fit_data *d)
{
    return d->s;
}

/*** set from samples ***/
latan_errno fit_data_set_covar_from_sample(fit_data *d, rs_sample * const *x,\
                                           rs_sample * const *data,          \
                                           const cor_flag flag,              \
                                           const bool *use_x_var)
{
    latan_errno status;
    mat *datavar,*xcovar,*xdatacovar;
    size_t ndata,nxdim,nydim;
    size_t k,k1,k2;
    bool is_data_cor,is_x_cor,is_xdata_cor;
    
    status       = LATAN_SUCCESS;
    ndata        = fit_data_get_ndata(d);
    nxdim        = fit_data_get_nxdim(d);
    nydim        = fit_data_get_nydim(d);
    is_data_cor  = ((flag & DATA_COR) != 0);
    is_x_cor     = ((flag & X_COR) != 0)&&(x != NULL);
    is_xdata_cor = ((flag & XDATA_COR) != 0)&&(x != NULL);
    
    datavar      = mat_create(ndata,is_data_cor ? ndata : 1);
    xcovar       = mat_create(ndata,(is_x_cor&&is_data_cor) ? ndata : 1);
    xdatacovar   = mat_create(ndata,is_data_cor ? ndata : 1);
    
    /* compute needed variances/covariances from samples */
    /** data **/
    for (k1=0;k1<nydim;k1++)
    for (k2=k1;k2<nydim;k2++)
    {
        if (is_data_cor)
        {
            USTAT(rs_sample_cov(datavar,data[k1],data[k2]));
            USTAT(fit_data_set_y_covar(d,k1,k2,datavar));
        }
        else if (k1 == k2)
        {
            USTAT(rs_sample_covp(datavar,data[k1],data[k2]));
            USTAT(fit_data_set_y_covar(d,k1,k2,datavar));
        }
    }
    for (k=0;k<nydim;k++)
    {
        USTAT(fit_data_set_y_k(d,k,rs_sample_pt_cent_val(data[k])));
    }
    /** points **/
    if (is_x_cor)
    {
        for (k1=0;k1<nxdim;k1++)
        for (k2=k1;k2<nxdim;k2++)
        {
            if (use_x_var[k1]&&use_x_var[k2])
            {
                if (is_data_cor)
                {
                    USTAT(rs_sample_cov(xcovar,x[k1],x[k2]));
                }
                else
                {
                    USTAT(rs_sample_covp(xcovar,x[k1],x[k2]));
                }
                USTAT(fit_data_set_x_covar(d,k1,k2,xcovar));
            }
        }
    }
    else if (x != NULL)
    {
        for (k=0;k<nxdim;k++)
        {
            if (use_x_var[k])
            {
                USTAT(rs_sample_varp(xcovar,x[k]));
                USTAT(fit_data_set_x_covar(d,k,k,xcovar));
            }
        }
    }
    /** data/points **/
    if (is_xdata_cor)
    {
        for (k1=0;k1<nydim;k1++)
        for (k2=0;k2<nxdim;k2++)
        {
            if (use_x_var[k2])
            {
                if (is_data_cor)
                {
                    USTAT(rs_sample_cov(xdatacovar,data[k1],x[k2]));
                }
                else
                {
                    USTAT(rs_sample_covp(xdatacovar,data[k1],x[k2]));
                }
                USTAT(fit_data_set_xy_covar(d,k1,k2,xdatacovar));
            }
        }
    }
    
    /* set data/points to central values */
    for (k=0;k<nydim;k++)
    {
        USTAT(fit_data_set_y_k(d,k,rs_sample_pt_cent_val(data[k])));
    }
    if (x != NULL)
    {
        for (k=0;k<nxdim;k++)
        {
            USTAT(fit_data_set_x_k(d,k,rs_sample_pt_cent_val(x[k])));
        }
    }
    
    mat_destroy(datavar);
    mat_destroy(xcovar);
    mat_destroy(xdatacovar);
    
    return status;
}

/*                          chi2 function                                   */
/****************************************************************************/
/** compute size of vectors for chi2 ***/
static size_t get_Xsize(const fit_data *d)
{
    size_t Xsize,npt;
    size_t k;

    Xsize = 0;
    npt   = fit_data_fit_point_num(d);
    for (k=0;k<d->nxdim;k++)
    {
        if (d->have_xy_covar[k]||d->have_x_covar[k])
        {
            Xsize += npt;
        }
    }

    return Xsize;
}

static size_t get_Ysize(const fit_data *d)
{
    return d->nydim*fit_data_fit_point_num(d);
}

/** set sub-matrix in variance matrix **/
static void set_var_subm(mat *var, const mat *subm, const size_t k1,\
                         const size_t k2, const fit_data *d)
{
    size_t ndata,npt;
    size_t i,j,k_i,l_j;
    
    ndata = fit_data_get_ndata(d);
    npt   = fit_data_fit_point_num(d);
    k_i = 0;
    for (i=0;i<ndata;i++)
    {
        if (fit_data_is_fit_point(d,i))
        {
            l_j = 0;
            for (j=0;j<ndata;j++)
            {
                if (fit_data_is_fit_point(d,j))
                {
                    mat_set(var,k1*npt+k_i,k2*npt+l_j,mat_get(subm,i,j));
                    l_j++;
                }
            }
            k_i++;
        }
    }
}

#define mat_debug(m,fname)\
{\
    FILE *f_;\
    f_=fopen(fname,"w");\
    mat_dump(f_,m,"% .15e");\
    fclose(f_);\
}

/** variance matrix inversion **/
static void pseudoinvert_var(mat *var, const bool is_corr)
{
    mat *sig;
    size_t i,j;
    double inv,fac_ij;
    
    if (is_corr)
    {
        sig = mat_create_from_dim(var);
    }
    else
    {
        sig = var;
    }
    for (i=0;i<nrow(var);i++)
    {
        inv = 1.0/mat_get(var,i,i);
        if (is_corr)
        {
            inv = sqrt(inv);
        }
        if (gsl_isinf(inv))
        {
            strbuf errmsg;
            sprintf(errmsg,"errorless data %lu excluded from fit",\
                    (long unsigned)i);
            LATAN_WARNING(errmsg,LATAN_EDOM);
            inv = 0.0;
        }
        mat_set(sig,i,i,inv);
    }
    if (is_corr)
    {
        FOR_VAL(var,i,j)
        {
            fac_ij = mat_get(sig,i,i)*mat_get(sig,j,j);
            mat_set(var,i,j,mat_get(var,i,j)*fac_ij);
        }
        mat_eqpseudoinv(var);
        FOR_VAL(var,i,j)
        {
            fac_ij = mat_get(sig,i,i)*mat_get(sig,j,j);
            mat_set(var,i,j,mat_get(var,i,j)*fac_ij);
        }
    }
    
    if (is_corr)
    {
        mat_destroy(sig);
    }
}
    
/* chi^2 function :
 * ----------------
 * 
 * chi^2 is defined this way :
 *
 * chi^2 = t(lX)*C^-1*lX  (t = transposition)
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
 * C is the covariance matrix defined by blocks :
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

/* (re)allocate chi2 buffers and compute C^-1 */
static void init_chi2(fit_data *d, const int thread, const int nthread)
{
    mat *tmp_covar;
    int t;
    size_t k1,k2;
    size_t ind,ndata,nydim,nxdim,lXsize,Ysize,Xsize,px_ind,px_ind1,px_ind2;
    bool have_xy_covar;
    
#ifdef _OPENMP
#pragma omp critical
#endif
    {
        tmp_covar      = NULL;
        have_xy_covar  = fit_data_have_xy_covar(d);
        ndata          = fit_data_get_ndata(d);
        nydim          = fit_data_get_nydim(d);
        nxdim          = fit_data_get_nxdim(d);
        Ysize          = get_Ysize(d);
        Xsize          = get_Xsize(d);
        lXsize         = Ysize + Xsize;
        
        /* (re)allocating buffers if necessary */
        if (nthread > d->nbuf)
        {
            REALLOC_NOERRET(d->buf,d->buf,chi2_buf *,nthread);
            for (t=d->nbuf;t<nthread;t++)
            {
                d->buf[t].x_f = mat_create(nxdim,1);
                d->buf[t].Y   = mat_create(Ysize,1);
                d->buf[t].CyY = mat_create(Ysize,1);
                if (Xsize > 0)
                {
                    d->buf[t].X              = mat_create(Xsize,1);
                    d->buf[t].CxX            = mat_create(Xsize,1);
                    d->buf[t].lX             = mat_create(lXsize,1);
                    d->buf[t].ClX            = mat_create(lXsize,1);
                    d->buf[t].is_xpart_alloc = true;
                }
                else
                {
                    d->buf[t].X              = NULL;
                    d->buf[t].CxX            = NULL;
                    d->buf[t].lX             = NULL;
                    d->buf[t].ClX            = NULL;
                    d->buf[t].is_xpart_alloc = false;
                }
            }
            d->nbuf = nthread;
        }
        else
        {
            if ((nrow(d->buf[thread].Y)     != Ysize)\
                ||(nrow(d->buf[thread].CyY) != Ysize))
            {
                mat_destroy(d->buf[thread].Y);
                mat_destroy(d->buf[thread].CyY);
                d->buf[thread].Y   = mat_create(Ysize,1);
                d->buf[thread].CyY = mat_create(Ysize,1);
            }
            if (Xsize > 0)
            {
                if ((d->buf[thread].is_xpart_alloc
                      &&((nrow(d->buf[thread].X)    != Xsize)  \
                        ||(nrow(d->buf[thread].CxX) != Xsize)  \
                        ||(nrow(d->buf[thread].lX)  != lXsize) \
                        ||(nrow(d->buf[thread].ClX) != lXsize)))
                    ||(!d->buf[thread].is_xpart_alloc))
                {
                    if (d->buf[thread].is_xpart_alloc)
                    {
                        mat_destroy(d->buf[thread].X);
                        mat_destroy(d->buf[thread].CxX);
                        mat_destroy(d->buf[thread].lX);
                        mat_destroy(d->buf[thread].ClX);
                    }
                    d->buf[thread].X              = mat_create(Xsize,1);
                    d->buf[thread].CxX            = mat_create(Xsize,1);
                    d->buf[thread].lX             = mat_create(lXsize,1);
                    d->buf[thread].ClX            = mat_create(lXsize,1);
                    d->buf[thread].is_xpart_alloc = true;
                }
            }
        }
        
        /* (re)allocating global inverse variance matrix if necessary */
        if (have_xy_covar)
        {
            if (d->var_inv == NULL)
            {
                d->var_inv = mat_create(lXsize,lXsize);
                mat_assume(d->var_inv,MAT_SYM|MAT_POS);
                mat_zero(d->var_inv);
                d->is_inverted = false;
            }
            else if (nrow(d->var_inv) != lXsize)
            {
                mat_destroy(d->var_inv);
                d->var_inv = mat_create(lXsize,lXsize);
                mat_assume(d->var_inv,MAT_SYM|MAT_POS);
                mat_zero(d->var_inv);
                d->is_inverted = false;
            }
        }
        
        /* (re)allocating vector for chi^2 composition if necessary */
        if (d->chi2_comp == NULL)
        {
            d->chi2_comp = mat_create(Xsize+Ysize+2,1);
        }
        else if (nrow(d->chi2_comp) != Xsize+Ysize+2)
        {
            mat_destroy(d->chi2_comp);
            d->chi2_comp = mat_create(Xsize+Ysize+2,1);
        }
        
        /* inverting covariance matrices if necessary */
        if (!d->is_inverted)
        {
            tmp_covar = mat_create(ndata,ndata);
            
            /** inverting data variance matrix **/
            /*** (re)allocating y inverse variance matrix if necessary ***/
            if (d->y_var_inv == NULL)
            {
                d->y_var_inv = mat_create(Ysize,Ysize);
                mat_assume(d->y_var_inv,MAT_SYM|MAT_POS);
            }
            else if (nrow(d->y_var_inv) != Ysize)
            {
                mat_destroy(d->y_var_inv);
                d->y_var_inv = mat_create(Ysize,Ysize);
                mat_assume(d->y_var_inv,MAT_SYM|MAT_POS);
            }
            mat_zero(d->y_var_inv);
            /*** building y variance matrix by blocks ***/
            for (k1=0;k1<nydim;k1++)
            for (k2=k1;k2<nydim;k2++)
            {
                ind = sym_rowmaj(k1,k2,nydim);
                mat_mulp(tmp_covar,d->y_covar[ind],d->cor_filter);
                set_var_subm(d->y_var_inv,tmp_covar,k1,k2,d);
                if (k1 != k2)
                {
                    mat_eqtranspose(tmp_covar);
                    set_var_subm(d->y_var_inv,tmp_covar,k2,k1,d);
                }
            }
            if (have_xy_covar)
            {
                mat_set_subm(d->var_inv,d->y_var_inv,0,0,Ysize-1,Ysize-1);
            }
            /*** inversion ***/
            latan_printf(DEBUG1,"Cd=\n");
            if (latan_get_verb() == DEBUG1)
            {
                mat_print(d->y_var_inv,"% 10e");
            }
            pseudoinvert_var(d->y_var_inv,fit_data_is_y_correlated(d));
            latan_printf(DEBUG1,"Cd^-1=\n");
            if (latan_get_verb() == DEBUG1)
            {
                mat_print(d->y_var_inv,"% 10e");
            }
            /** inverting x variance matrix **/
            if (fit_data_have_x_var(d)||have_xy_covar)
            {
                /*** (re)allocating x inverse variance matrix if necessary ***/
                if (d->x_var_inv == NULL)
                {
                    d->x_var_inv = mat_create(Xsize,Xsize);
                    mat_assume(d->x_var_inv,MAT_SYM|MAT_POS);
                }
                else if (nrow(d->x_var_inv) != Xsize)
                {
                    mat_destroy(d->x_var_inv);
                    d->x_var_inv = mat_create(Xsize,Xsize);
                    mat_assume(d->x_var_inv,MAT_SYM|MAT_POS);
                }
                mat_zero(d->x_var_inv);
                /*** building x variance matrix by blocks ***/
                px_ind1 = 0;
                for (k1=0;k1<nxdim;k1++)
                {
                    if (fit_data_have_x_covar(d,k1))
                    {
                        px_ind2 = px_ind1;
                        for (k2=k1;k2<nxdim;k2++)
                        {
                            if (fit_data_have_x_covar(d,k2))
                            {
                                ind = sym_rowmaj(k1,k2,nxdim);
                                mat_mulp(tmp_covar,d->x_covar[ind],\
                                         d->cor_filter);
                                set_var_subm(d->x_var_inv,tmp_covar,px_ind1,\
                                             px_ind2,d);
                                if (k1 != k2)
                                {
                                    mat_eqtranspose(tmp_covar);
                                    set_var_subm(d->x_var_inv,tmp_covar,\
                                                 px_ind2,px_ind1,d);
                                }
                                px_ind2++;
                            }
                            
                        }
                        px_ind1++;
                    }
                }
                if (have_xy_covar)
                {
                    mat_set_subm(d->var_inv,d->x_var_inv,Ysize,Ysize,\
                                 lXsize-1,lXsize-1);
                }
                /*** inversion ***/
                latan_printf(DEBUG1,"Cx=\n");
                if (latan_get_verb() == DEBUG1)
                {
                    mat_print(d->x_var_inv,"% 10e");
                }
                pseudoinvert_var(d->x_var_inv,fit_data_is_x_correlated(d));
                latan_printf(DEBUG1,"Cx^-1=\n");
                if (latan_get_verb() == DEBUG1)
                {
                    mat_print(d->x_var_inv,"% 10e");
                }
                
            }
            /** inverting global variance matrix **/
            if (have_xy_covar)
            {
                for (k1=0;k1<nydim;k1++)
                {
                    px_ind = 0;
                    for (k2=0;k2<nxdim;k2++)
                    {
                        ind = rowmaj(k1,k2,nydim,nxdim);
                        if (d->have_xy_covar[k2])
                        {
                            mat_mulp(tmp_covar,d->xy_covar[ind],d->cor_filter);
                            set_var_subm(d->var_inv,tmp_covar,k1,px_ind+nydim,\
                                         d);
                            mat_eqtranspose(tmp_covar);
                            set_var_subm(d->var_inv,tmp_covar,px_ind+nydim,k1,\
                                         d);
                            px_ind++;
                        }
                    }
                }
                latan_printf(DEBUG1,"C=\n");
                if (latan_get_verb() == DEBUG1)
                {
                    mat_print(d->var_inv,"% 10e");
                }
                pseudoinvert_var(d->var_inv,true);
                latan_printf(DEBUG1,"C^-1=\n");
                if (latan_get_verb() == DEBUG1)
                {
                    mat_print(d->var_inv,"% 10e");
                }
            }
            d->is_inverted = true;
            
            mat_destroy(tmp_covar);
        }
    }
}

/* set X and Y */
static void set_X_Y(mat* X, mat *Y, mat *x_buf, const mat *p, const fit_data *d)
{
    size_t ndata,npt,nxdim,nydim,npar,px_ind;
    size_t i,k,k_i;
    double buf,eval;
    
    ndata   = fit_data_get_ndata(d);
    npt     = fit_data_fit_point_num(d);
    nxdim   = fit_data_get_nxdim(d);
    nydim   = fit_data_get_nydim(d);
    npar    = fit_data_get_npar(d);
    k_i     = 0;
    
    /* setting X and Y in case of covariance in x */
    if (fit_data_have_x_var(d))
    {
        for (i=0;i<ndata;i++)
        {
            px_ind = 0;
            if (fit_data_is_fit_point(d,i))
            {
                for (k=0;k<nxdim;k++)
                {
                    if (fit_data_have_x_covar(d,k))
                    {
                        buf = mat_get(p,npar+px_ind*npt+k_i,0);
                        mat_set(X,px_ind*npt+k_i,0,buf-fit_data_get_x(d,i,k));
                        px_ind++;
                    }
                    else
                    {
                        buf = fit_data_get_x(d,i,k);
                    }
                    mat_set(x_buf,k,0,buf);
                }
                for (k=0;k<nydim;k++)
                {
                    eval = fit_model_eval(d->model,k,x_buf,p,d->model_param);
                    mat_set(Y,k*ndata+k_i,0,eval-fit_data_get_y(d,i,k));
                }
                k_i++;
            }
        }
    }
    /* setting Y in case of no covariance in x */
    else
    {
        for (i=0;i<ndata;i++)
        {
            if (fit_data_is_fit_point(d,i))
            {
                for (k=0;k<nydim;k++)
                {
                    mat_set(Y,k*ndata+k_i,0,fit_data_model_eval(d,k,i,p)\
                            -fit_data_get_y(d,i,k));
                }
                k_i++;
            }
        }
    }
}

/* compute t(lX)*C^-1*lX */
static double chi2_base(const mat *p, void *vd)
{
    fit_data *d;
    int nthread,thread;
    mat *x_f,*Y,*CyY,*Cy,*X,*CxX,*Cx,*lX,*ClX,*C;
    double res,buf;
    
    d       = (fit_data *)vd;
#ifdef _OPENMP
    nthread = omp_get_num_threads();
    thread  = omp_get_thread_num();
#else
    nthread = 1;
    thread  = 0;
#endif

    /* buffers and inverse variance matrices initialization */
    init_chi2(d,thread,nthread);
    x_f = d->buf[thread].x_f;
    Y   = d->buf[thread].Y;
    CyY = d->buf[thread].CyY;
    Cy  = d->y_var_inv;
    X   = d->buf[thread].X;
    CxX = d->buf[thread].CxX;
    Cx  = d->x_var_inv;
    lX  = d->buf[thread].lX;
    ClX = d->buf[thread].ClX;
    C   = d->var_inv;

    /* setting X and Y */
    set_X_Y(X,Y,x_f,p,d);
    
    /* setting lX and computing chi^2 in case of data/x covariance */
    if (fit_data_have_xy_covar(d))
    {
        mat_set_subm(lX,Y,0,0,nrow(Y)-1,0);
        mat_set_subm(lX,X,nrow(Y),0,nrow(lX)-1,0);
        mat_mul(ClX,C,'n',lX,'n');
        latan_blas_ddot(ClX,lX,&res);
    }
    /* computing chi^2 by blocks in case of no data/x covariance */
    else
    {
        mat_mul(CyY,Cy,'n',Y,'n');
        latan_blas_ddot(CyY,Y,&res);
        if (fit_data_have_x_var(d))
        {
            mat_mul(CxX,Cx,'n',X,'n');
            latan_blas_ddot(CxX,X,&buf);
            res += buf;
        }
    }
    
    return res;
}

/* compute the total chi^2 (with an eventual custom extension) */
double chi2(const mat *p, void *vd)
{    
    return chi2_base(p,vd) + ((fit_data *)vd)->chi2_ext(p,vd);
}

/* compute chi^2 composition */
latan_errno chi2_get_comp(mat *comp, mat *p, fit_data *d)
{
    int nthread,thread;
    size_t Ysize,Xsize;
    size_t i;
    double base,uncor,el;
    mat *x,*Y,*X;
    
    Ysize = get_Ysize(d);
    Xsize = get_Xsize(d);
    
    if (nrow(comp) != Xsize + Ysize + 2)
    {
        LATAN_ERROR("wrong result matrix dimensions",LATAN_EBADLEN);
    }
    
    /* compute base chi^2 */
    base = chi2_base(p,d);
    
    /* compute ext chi^2 */
    mat_set(comp,Xsize+Ysize+1,0,d->chi2_ext(p,d));
    
    /* compute diagonal chi^2 elements */
#ifdef _OPENMP
    nthread = omp_get_num_threads();
    thread  = omp_get_thread_num();
#else
    nthread = 1;
    thread  = 0;
#endif
    /** buffers and inverse variance matrices initialization **/
    init_chi2(d,thread,nthread);
    x   = d->buf[thread].x_f;
    Y   = d->buf[thread].Y;
    X   = d->buf[thread].X;
    /** setting X and Y **/
    set_X_Y(X,Y,x,p,d);
    /** diagonal y elements **/
    uncor = 0.0;
    for (i=0;i<Ysize;i++)
    {
        el     = mat_get(Y,i,0)*sqrt(mat_get(d->y_var_inv,i,i));
        uncor += el;
        mat_set(comp,i,0,el);
    }
    /** diagonal x elements **/
    for (i=0;i<Xsize;i++)
    {
        el     = mat_get(X,i,0)*sqrt(mat_get(d->x_var_inv,i,i));
        uncor += el;
        mat_set(comp,i+Ysize,0,el);
    }
    
    /* compute correlation contribution */
    mat_set(comp,Xsize+Ysize,0,base-uncor);
    
    return LATAN_SUCCESS;
}

/*                          fit functions                                   */
/****************************************************************************/
latan_errno data_fit(mat *p, fit_data *d)
{
    latan_errno status;
    strbuf cor_status;
    double chi2_min;
    
    strbufcpy(cor_status,"correlations :");
    if (fit_data_is_y_correlated(d))
    {
        strcat(cor_status," data/data");
    }
    if (fit_data_is_x_correlated(d))
    {
        strcat(cor_status," x/x");
    }
    if (fit_data_have_xy_covar(d))
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
        d->chi2_val = chi2_min;
        chi2_get_comp(d->chi2_comp,p,d);
    }
    
    return status;
}

latan_errno rs_data_fit(rs_sample *p, rs_sample * const *x,         \
                          rs_sample * const *data, fit_data *d,       \
                          const cor_flag flag, const bool *use_x_var)
{
    latan_errno status;
    mat *pbuf,*comp_backup,*x_k;
    size_t npt,ndata,nxdim,nydim,npar,nsample,Xsize,Ysize,px_ind;
    size_t i,k,k_i,s;
    int verb_backup;
    double chi2_backup;
    
    verb_backup  = latan_get_verb();
    status       = LATAN_SUCCESS;
    npt          = fit_data_fit_point_num(d);
    ndata        = fit_data_get_ndata(d);
    nxdim        = fit_data_get_nxdim(d);
    nydim        = fit_data_get_nydim(d);
    npar         = fit_data_get_npar(d);
    nsample      = rs_sample_get_nsample(data[0]);
    
    /* compute needed variances/covariances from samples */
    fit_data_set_covar_from_sample(d,x,data,flag,use_x_var);    
    Xsize = get_Xsize(d);
    Ysize = get_Ysize(d);

    /* central value fit */
    d->s        = 0;
    pbuf        = mat_create(npar+Xsize,1);
    comp_backup = mat_create(Xsize+Ysize+2,1);
    /** setting data and initial parameters **/
    USTAT(mat_set_subm(pbuf,rs_sample_pt_cent_val(p),0,0,npar-1,0));
    for (k=0;k<nydim;k++)
    {
        USTAT(fit_data_set_y_k(d,k,rs_sample_pt_cent_val(data[k])));
    }
    if (x != NULL)
    {
        px_ind = 0;
        for (k=0;k<nxdim;k++)
        {
            USTAT(fit_data_set_x_k(d,k,rs_sample_pt_cent_val(x[k])));
            if (use_x_var[k])
            {
                k_i = 0;
                x_k = rs_sample_pt_cent_val(x[k]);
                for (i=0;i<ndata;i++)
                {
                    if (fit_data_is_fit_point(d,i))
                    {
                        mat_set(pbuf,npar+px_ind*npt+k_i,0,mat_get(x_k,i,0));
                        k_i++;
                    }
                }
                px_ind++;
            }
        }
    }
    /** looking at starting chi^2/dof **/
    latan_printf(VERB,"starting chi^2/dof = %f\n",         \
                 chi2(pbuf,d)/((double)fit_data_get_dof(d)));
    /**  fit **/
    USTAT(data_fit(pbuf,d));
    USTAT(mat_get_subm(rs_sample_pt_cent_val(p),pbuf,0,0,npar-1,0));
    chi2_backup = fit_data_get_chi2(d);
    fit_data_get_chi2_comp(comp_backup,d);
    latan_printf(VERB,"fit: central value chi^2/dof = %e -- p-value = %e\n",\
                 fit_data_get_chi2pdof(d),fit_data_get_pvalue(d));
    
    /* sample fits */
    for (s=0;s<nsample;s++)
    {
        (d->s)++;
        /** setting data and initial parameters **/
        USTAT(mat_set_subm(pbuf,rs_sample_pt_cent_val(p),0,0,npar-1,0));
        for (k=0;k<nydim;k++)
        {
            USTAT(fit_data_set_y_k(d,k,rs_sample_pt_sample(data[k],s)));
        }
        if (x != NULL)
        {
            px_ind = 0;
            for (k=0;k<nxdim;k++)
            {
                USTAT(fit_data_set_x_k(d,k,rs_sample_pt_sample(x[k],s)));
                if (use_x_var[k])
                {
                    k_i = 0;
                    x_k = rs_sample_pt_sample(x[k],s);
                    for (i=0;i<ndata;i++)
                    {
                        if (fit_data_is_fit_point(d,i))
                        {
                            mat_set(pbuf,npar+px_ind*npt+k_i,0,\
                                    mat_get(x_k,i,0));
                            k_i++;
                        }
                    }
                    px_ind++;
                }
            }
        }
        /** fit **/
        if (verb_backup != DEBUG2)
        {
            USTAT(latan_set_verb(QUIET));
        }
        USTAT(data_fit(pbuf,d));
        USTAT(latan_set_verb(verb_backup));
        if (latan_get_verb() == VERB)
        {
            printf("[");
            for (i=0;i<60*(s+1)/nsample;i++)
            {
                printf("=");
            }
            for (i=60*(s+1)/nsample;i<60;i++)
            {
                printf(" ");
            }
            printf("]  %d/%d\r",(int)s+1,(int)nsample);
            fflush(stdout);
        }
        latan_printf(DEBUG2,"fit: sample %d/%d chi^2/dof = %e\n",(int)s+1,\
                     (int)nsample,fit_data_get_chi2pdof(d));
        USTAT(mat_get_subm(rs_sample_pt_sample(p,s),pbuf,0,0,npar-1,0));
    }
    if (latan_get_verb() == VERB)
    {
        printf("\n");
    }
    d->chi2_val = chi2_backup;
    mat_cp(d->chi2_comp,comp_backup);
    for (k=0;k<nydim;k++)
    {
        USTAT(fit_data_set_y_k(d,k,rs_sample_pt_cent_val(data[k])));
    }
    if (x != NULL)
    {
        for (k=0;k<nxdim;k++)
        {
            USTAT(fit_data_set_x_k(d,k,rs_sample_pt_cent_val(x[k])));
        }
    }
    
    mat_destroy(pbuf);
    mat_destroy(comp_backup);
    
    return status;
}

void fit_residual(mat *res, const fit_data *d, const size_t ky, const mat *p)
{
    size_t i,ndata;
    double res_i;

    ndata = fit_data_get_ndata(d);
    
    for (i=0;i<ndata;i++)
    {
        res_i = fit_data_get_y(d,i,ky) - fit_data_model_eval(d,ky,i,p);
        mat_set(res,i,0,res_i);
    }
}

void fit_partresidual(mat *res, const fit_data *d, const size_t ky, \
                      const mat *x_ex, const size_t kx, const mat *p)
{
    size_t i,ndata,nxdim;
    double res_i;
    mat *x_buf;
    
    nxdim  = fit_data_get_nxdim(d);
    
    if (nxdim > 1)
    {
        ndata = fit_data_get_ndata(d);
        x_buf = mat_create(nxdim,1);
        for (i=0;i<ndata;i++)
        {
            mat_cp(x_buf,x_ex);
            mat_set(x_buf,kx,0,fit_data_get_x(d,i,kx));
            res_i = fit_data_get_y(d,i,ky) - fit_data_model_eval(d,ky,i,p)\
                    + fit_data_model_xeval(d,ky,x_buf,p);
            mat_set(res,i,0,res_i);
        }
        mat_destroy(x_buf);
    }
    else
    {
        fit_data_get_y_k(res,d,ky);
    }
}
