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
#include <latan/latan_math.h>
#include <latan/latan_minimizer.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <latan/latan_io.h>

/*                          fit model structure                             */
/****************************************************************************/
/** some useful constant npar_func **/
#define DEFINE_CST_NPAR_FUNC(n)\
size_t npar_##n(const unsigned int stage_flag, void *nothing)\
{\
    unsigned int dumb;\
    \
    dumb    = stage_flag;\
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

size_t fit_model_get_npar(const fit_model *model,                     \
                          const unsigned int stage_flag, void *model_param)
{
    size_t npar;

    npar = model->npar(stage_flag,model_param);

    return npar;
}

double fit_model_eval(const fit_model *model, mat *x, mat *p,    \
                      const unsigned int stage_flag, void *model_param)
{
    unsigned int i;
    double res;
    gsl_matrix_view p_view;
    mat subp;
    size_t ind_i,ind_f;
    
    res           = 0.0;
    ind_i         = 0;
    ind_f         = 0;
    
    for (i=0;i<MAX_STAGE;i++)
    {
        ind_f = ind_i + model->npar(STAGE(i),model_param) - 1;
        if (stage_flag & STAGE(i))
        {
            p_view = gsl_matrix_submatrix(p->data_cpu,ind_i,0,ind_f-ind_i+1,1);
            subp.data_cpu = &(p_view.matrix);
            res += model->func[i](x,&subp,model_param);   
        }
        ind_i = ind_f + 1;
    }
    
    return res;
}

/*                          fit data structure                              */
/****************************************************************************/
#define NCHI2BUF 9
enum
{
    i_X   = 0,
    i_CdX = 1,
    i_Cd  = 2,
    i_Y   = 3,
    i_CxY = 4,
    i_Cx  = 5,
    i_lX  = 6,
    i_lCX = 7,
    i_lC  = 8
};

/** allocation **/
fit_data *fit_data_create(const size_t ndata, const size_t ndim)
{
    fit_data *d;
    size_t i;
    
    MALLOC_ERRVAL(d,fit_data *,1,NULL);
    d->x           = mat_create(ndata*ndim,1);
    d->x_var       = mat_create(ndata*ndim,ndata*ndim);
    mat_id(d->x_var);
    mat_assume_sym(d->x_var,true);
    d->data        = mat_create(ndata,1);
    d->data_var    = mat_create(ndata,ndata);
    mat_id(d->data_var);
    mat_assume_sym(d->data_var,true);
    d->xdata_covar = mat_create(ndata,ndata*ndim);
    mat_zero(d->xdata_covar);
    d->buf_chi2[i_X]   = mat_create(ndata,1);
    d->buf_chi2[i_CdX] = mat_create(ndata,1);
    d->buf_chi2[i_Cd]  = mat_create(ndata,ndata);
    mat_assume_sym(d->buf_chi2[i_Cd],true);
    d->buf_chi2[i_Y]   = mat_create(ndata*ndim,1);
    d->buf_chi2[i_CxY] = mat_create(ndata*ndim,1);
    d->buf_chi2[i_Cx]  = mat_create(ndata*ndim,ndata*ndim);
    mat_assume_sym(d->buf_chi2[i_Cx],true);
    d->buf_chi2[i_lX]  = mat_create(ndata*(ndim+1),1);
    d->buf_chi2[i_lCX] = mat_create(ndata*(ndim+1),1);
    d->buf_chi2[i_lC]  = mat_create(ndata*(ndim+1),ndata*(ndim+1));
    mat_assume_sym(d->buf_chi2[i_lC],true);
    d->is_inverted     = false;
    MALLOC_ERRVAL(d->to_fit,bool*,ndata,NULL);
    for (i=0;i<ndata;i++)
    {
        d->to_fit[i] = false;
    }
    d->model              = NULL;
    d->model_param        = NULL;
    d->stage_flag         = 1;
    d->ndata              = ndata;
    d->ndim               = ndim;
    d->is_data_correlated = false;
    d->is_x_correlated    = false;
    d->have_x_var         = false;
    d->have_xdata_covar   = false;
    d->save_chi2pdof      = true;
    
    return d;
}

void fit_data_destroy(fit_data *d)
{
    int i;

    mat_destroy(d->x);
    mat_destroy(d->x_var);
    mat_destroy(d->data);
    mat_destroy(d->data_var);
    for(i=0;i<NCHI2BUF;i++)
    {
        mat_destroy(d->buf_chi2[i]);
    }
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
void fit_data_set_x(fit_data *d, const size_t i, const size_t j,\
                    const double x_ij)
{
    mat_set(d->x,i*d->ndim+j,0,x_ij);
}

double fit_data_get_x(fit_data *d, const size_t i, const size_t j)
{
    return mat_get(d->x,i*d->ndim+j,0);
}

mat *fit_data_pt_x(const fit_data *d)
{
    return d->x;
}

latan_errno fit_data_set_x_var(fit_data *d, mat *var)
{
    latan_errno status;
    
    d->have_x_var      = true;
    d->is_x_correlated = mat_is_square(var);
    d->is_inverted     = false;
    
    if (d->is_x_correlated)
    {
        status = mat_cp(d->x_var,var);
    }
    else
    {
        status = mat_set_diag(d->x_var,var);
    }

    return status;
}

bool fit_data_have_x_var(const fit_data *d)
{
    return d->have_x_var;
}

bool fit_data_is_x_correlated(const fit_data *d)
{
    return d->is_x_correlated;
}

void fit_data_fit_all_points(fit_data *d, bool fit)
{
    size_t i;
    
    for (i=0;i<d->ndata;i++)
    {
        d->to_fit[i] = fit;
    }
}

void fit_data_fit_point(fit_data *d, size_t i, bool fit)
{
    if (i>=d->ndata)
    {
        LATAN_ERROR_VOID("index out of range",LATAN_EBADLEN);
    }
    
    d->to_fit[i] = fit;
}

void fit_data_fit_range(fit_data *d, size_t start, size_t end, bool fit)
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

bool fit_data_is_fit_point(const fit_data *d, size_t i)
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

double fit_data_get_data(fit_data *d, const size_t i)
{
    return mat_get(d->data,i,0);
}

mat *fit_data_pt_data(const fit_data *d)
{
    return d->data;
}

latan_errno fit_data_set_data_var(fit_data *d, mat *var)
{
    latan_errno status;
    
    d->is_data_correlated = mat_is_square(var);
    d->is_inverted        = false;
    
    if (d->is_data_correlated)
    {
        status = mat_cp(d->data_var,var);
    }
    else
    {
        status = mat_set_diag(d->data_var,var);
    }
    
    return status;
}

bool fit_data_is_data_correlated(const fit_data *d)
{
    return d->is_data_correlated;
}

latan_errno fit_data_set_xdata_covar(fit_data *d, mat *covar)
{
    latan_errno status;

    status              = mat_cp(d->xdata_covar,covar);
    d->have_xdata_covar = true;
    d->is_inverted      = false;

    return status;
}

bool fit_data_have_xdata_covar(const fit_data *d)
{
    return d->have_xdata_covar;
}

/*** model ***/
latan_errno fit_data_set_model(fit_data *d, const fit_model *model)
{
    if (model->ndim != d->ndim)
    {
        LATAN_ERROR("fit model and fit data dimension mismatch",LATAN_EBADLEN);
    }
    
    d->model = model;
    
    return LATAN_SUCCESS;
}

const fit_model *fit_data_pt_model(fit_data *d)
{
    return d->model;
}

void fit_data_set_model_param(fit_data *d, void *model_param)
{
    d->model_param = model_param;
}

/* CRITICALLY CALLED FUNCTION
 * fit_data_model_eval is heavily called during one call of minimize function,
 * optimization is done using MAT_PT_SUBM macro
 */
double fit_data_model_eval(const fit_data *d, const size_t i,\
                           mat *p)
{
    mat x_i;
    gsl_matrix_view x_view;
    double res;
    
    x_view       = gsl_matrix_submatrix(d->x->data_cpu,i*d->ndim,0,d->ndim,1);
    x_i.data_cpu = &(x_view.matrix);
    res          = fit_model_eval(d->model,&x_i,p,d->stage_flag,d->model_param);
    
    return res;
}

/*** stages ***/
void fit_data_set_stage_flag(fit_data *d, const unsigned int stage_flag)
{
    d->stage_flag = stage_flag;
}

unsigned int fit_data_get_stage_flag(const fit_data *d)
{
    return d->stage_flag;
}

size_t fit_data_get_npar(const fit_data *d)
{
    size_t npar;

    npar = fit_model_get_npar(d->model,d->stage_flag,d->model_param);

    return npar;
}

void fit_data_set_stages(fit_data *d, const stage_ar s)
{
    unsigned int i;
    
    d->stage_flag = 0;
    for (i=0;i<MAX_STAGE;i++)
    {
        if (s[i])
        {
            d->stage_flag |= STAGE(i);
        }
    }
}

void fit_data_get_stages(stage_ar s, const fit_data *d)
{
    unsigned int i;
    
    for (i=0;i<MAX_STAGE;i++)
    {
        s[i] = ((d->stage_flag & STAGE(i)) != 0);
    }
}

/*** dof ***/
int fit_data_get_dof(const fit_data *d)
{
    return fit_data_fit_point_num(d) - d->model->npar(d->stage_flag,\
                                                      d->model_param);
}

/*                          chi2 function                                   */
/****************************************************************************/
static void invert_var(fit_data *d)
{
    mat *lC,*txdata_covar;
    size_t i;

    lC           = d->buf_chi2[i_lC];
    txdata_covar = NULL;

    if (d->have_xdata_covar)
    {
        txdata_covar = mat_create_from_trdim(d->xdata_covar);
        
        mat_set_subm(lC,d->data_var,0,0,d->ndata-1,d->ndata-1);
        mat_set_subm(lC,d->x_var,d->ndata,d->ndata,nrow(lC)-1,ncol(lC)-1);
        mat_set_subm(lC,d->xdata_covar,0,d->ndata,d->ndata-1,ncol(lC)-1);
        mat_transpose(txdata_covar,d->xdata_covar);
        mat_set_subm(lC,txdata_covar,d->ndata,0,nrow(lC)-1,d->ndata-1);
        mat_eqinv(lC);

        mat_destroy(txdata_covar);
    }
    else
    {
        if (d->have_x_var)
        {
            if (d->is_x_correlated)
            {
                mat_inv(d->buf_chi2[i_Cx],d->x_var);
            }
            else
            {
                mat_zero(d->buf_chi2[i_Cx]);
                for (i=0;i<d->ndata*d->ndim;i++)
                {
                    mat_set(d->buf_chi2[i_Cx],i,i,1.0/mat_get(d->x_var,i,i));
                }
            }
        }
        if (d->is_data_correlated)
        {
            mat_inv(d->buf_chi2[i_Cd],d->data_var);
        }
        else
        {
            mat_zero(d->buf_chi2[i_Cd]);
            for (i=0;i<d->ndata;i++)
            {
                mat_set(d->buf_chi2[i_Cd],i,i,1.0/mat_get(d->data_var,i,i));
            }
        }
    }
    d->is_inverted = true;
}

/* CRITICALLY CALLED FUNCTION
 * chi2 is heavily called during one call of minimize function,
 * optimization is done using vector/matrix views from GSL, matrix buffers
 * in fit_data structure, and ddot BLAS operation
 */
double chi2(mat *p, void *d)
{
    fit_data *dt;
    size_t ndata,npar;
    size_t i,j;
    mat *X,*CdX,*Cd,*Y,*CxY,*Cx,*lX,*lCX,*lC,x_i;
    gsl_matrix_view x_view;
    double res,buf,eval;
    
    dt    = (fit_data *)d;
    ndata = fit_data_get_ndata(dt);
    npar  = fit_data_get_npar(d);
    
#ifdef _OPENMP
    #pragma omp critical
#endif
    {
        X   = dt->buf_chi2[i_X];
        CdX = dt->buf_chi2[i_CdX];
        Cd  = dt->buf_chi2[i_Cd];
        Y   = dt->buf_chi2[i_Y];
        CxY = dt->buf_chi2[i_CxY];
        Cx  = dt->buf_chi2[i_Cx];
        lX  = dt->buf_chi2[i_lX];
        lCX = dt->buf_chi2[i_lCX];
        lC  = dt->buf_chi2[i_lC];

        if (!dt->is_inverted)
        {
            latan_printf(VERB,"inverting variance matrix...\n");
            invert_var(dt);
        }
        if (dt->have_x_var)
        {
            for (i=0;i<ndata;i++)
            {
                if (fit_data_is_fit_point(dt,i))
                {
                    x_view       = gsl_matrix_submatrix(p->data_cpu,    \
                                                        i*dt->ndim+npar,\
                                                        0,dt->ndim,1);
                    x_i.data_cpu = &(x_view.matrix);
                    eval         = fit_model_eval(dt->model,&x_i,p,     \
                                                  dt->stage_flag,       \
                                                  dt->model_param);
                    mat_set(X,i,0,eval - fit_data_get_data(dt,i));
                    for (j=0;j<dt->ndim;j++)
                    {
                        mat_set(Y,i*dt->ndim+j,0,                       \
                                mat_get(&x_i,j,0)-fit_data_get_x(dt,i,j));
                    }
                }
                else
                {
                    for (j=0;j<dt->ndim;j++)
                    {
                        mat_set(Y,i*dt->ndim+j,0,0.0);
                    }
                }
            }
        }
        else
        {
            for (i=0;i<ndata;i++)
            {
                if (fit_data_is_fit_point(dt,i))
                {
                    mat_set(X,i,0,fit_data_model_eval(dt,i,p)
                            - fit_data_get_data(dt,i));
                }
                else
                {
                    mat_set(X,i,0,0.0);
                }
            }
        }
        if (dt->have_xdata_covar)
        {
            mat_set_subm(lX,X,0,0,dt->ndata-1,0);
            mat_set_subm(lX,Y,dt->ndata,0,nrow(lX)-1,0);
            mat_mul(lCX,lC,'n',lX,'n');
            latan_blas_ddot(lCX,lX,&res);
        }
        else
        {
            mat_mul(CdX,Cd,'n',X,'n');
            latan_blas_ddot(CdX,X,&res);
            if (dt->have_x_var)
            {
                mat_mul(CxY,Cx,'n',Y,'n');
                latan_blas_ddot(CxY,Y,&buf);
                res += buf;
            }
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
    
    fit_data_is_data_correlated(d) ? \
    (strbufcpy(cor_status,"correlated")) : (strbufcpy(cor_status,"uncorrelated"));
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

latan_errno rs_data_fit(rs_sample *p, rs_sample *data, fit_data *d)
{
    latan_errno status;
    mat *datavar;
    size_t i;
    int verb_backup;
    double chi2pdof_backup;
        
    verb_backup = latan_get_verb();
    status = LATAN_SUCCESS;
    
    datavar = mat_create(nrow(rs_sample_pt_cent_val(data)),nrow(rs_sample_pt_cent_val(data)));
    
    LATAN_UPDATE_STATUS(status,rs_sample_varp(datavar,data));
    LATAN_UPDATE_STATUS(status,fit_data_set_data_var(d,datavar));
    mat_cp(fit_data_pt_data(d),rs_sample_pt_cent_val(data));
    LATAN_UPDATE_STATUS(status,data_fit(rs_sample_pt_cent_val(p),d));
    chi2pdof_backup = fit_data_get_chi2pdof(d);
    latan_printf(DEBUG,"central value chi^2/dof = %e\n",chi2pdof_backup);
    if (verb_backup != DEBUG)
    {
        latan_set_verb(QUIET);
    }
    for (i=0;i<rs_sample_get_nsample(data);i++)
    {
        mat_cp(rs_sample_pt_sample(p,i),\
               rs_sample_pt_cent_val(p));
        mat_cp(fit_data_pt_data(d),rs_sample_pt_sample(data,i));
        LATAN_UPDATE_STATUS(status,data_fit(rs_sample_pt_sample(p,i),d));
        latan_printf(DEBUG,"sample %lu chi^2/dof = %e\n",(long unsigned)i,\
                     fit_data_get_chi2pdof(d));
    }
    d->chi2pdof = chi2pdof_backup;
    latan_set_verb(verb_backup);
    
    mat_destroy(datavar);
    
    return status;
}

latan_errno rs_x_data_fit(rs_sample *p, rs_sample *x, rs_sample *data,\
                          fit_data *d)
{
    latan_errno status;
    mat *datavar;
    size_t i;
    int verb_backup;
    double chi2pdof_backup;
    
    verb_backup = latan_get_verb();
    status = LATAN_SUCCESS;
    
    datavar = mat_create(nrow(rs_sample_pt_cent_val(data)),1);
    
    LATAN_UPDATE_STATUS(status,rs_sample_varp(datavar,data));
    LATAN_UPDATE_STATUS(status,fit_data_set_data_var(d,datavar));
    mat_cp(fit_data_pt_data(d),rs_sample_pt_cent_val(data));
    mat_cp(fit_data_pt_x(d),rs_sample_pt_cent_val(x));
    LATAN_UPDATE_STATUS(status,data_fit(rs_sample_pt_cent_val(p),d));
    chi2pdof_backup = fit_data_get_chi2pdof(d);
    latan_printf(DEBUG,"central value chi^2/dof = %e\n",chi2pdof_backup);
    if (verb_backup != DEBUG)
    {
        latan_set_verb(QUIET);
    }
    for (i=0;i<rs_sample_get_nsample(data);i++)
    {
        mat_cp(rs_sample_pt_sample(p,i),\
               rs_sample_pt_cent_val(p));
        mat_cp(fit_data_pt_data(d),rs_sample_pt_sample(data,i));
        mat_cp(fit_data_pt_x(d),rs_sample_pt_sample(x,i));
        LATAN_UPDATE_STATUS(status,data_fit(rs_sample_pt_sample(p,i),d));
        latan_printf(DEBUG,"sample %lu chi^2/dof = %e\n",(long unsigned)i,\
                     fit_data_get_chi2pdof(d));
    }
    d->chi2pdof = chi2pdof_backup;
    latan_set_verb(verb_backup);
    
    mat_destroy(datavar);
    
    return status;
}

void fit_residual(mat *res, mat *p, fit_data *d)
{
    size_t i,ndata;
    double res_i;

    ndata = fit_data_get_ndata(d);
    
    for (i=0;i<ndata;i++)
    {
        res_i = fit_data_get_data(d,i)/fit_data_model_eval(d,i,p) - 1.0;
        mat_set(res,i,0,res_i);
    }
}

void rs_fit_residual(rs_sample *res, rs_sample *p, rs_sample *data,\
                     fit_data *d)
{
    size_t i,nsample;
    
    nsample = rs_sample_get_nsample(data);

    mat_cp(fit_data_pt_data(d),rs_sample_pt_cent_val(data));
    fit_residual(rs_sample_pt_cent_val(res),rs_sample_pt_cent_val(p),d);
    for (i=0;i<nsample;i++)
    {
        mat_cp(fit_data_pt_data(d),rs_sample_pt_sample(data,i));
        fit_residual(rs_sample_pt_sample(res,i),rs_sample_pt_sample(p,i),d);
    }
}

void rs_x_fit_residual(rs_sample *res, rs_sample *p, rs_sample *x,\
                       rs_sample *data, fit_data *d)
{
    size_t i,nsample;

    nsample = rs_sample_get_nsample(data);

    mat_cp(fit_data_pt_x(d),rs_sample_pt_cent_val(x));
    mat_cp(fit_data_pt_data(d),rs_sample_pt_cent_val(data));
    fit_residual(rs_sample_pt_cent_val(res),rs_sample_pt_cent_val(p),d);
    for (i=0;i<nsample;i++)
    {
        mat_cp(fit_data_pt_x(d),rs_sample_pt_sample(x,i));
        mat_cp(fit_data_pt_data(d),rs_sample_pt_sample(data,i));
        fit_residual(rs_sample_pt_sample(res,i),rs_sample_pt_sample(p,i),d);
    }
}
