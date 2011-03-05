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

size_t fit_model_get_npar(const fit_model *model,                     \
                          const unsigned int stage_flag, void *model_param)
{
    size_t npar;
    size_t i;

    npar = 0;

    for (i=0;i<model->nstage;i++)
    {
        if (stage_flag & STAGE(i))
        {
            npar += model->npar[i](model_param);
        }
    }

    return npar;
}

static double fit_model_eval_ker(const fit_model *model, const mat *x, \
                                 const mat *p, const stage_ar do_stage,\
                                 const size_t stage_npar[MAX_STAGE],   \
                                 void *model_param)
{
    size_t i;
    double res;
    size_t ind_i,ind_f;
    gsl_matrix_view p_view;
    mat subp;

    res   = 0.0;
    ind_i = 0;
    ind_f = 0;
    
    for (i=0;i<model->nstage;i++)
    {
        if (do_stage[i])
        {
            ind_f  = ind_i + stage_npar[i] - 1;
            p_view = gsl_matrix_submatrix(p->data_cpu,ind_i,0,ind_f-ind_i+1,1);
            subp.data_cpu = &(p_view.matrix);
            res   += model->func[i](x,&subp,model_param);
            ind_i  = ind_f + 1;
        }
    }

    return res;
}

static void init_npar(stage_ar do_stage, size_t stage_npar[MAX_STAGE],      \
                      const fit_model *model, const unsigned int stage_flag,\
                      void *model_param)
{
    size_t i;

    for (i=0;i<model->nstage;i++)
    {
        do_stage[i]   = ((stage_flag & STAGE(i)) != 0);
        stage_npar[i] = model->npar[i](model_param);
    }
}

double fit_model_eval(const fit_model *model, mat *x, mat *p,    \
                      const unsigned int stage_flag, void *model_param)
{
    double res;
    stage_ar do_stage;
    size_t stage_npar[MAX_STAGE];

    init_npar(do_stage,stage_npar,model,stage_flag,model_param);
    res = fit_model_eval_ker(model,x,p,do_stage,stage_npar,model_param);

    return res;
}

/*                          fit data structure                              */
/****************************************************************************/
/** allocation **/
fit_data *fit_data_create(const size_t ndata, const size_t ndim)
{
    fit_data *d;
    size_t i;
    
    MALLOC_ERRVAL(d,fit_data *,1,NULL);
    MALLOC_ERRVAL(d->to_fit,bool*,ndata,NULL);

    d->ndata              = ndata;
    d->ndim               = ndim;
    d->npar               = 0;
    d->x                  = mat_create(ndata*ndim,1);
    d->x_var              = mat_create(ndata*ndim,ndata*ndim);
    d->x_var_inv          = mat_create(ndata*ndim,ndata*ndim);
    d->is_x_correlated    = false;
    d->have_x_var         = false;
    for (i=0;i<ndata;i++)
    {
        d->to_fit[i] = false;
    }
    d->data               = mat_create(ndata,1);
    d->data_var           = mat_create(ndata,ndata);
    d->data_var_inv       = mat_create(ndata,ndata);
    d->is_data_correlated = false;
    d->xdata_covar        = mat_create(ndata,ndata*ndim);
    d->have_xdata_covar   = false;
    d->var_inv            = mat_create(ndata*(ndim+1),ndata*(ndim+1));
    d->is_init            = false;
    d->model              = NULL;
    d->model_param        = NULL;
    d->stage_flag         = 1;
    for (i=0;i<MAX_STAGE;i++)
    {
        d->do_stage[i]   = false;
        d->stage_npar[i] = 0;
    }
    d->chi2pdof           = -1.0;
    d->save_chi2pdof      = true;
    d->buf                = NULL;
    d->nbuf               = 0;

    mat_id(d->x_var);
    mat_assume_sym(d->x_var,true);
    mat_id(d->data_var);
    mat_assume_sym(d->data_var,true);
    mat_zero(d->xdata_covar);

    return d;
}

void fit_data_destroy(fit_data *d)
{
    int i;

    mat_destroy(d->x);
    mat_destroy(d->x_var);
    mat_destroy(d->x_var_inv);
    mat_destroy(d->data);
    mat_destroy(d->data_var);
    mat_destroy(d->data_var_inv);
    mat_destroy(d->xdata_covar);
    mat_destroy(d->var_inv);
    for(i=0;i<d->nbuf;i++)
    {
        mat_destroy(d->buf[i].X);
        mat_destroy(d->buf[i].CdX);
        mat_destroy(d->buf[i].Y);
        mat_destroy(d->buf[i].CxY);
        mat_destroy(d->buf[i].lX);
        mat_destroy(d->buf[i].ClX);
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
    d->is_init     = false;
    
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
    d->is_init        = false;
    
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
    d->is_init          = false;

    return status;
}

bool fit_data_have_xdata_covar(const fit_data *d)
{
    return d->have_xdata_covar;
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
    d->npar        = fit_model_get_npar(d->model,d->stage_flag,d->model_param);
    init_npar(d->do_stage,d->stage_npar,d->model,d->stage_flag,d->model_param);

    return LATAN_SUCCESS;
}

/* CRITICALLY CALLED FUNCTION
 * fit_data_model_eval is heavily called during one call of minimize function,
 * optimization is done using GSL matrix view
 */
double fit_data_model_eval(const fit_data *d, const size_t i,\
                           mat *p)
{
    mat x_i;
    gsl_matrix_view x_view;
    double res;
    
    x_view       = gsl_matrix_submatrix(d->x->data_cpu,i*d->ndim,0,d->ndim,1);
    x_i.data_cpu = &(x_view.matrix);
    res          = fit_model_eval_ker(d->model,&x_i,p,d->do_stage,\
                                      d->stage_npar,d->model_param);
    
    return res;
}

/*** stages ***/
void fit_data_set_stage_flag(fit_data *d, const unsigned int stage_flag)
{
    d->stage_flag = stage_flag;
    if (d->model != NULL)
    {
        d->npar = fit_model_get_npar(d->model,d->stage_flag,d->model_param);
        init_npar(d->do_stage,d->stage_npar,d->model,d->stage_flag,\
                  d->model_param);
    }
}

unsigned int fit_data_get_stage_flag(const fit_data *d)
{
    return d->stage_flag;
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
    if (d->model != NULL)
    {
        d->npar = fit_model_get_npar(d->model,d->stage_flag,d->model_param);
        init_npar(d->do_stage,d->stage_npar,d->model,d->stage_flag,\
                  d->model_param);
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
    return fit_data_fit_point_num(d)                                 \
           - fit_model_get_npar(d->model,d->stage_flag,d->model_param);
}

/*                          chi2 function                                   */
/****************************************************************************/
static void init_chi2(fit_data *d, int nthread)
{
    mat *txdata_covar;
    int t;
    size_t i;
    size_t ndata,ndim,ldim;
    double inv;

    txdata_covar = NULL;
    ndata        = d->ndata;
    ndim         = d->ndim;
    ldim         = ndata*(ndim+1);

#ifdef _OPENMP
    #pragma omp critical
#endif
    {
        if (nthread > d->nbuf)
        {
            REALLOC_NOERRET(d->buf,d->buf,chi2_buf *,nthread);
            for (t=d->nbuf;t<nthread;t++)
            {
                d->buf[t].X   = mat_create(ndata,1);
                d->buf[t].CdX = mat_create(ndata,1);
                d->buf[t].Y   = mat_create(ndata*ndim,1);
                d->buf[t].CxY = mat_create(ndata*ndim,1);
                d->buf[t].lX  = mat_create(ldim,1);
                d->buf[t].ClX = mat_create(ldim,1);
            }
            d->nbuf = nthread;
        }
        if (!d->is_init)
        {
            if (d->have_xdata_covar)
            {
                txdata_covar = mat_create_from_trdim(d->xdata_covar);

                mat_set_subm(d->var_inv,d->data_var,0,0,ndata-1,ndata-1);
                mat_set_subm(d->var_inv,d->x_var,ndata,ndata,ldim-1,ldim-1);
                mat_set_subm(d->var_inv,d->xdata_covar,0,ndata,ndata-1,ldim-1);
                mat_transpose(txdata_covar,d->xdata_covar);
                mat_set_subm(d->var_inv,txdata_covar,ndata,0,ldim-1,ndata-1);
                mat_eqinv(d->var_inv);
                latan_printf(DEBUG,"C^-1 = \n");
                if (latan_get_verb() == DEBUG)
                {
                    mat_print(d->var_inv,"%6.1e");
                }

                mat_destroy(txdata_covar);
            }
            else
            {
                if (d->have_x_var)
                {
                    if (d->is_x_correlated)
                    {
                        mat_inv(d->x_var_inv,d->x_var);
                    }
                    else
                    {
                        mat_zero(d->x_var_inv);
                        for (i=0;i<d->ndata*d->ndim;i++)
                        {
                            inv = 1.0/mat_get(d->x_var,i,i);
                            if (gsl_isinf(inv))
                            {
                                mat_set(d->x_var_inv,i,i,0.0);
                            }
                            else
                            {
                                mat_set(d->x_var_inv,i,i,inv);
                            }
                        }
                    }
                    latan_printf(DEBUG,"Cx^-1 = \n");
                    if (latan_get_verb() == DEBUG)
                    {
                        mat_print(d->x_var_inv,"%6.1e");
                    }
                }
                if (d->is_data_correlated)
                {
                    mat_inv(d->data_var_inv,d->data_var);
                }
                else
                {
                    mat_zero(d->data_var_inv);
                    for (i=0;i<d->ndata;i++)
                    {
                        inv = 1.0/mat_get(d->data_var,i,i);
                        if (gsl_isinf(inv))
                        {
                            mat_set(d->data_var_inv,i,i,0.0);
                        }
                        else
                        {
                            mat_set(d->data_var_inv,i,i,inv);
                        }
                    }
                }
                latan_printf(DEBUG,"Cd^-1 = \n");
                if (latan_get_verb() >= DEBUG)
                {
                    mat_print(d->data_var_inv,"%6.1e");
                }
            }
            d->is_init = true;
        }
    }
}

/* CRITICALLY CALLED FUNCTION
 * chi2 is heavily called during one call of minimize function,
 * optimization is done using vector/matrix views from GSL, matrix buffers
 * in fit_data structure, and ddot BLAS operation
 */
double chi2(mat *p, void *vd)
{
    fit_data *d;
    int nthread,thread;
    size_t ndata,npar;
    size_t i,j;
    mat *X,*CdX,*Cd,*Y,*CxY,*Cx,*lX,*ClX,*C,x_i;
    gsl_matrix_view x_view;
    double res,buf,eval;
    
    d       = (fit_data *)vd;
    ndata   = fit_data_get_ndata(d);
    npar    = fit_data_get_npar(d);
#ifdef _OPENMP
    nthread = omp_get_num_threads();
    thread  = omp_get_thread_num();
#else
    nthread = 1;
    thread  = 0;
#endif

    init_chi2(d,nthread);
    X      = d->buf[thread].X;
    CdX    = d->buf[thread].CdX;
    Cd     = d->data_var_inv;
    Y      = d->buf[thread].Y;
    CxY    = d->buf[thread].CxY;
    Cx     = d->x_var_inv;
    lX     = d->buf[thread].lX;
    ClX    = d->buf[thread].ClX;
    C      = d->var_inv;
    if (d->have_x_var)
    {
        for (i=0;i<ndata;i++)
        {
            if (fit_data_is_fit_point(d,i))
            {
                x_view       = gsl_matrix_submatrix(p->data_cpu,    \
                                                    i*d->ndim+npar,\
                                                    0,d->ndim,1);
                x_i.data_cpu = &(x_view.matrix);
                eval         = fit_model_eval(d->model,&x_i,p,     \
                                              d->stage_flag,       \
                                              d->model_param);
                mat_set(X,i,0,eval - fit_data_get_data(d,i));
                for (j=0;j<d->ndim;j++)
                {
                    mat_set(Y,i*d->ndim+j,0,                       \
                            mat_get(&x_i,j,0)-fit_data_get_x(d,i,j));
                }
            }
            else
            {
                for (j=0;j<d->ndim;j++)
                {
                    mat_set(Y,i*d->ndim+j,0,0.0);
                }
            }
        }
    }
    else
    {
        for (i=0;i<ndata;i++)
        {
            if (fit_data_is_fit_point(d,i))
            {
                mat_set(X,i,0,fit_data_model_eval(d,i,p)
                        - fit_data_get_data(d,i));
            }
            else
            {
                mat_set(X,i,0,0.0);
            }
        }
    }
    if (d->have_xdata_covar)
    {
        mat_set_subm(lX,X,0,0,d->ndata-1,0);
        mat_set_subm(lX,Y,d->ndata,0,nrow(lX)-1,0);
        mat_mul(ClX,C,'n',lX,'n');
        latan_blas_ddot(ClX,lX,&res);
    }
    else
    {
        mat_mul(CdX,Cd,'n',X,'n');
        latan_blas_ddot(CdX,X,&res);
        if (d->have_x_var)
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
        LATAN_UPDATE_STATUS(status,rs_sample_var(datavar,data));
    }
    else
    {
        LATAN_UPDATE_STATUS(status,rs_sample_varp(datavar,data));
    }
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

latan_errno rs_x_data_fit(rs_sample *p, const rs_sample *x,  \
                          const rs_sample *data, fit_data *d,\
                          const cor_flag flag)
{
    latan_errno status;
    mat *datavar,*xvar,*xdatacovar,*pbuf;
    size_t ndata,ndim,npar;
    size_t i,j,k;
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
    xvar         = mat_create(ndata*ndim,is_x_cor ? ndata*ndim : 1);
    xdatacovar   = mat_create(ndata,ndata*ndim);
    pbuf         = mat_create(npar + ndata*ndim,1);

    if (is_data_cor)
    {
        LATAN_UPDATE_STATUS(status,rs_sample_var(datavar,data));
    }
    else
    {
        LATAN_UPDATE_STATUS(status,rs_sample_varp(datavar,data));
    }
    LATAN_UPDATE_STATUS(status,fit_data_set_data_var(d,datavar));
    if (is_x_cor)
    {
        LATAN_UPDATE_STATUS(status,rs_sample_var(xvar,x));
        if (!is_data_cor)
        {
            for (i=0;i<ndata;i++)
            {
                for (j=0;j<ndim;j++)
                {
                    for (k=0;k<(ndata-1)*ndim;k++)
                    {
                        mat_set(xvar,((i+1)*ndim+k)%(ndim*ndata),i*ndim+j,0.0);
                    }
                }
            }
        }
    }
    else
    {
        LATAN_UPDATE_STATUS(status,rs_sample_varp(xvar,x));
    }
    LATAN_UPDATE_STATUS(status,fit_data_set_x_var(d,xvar));
    if (is_xdata_cor)
    {
        LATAN_UPDATE_STATUS(status,rs_sample_cov(xdatacovar,data,x));
        if (!is_data_cor)
        {
            for (i=0;i<ndata;i++)
            {
                for (k=0;k<(ndata-1)*ndim;k++)
                {
                    mat_set(xdatacovar,i,((i+1)*ndim+k)%(ndim*ndata),0.0);
                }
            }
        }
        LATAN_UPDATE_STATUS(status,fit_data_set_xdata_covar(d,xdatacovar));
    }
    mat_cp(fit_data_pt_data(d),rs_sample_pt_cent_val(data));
    mat_cp(fit_data_pt_x(d),rs_sample_pt_cent_val(x));
    mat_set_subm(pbuf,rs_sample_pt_cent_val(p),0,0,npar-1,0);
    mat_set_subm(pbuf,fit_data_pt_x(d),npar,0,nrow(pbuf)-1,0);
    latan_printf(VERB,"starting chi^2/dof = %f\n",         \
                 chi2(pbuf,d)/((double)fit_data_get_dof(d)));
    LATAN_UPDATE_STATUS(status,data_fit(pbuf,d));
    mat_get_subm(rs_sample_pt_cent_val(p),pbuf,0,0,npar-1,0);
    chi2pdof_backup = fit_data_get_chi2pdof(d);
    latan_printf(DEBUG,"central value chi^2/dof = %e\n",chi2pdof_backup);
    if (verb_backup != DEBUG)
    {
        latan_set_verb(QUIET);
    }
    for (i=0;i<rs_sample_get_nsample(data);i++)
    {
        mat_cp(fit_data_pt_data(d),rs_sample_pt_sample(data,i));
        mat_cp(fit_data_pt_x(d),rs_sample_pt_sample(x,i));
        mat_set_subm(pbuf,fit_data_pt_x(d),npar,0,nrow(pbuf)-1,0);
        LATAN_UPDATE_STATUS(status,data_fit(pbuf,d));
        mat_get_subm(rs_sample_pt_sample(p,i),pbuf,0,0,npar-1,0);
        latan_printf(DEBUG,"sample %lu chi^2/dof = %e\n",(long unsigned)i,\
                     fit_data_get_chi2pdof(d));
    }
    d->chi2pdof = chi2pdof_backup;
    latan_set_verb(verb_backup);
    
    mat_destroy(datavar);
    mat_destroy(xvar);
    mat_destroy(xdatacovar);
    mat_destroy(pbuf);
    
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
