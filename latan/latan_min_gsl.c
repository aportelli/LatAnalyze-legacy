/* latan_min_gsl.c, part of LatAnalyze library
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

#include <latan/latan_min_gsl.h>
#include <latan/latan_includes.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

#define DIFF_PREC    1.0e-8
#define FDF_CVG_PREC 1.0e-3
#define F_CVG_PREC   1.0e-6
#define INIT_STEP    1.0e-2
#define DIR_TOL      1.0e-2

static double gsl_f(const gsl_vector *v, void *gparam);
static double gsl_f_i(double v_i, void *param);
static void gsl_df(const gsl_vector *v, void *gparam, gsl_vector *df);
static void gsl_fdf(const gsl_vector *v, void *v_gf_param, double *f,\
                    gsl_vector *df);

typedef struct
{
    void *param;
    min_func *f;
    gsl_vector *scale;
    gsl_vector *buf_gsl_f;
} gsl_f_eval_param;

typedef struct
{
    gsl_f_eval_param *gf_param;
    const gsl_vector *v;
    size_t i;
} gsl_f_i_param;

double gsl_f(const gsl_vector *v, void *v_gf_param)
{
    gsl_f_eval_param *gf_param;
    mat x;
    gsl_matrix_view x_mview;
    double f_val;
    
    gf_param = (gsl_f_eval_param*)v_gf_param;
    gsl_vector_memcpy(gf_param->buf_gsl_f,v);
    gsl_vector_mul(gf_param->buf_gsl_f,gf_param->scale);
    x_mview    = gsl_matrix_view_vector(gf_param->buf_gsl_f,v->size,1);
    x.data_cpu = &(x_mview.matrix);
    f_val      = gf_param->f(&x,gf_param->param);
    
    return f_val;
}

double gsl_f_i(double v_i, void *v_gfi_param)
{
    gsl_f_i_param *gfi_param;
    gsl_vector *v_mod;
    double f_val;
    
    gfi_param = (gsl_f_i_param*)v_gfi_param;
    
    v_mod = gsl_vector_alloc(gfi_param->v->size);
    
    gsl_vector_memcpy(v_mod,gfi_param->v);
    gsl_vector_set(v_mod,gfi_param->i,v_i);
    f_val = gsl_f(v_mod,gfi_param->gf_param);
    
    gsl_vector_free(v_mod);
    
    return f_val;
}

void gsl_df(const gsl_vector *v, void *v_gf_param, gsl_vector *df)
{
    size_t i;
    double dfodvi, dummy;
    gsl_function s_gsl_f_i;
    gsl_f_i_param gfi_param;
    
    gfi_param.gf_param = (gsl_f_eval_param*)v_gf_param;
    gfi_param.v        = v;
    s_gsl_f_i.function = &gsl_f_i;
    s_gsl_f_i.params   = &gfi_param;
    
    for (i=0;i<v->size;i++)
    {
        gfi_param.i = i;
        gsl_deriv_central(&s_gsl_f_i,gsl_vector_get(v,i),DIFF_PREC,&dfodvi,\
                          &dummy);
        gsl_vector_set(df,i,dfodvi);
    }
}

void gsl_fdf(const gsl_vector *v, void *v_gf_param, double *f, gsl_vector *df)
{
    *f = gsl_f(v,v_gf_param);
    gsl_df(v,v_gf_param,df);
}

latan_errno minimize_gsl(mat *x, double *f_min, min_func *f, void *param)
{
    latan_errno status;
    size_t n;
    size_t i;
    unsigned int iter, max_iteration;
    gsl_vector *gsl_x,*step_size,*one;
    gsl_f_eval_param gf_param;
    gsl_multimin_function_fdf gsl_min_func_fdf;
    gsl_multimin_function gsl_min_func_f;
    bool need_df;
    const gsl_multimin_fdfminimizer_type *minimizer_fdf_t;
    gsl_multimin_fdfminimizer *minimizer_fdf;
    const gsl_multimin_fminimizer_type *minimizer_f_t;
    gsl_multimin_fminimizer *minimizer_f;
    strbuf buf, x_dump, war_msg;
    
    n                       = nrow(x);
    iter                    = 0u;
    max_iteration           = minimizer_get_max_iteration();
    gf_param.f              = f;
    gf_param.param          = param;
    gsl_min_func_fdf.f      = &gsl_f;
    gsl_min_func_fdf.df     = &gsl_df;
    gsl_min_func_fdf.fdf    = &gsl_fdf;
    gsl_min_func_fdf.params = &gf_param;
    gsl_min_func_fdf.n      = n;
    gsl_min_func_f.f        = &gsl_f;
    gsl_min_func_f.params   = &gf_param;
    gsl_min_func_f.n        = n;
    minimizer_fdf_t         = NULL;
    minimizer_fdf           = NULL;
    minimizer_f_t           = NULL;
    minimizer_f             = NULL;
    
    switch (minimizer_get_alg())
    {
        case GSL_GRAD_FR:
            need_df = true;
            minimizer_fdf_t = gsl_multimin_fdfminimizer_conjugate_fr;
            break;
        case GSL_GRAD_PR:
            need_df = true;
            minimizer_fdf_t = gsl_multimin_fdfminimizer_conjugate_pr;
            break;
        case GSL_VEC_BFGS:
            need_df = true;
            minimizer_fdf_t = gsl_multimin_fdfminimizer_vector_bfgs;
            break;
        case GSL_SIMPLEX_NM:
            need_df = false;
            minimizer_f_t = gsl_multimin_fminimizer_nmsimplex2;
            break;
        default:
            LATAN_ERROR("invalid GSL minimization algorithm flag",\
                        LATAN_EINVAL);
    }
    
    gf_param.scale     = gsl_vector_alloc(n);
    gf_param.buf_gsl_f = gsl_vector_alloc(n);
    gsl_x              = gsl_vector_alloc(n);
    step_size          = gsl_vector_alloc(n);
    one                = gsl_vector_alloc(n);
    
    gsl_matrix_get_col(gsl_x,x->data_cpu,0);
    gsl_vector_memcpy(gf_param.scale,gsl_x);
    gsl_vector_set_all(one,1.0);
    if (need_df)
    {
        minimizer_fdf = gsl_multimin_fdfminimizer_alloc(minimizer_fdf_t,n);
        gsl_multimin_fdfminimizer_set(minimizer_fdf,&gsl_min_func_fdf,one,\
                                      INIT_STEP,DIR_TOL);
    }
    else
    {
        minimizer_f = gsl_multimin_fminimizer_alloc(minimizer_f_t,n);
        gsl_vector_set_all(step_size,INIT_STEP);
        gsl_multimin_fminimizer_set(minimizer_f,&gsl_min_func_f,one,step_size);
    }
    do
    {
        iter++;
        if (need_df)
        {
            status = gsl_multimin_fdfminimizer_iterate(minimizer_fdf);
        }
        else
        {
            status = gsl_multimin_fminimizer_iterate(minimizer_f);
        }
        if (status)
        {
            break;
        }
        latan_printf(DEBUG2,"------ iteration %d\n",iter);
        if (need_df)
        {
            status = gsl_multimin_test_gradient(minimizer_fdf->gradient,\
                                                FDF_CVG_PREC);
            latan_printf(DEBUG2,"f(x)= %10.6f gradf(x)= %10.4e\n",\
                         minimizer_fdf->f,                       \
                         gsl_blas_dnrm2(minimizer_fdf->gradient));
            gsl_vector_memcpy(gsl_x,minimizer_fdf->x);
        }
        else
        {
            status = gsl_multimin_test_size(minimizer_f->size,F_CVG_PREC);
            latan_printf(DEBUG2,"f(x)= %10.6f size= %10.4e\n",minimizer_f->fval,\
                         minimizer_f->size);
            gsl_vector_memcpy(gsl_x,minimizer_f->x);
        }
        gsl_vector_mul(gsl_x,gf_param.scale);
        if (latan_get_verb() == DEBUG2)
        {
            strbufcpy(x_dump,"");
            for (i=0;i<n;i++)
            {
                sprintf(buf,"x_%lu= %10.6f ",(long unsigned)i,\
                gsl_vector_get(gsl_x,i));
                strbufcat(x_dump,buf);
            }
        }
        latan_printf(DEBUG2,"%s\n",x_dump);
    } while ((status == (latan_errno)GSL_CONTINUE) && (iter <= max_iteration));
    if (status)
    {
        sprintf(war_msg,"GSL warning : %s",\
                gsl_strerror(status));
    }
    *f_min = need_df ? minimizer_fdf->f : minimizer_f->fval;
    gsl_matrix_set_col(x->data_cpu,0,gsl_x);
    
    gsl_vector_free(gf_param.scale);
    gsl_vector_free(gf_param.buf_gsl_f);
    gsl_vector_free(gsl_x);
    gsl_vector_free(step_size);
    gsl_vector_free(one);
    if (need_df)
    {
        gsl_multimin_fdfminimizer_free(minimizer_fdf);
    }
    else
    {
        gsl_multimin_fminimizer_free(minimizer_f);
    }
    
    return status;
}
