/* latan_fit.h, part of LatAnalyze library
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

#ifndef LATAN_FIT_H_
#define LATAN_FIT_H_

#include <latan/latan_globals.h>
#include <latan/latan_minimizer.h>
#include <latan/latan_statistics.h>

#ifndef MAX_YDIM
#define MAX_YDIM 16
#endif

__BEGIN_DECLS

/* fit model structure */
typedef double model_func(const mat *x, const mat *p, void *model_param);
typedef size_t npar_func(void *model_param);
typedef void plot2dstr_func(strbuf str, const size_t kx, const mat *x,\
                            const mat *p, void *model_param);

typedef struct
{
    strbuf name;
    model_func *func[MAX_YDIM];
    npar_func *npar;
    plot2dstr_func *plot2dstr[MAX_YDIM];
    size_t nxdim;
    size_t nydim;
} fit_model;

/** some useful constant npar_func **/
npar_func npar_1;
npar_func npar_2;
npar_func npar_3;
npar_func npar_4;
npar_func npar_5;
npar_func npar_6;
npar_func npar_7;
npar_func npar_8;
npar_func npar_9;
npar_func npar_10;

/** access **/
void fit_model_get_name(strbuf name, const fit_model *model);
size_t fit_model_get_npar(const fit_model *model, void *model_param);
double fit_model_eval(const fit_model *model, const size_t k, const mat *x,\
                      const mat *p, void *model_param);
void fit_model_plot2dstr(strbuf str, const fit_model *model, const size_t kx,\
                         const size_t ky, const mat *x, const mat *p,        \
                         void *model_param);

/* fit data structure */
/** chi^2 buffer **/
typedef struct chi2_buf_s
{
    mat *x;
    mat *Y;
    mat *CyY;
    mat *X;
    mat *CxX;
    mat *lX;
    mat *ClX;
    bool is_xpart_alloc;
} chi2_buf;

/** the main structure **/
typedef struct fid_data_s
{
    /* sizes */
    size_t ndata;
    size_t nxdim;
    size_t nydim;
    size_t npar;
    size_t ndpar;
    /* point matrices */
    mat *x;
    mat **x_covar;
    mat *x_var_inv;
    bool is_x_correlated;
    bool *have_x_covar;
    bool *to_fit;
    /* data matrices */
    mat *y;
    mat **y_covar;
    mat *y_var_inv;
    bool is_y_correlated;
    /* data/point covariance matrix */
    mat **xy_covar;
    bool *have_xy_covar;
    /* correlation filter */
    mat *cor_filter;
    /* inverse variance matrix */
    mat *var_inv;
    /* is everything ready to perform a fit ? */
    bool is_inverted;
    /* fit model */
    fit_model const *model;
    void *model_param;
    /* chi^2 */
    min_func *chi2_ext;
    double chi2pdof;
    bool save_chi2pdof;
    chi2_buf *buf;
    int nbuf;
    /* sample counter */
    size_t s;
} fit_data;

/** allocation **/
fit_data *fit_data_create(const size_t ndata, const size_t nxdim,\
                          const size_t nydim);
void fit_data_destroy(fit_data *d);

/** access **/
/*** sizes ***/
size_t fit_data_get_ndata(const fit_data *d);
size_t fit_data_get_nydim(const fit_data *d);
size_t fit_data_get_nxdim(const fit_data *d);
size_t fit_data_get_npar(const fit_data *d);
size_t fit_data_get_ndumbpar(const fit_data *d);
void fit_data_set_ndumbpar(fit_data *d, const size_t ndpar);

/*** chi^2 ***/
void fit_data_save_chi2pdof(fit_data *d, bool save);
double fit_data_get_chi2pdof(const fit_data *d);
void fit_data_set_chi2_ext(fit_data *d, min_func *f);

/*** data ***/
void fit_data_set_y(fit_data *d, const size_t i, const size_t k,\
                    const double y_ik);
latan_errno fit_data_set_y_k(fit_data *d, const size_t i, const mat *y_i);
double fit_data_get_y(const fit_data *d, const size_t i, const size_t k);
latan_errno fit_data_get_y_k(mat *y_i, const fit_data *d, const size_t i);
mat * fit_data_pt_y(const fit_data *d);
latan_errno fit_data_set_y_covar(fit_data *d, const size_t k1,  \
                                 const size_t k2, const mat *var);
const mat * fit_data_pt_y_covar(const fit_data *d, const size_t k1,\
                                const size_t k2);
bool fit_data_is_y_correlated(const fit_data *d);
latan_errno fit_data_set_xy_covar(fit_data *d, const size_t ky, \
                                  const size_t kx, const mat *covar);
bool fit_data_have_xy_covar(const fit_data *d);

/*** fit points ***/
void fit_data_set_x(fit_data *d, const size_t i, const size_t k,\
                    const double x_ik);
latan_errno fit_data_set_x_k(fit_data *d, const size_t i, const mat *x_i);
double fit_data_get_x(const fit_data *d, const size_t i, const size_t k);
latan_errno fit_data_get_x_k(mat *x, const fit_data *d, const size_t i);
mat * fit_data_pt_x(const fit_data *d);
latan_errno fit_data_set_x_covar(fit_data *d, const size_t k1,  \
                                 const size_t k2, const mat *var);
const mat * fit_data_pt_x_covar(const fit_data *d, const size_t k1,\
                                const size_t k2);
bool fit_data_have_x_covar(const fit_data *d, const size_t j);
bool fit_data_have_x_var(const fit_data *d);
bool fit_data_is_x_correlated(const fit_data *d);
void fit_data_fit_all_points(fit_data *d, bool fit);
void fit_data_fit_point(fit_data *d, size_t i, bool fit);
void fit_data_fit_range(fit_data *d, size_t start, size_t end, bool fit);
void fit_data_fit_region(fit_data *d, double **xb);
bool fit_data_is_fit_point(const fit_data *d, size_t i);
size_t fit_data_fit_point_num(const fit_data *d);

/*** correlation filter ***/
latan_errno fit_data_set_data_cor(fit_data *d, const size_t i, const size_t j,\
                                  bool is_cor);
bool fit_data_is_data_cor(const fit_data *d, const size_t i, const size_t j);

/*** fit model ***/
latan_errno fit_data_set_model(fit_data *d, const fit_model *model,\
                               void *model_param);
const void * fit_data_pt_model_param(const fit_data *d);
double fit_data_model_xeval(const fit_data *d, const size_t k, const mat *x,\
                            const mat *p);
double fit_data_model_eval(const fit_data *d, const size_t k, const size_t i,\
                           const mat *p);
void fit_data_plot2dstr(strbuf str, const fit_data *d, const size_t ky,\
                        const mat *x, const size_t kx, const mat *p);

/*** dof ***/
size_t fit_data_get_dof(const fit_data *d);

/*** sample counter ***/
size_t fit_data_get_sample_counter(const fit_data *d);

/*** set from samples ***/
typedef enum
{
    NO_COR    = 0,
    DATA_COR  = 1 << 0,
    X_COR     = 1 << 1,
    XDATA_COR = 1 << 2
} cor_flag;

latan_errno fit_data_set_covar_from_sample(fit_data *d, rs_sample * const *x,\
                                           rs_sample * const *data,          \
                                           const cor_flag flag,              \
                                           const bool *use_x_var);

/* chi2 functions, have min_func type */
double chi2(mat *p, void *vd);

/* fit functions */
latan_errno data_fit(mat *p, fit_data *d);
latan_errno rs_data_fit(rs_sample *p, rs_sample * const *x,         \
                          rs_sample * const *data, fit_data *d,       \
                          const cor_flag flag, const bool *use_x_var);
void fit_residual(mat *res, const fit_data *d, const size_t ky, const mat *p);
void fit_partresidual(mat *res, const fit_data *d, const size_t ky, \
                      const mat *x_ex, const size_t kx, const mat *p);

__END_DECLS

#endif
