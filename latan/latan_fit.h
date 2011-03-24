/* latan_fit.h, part of LatAnalyze library
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

#ifndef LATAN_FIT_H_
#define LATAN_FIT_H_

#include <latan/latan_globals.h>
#include <latan/latan_statistics.h>

#define MAX_STAGE (sizeof(unsigned int)*8)
#define STAGE(s) (1 << (s))

__BEGIN_DECLS

/* fit model structure */
typedef double model_func(const mat *x, const mat *p, void *model_param);
typedef size_t npar_func(void *model_param);

typedef struct
{
    strbuf name;
    model_func *func[MAX_STAGE];
    npar_func *npar[MAX_STAGE];
    size_t ndim;
    size_t nstage;
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
size_t fit_model_get_npar(const fit_model *model,                     \
                          const unsigned int stage_flag, void *model_param);
double fit_model_eval(const fit_model *model, const mat *x,       \
                      const mat *p, const unsigned int stage_flag,\
                      void *model_param);

/* fit data structure */
/** chi^2 buffer **/
typedef struct
{
    mat *x;
    mat *X;
    mat *CdX;
    mat *Y;
    mat *CxY;
    mat *lX;
    mat *ClX;
    bool is_xpart_alloc;
} chi2_buf;

/** boolean array for stages **/
typedef bool stage_ar[MAX_STAGE];

/** the main structure **/
typedef struct
{
    /* sizes */
    size_t ndata;
    size_t ndim;
    size_t npar;
    /* point matrices */
    mat *x;
    mat **x_covar;
    mat *x_var_inv;
    bool is_x_correlated;
    bool *have_x_covar;
    bool *to_fit;
    /* data matrices */
    mat *data;
    mat *data_var;
    mat *data_var_inv;
    bool is_data_correlated;
    /* data/point covariance matrix */
    mat **xdata_covar;
    bool *have_xdata_covar;
    /* inverse variance matrix */
    mat *var_inv;
    /* is everything ready to perform a fit ? */
    bool is_inverted;
    /* fit model */
    const fit_model *model;
    void *model_param;
    /* stages */
    unsigned int stage_flag;
    /* buffers for chi^2 computation */
    double chi2pdof;
    bool save_chi2pdof;
    chi2_buf *buf;
    int nbuf;
} fit_data;



/** allocation **/
fit_data *fit_data_create(const size_t ndata, const size_t ndim);
void fit_data_destroy(fit_data *d);

/** access **/
/*** sizes ***/
size_t fit_data_get_ndata(const fit_data *d);
size_t fit_data_get_ndim(const fit_data *d);

/*** chi2 value ***/
void fit_data_save_chi2pdof(fit_data *d, bool save);
double fit_data_get_chi2pdof(const fit_data *d);

/*** fit points ***/
void fit_data_set_x(fit_data *d, const size_t i, const size_t j,\
                    const double x_i);
latan_errno fit_data_set_x_from_mat(fit_data *d, const size_t j, const mat *m);
double fit_data_get_x(const fit_data *d, const size_t i, const size_t j);
mat *fit_data_pt_x(const fit_data *d);
latan_errno fit_data_set_x_covar(fit_data *d, const size_t i1,  \
                                 const size_t i2, const mat *var);
bool fit_data_have_x_covar(const fit_data *d, const size_t j);
bool fit_data_have_x_var(const fit_data *d);
bool fit_data_is_x_correlated(const fit_data *d);
void fit_data_fit_all_points(fit_data *d, bool fit);
void fit_data_fit_point(fit_data *d, size_t i, bool fit);
void fit_data_fit_range(fit_data *d, size_t start, size_t end, bool fit);
bool fit_data_is_fit_point(const fit_data *d, size_t i);
size_t fit_data_fit_point_num(const fit_data *d);

/*** data ***/
void fit_data_set_data(fit_data *d, const size_t i, const double data_i);
double fit_data_get_data(const fit_data *d, const size_t i);
mat *fit_data_pt_data(const fit_data *d);
latan_errno fit_data_set_data_var(fit_data *d, const mat *var);
bool fit_data_is_data_correlated(const fit_data *d);
latan_errno fit_data_set_xdata_covar(fit_data *d, const size_t j,\
                                     const mat *covar);
bool fit_data_have_xdata_covar(const fit_data *d);

/*** fit model ***/
latan_errno fit_data_set_model(fit_data *d, const fit_model *model,\
                               void *model_param);
void fit_data_set_model_param(fit_data *d, void *model_param);
double fit_data_model_eval(const fit_data *d, const size_t i, const mat *p);

/*** stages ***/
void fit_data_set_stage_flag(fit_data *d, const unsigned int stage_flag);
unsigned int fit_data_get_stage_flag(const fit_data *d);
size_t fit_data_get_npar(const fit_data *d);
void fit_data_set_stages(fit_data *d, const stage_ar s);
void fit_data_get_stages(stage_ar s, const fit_data *d);

/*** dof ***/
size_t fit_data_get_dof(const fit_data *d);

/* chi2 functions, have min_func type */
double chi2(mat *var, void *d);

/* fit functions */
typedef enum
{
    NO_COR    = 0,
    DATA_COR  = 1 << 0,
    X_COR     = 1 << 1,
    XDATA_COR = 1 << 2
} cor_flag;

latan_errno data_fit(mat *p, fit_data *d);
latan_errno rs_data_fit(rs_sample *p, const rs_sample *data, fit_data *d,\
                        const cor_flag flag);
latan_errno rs_x_data_fit(rs_sample *p, rs_sample * const *x,         \
                          const rs_sample *data, fit_data *d,         \
                          const cor_flag flag, const bool *use_x_var);
void fit_residual(mat *res, const mat *p, const fit_data *d);
void rs_fit_residual(rs_sample *res, const rs_sample *p, const rs_sample *data,\
                     fit_data *d);
void rs_x_fit_residual(rs_sample *res, const rs_sample *p,         \
                       rs_sample * const *x, const rs_sample *data,\
                       fit_data *d);
__END_DECLS

#endif
