#ifndef LATAN_FIT_H_
#define LATAN_FIT_H_

#include <latan/latan_globals.h>
#include <latan/latan_statistics.h>

__BEGIN_DECLS

/* fit model structure */
typedef double model_func(mat *x, mat *p, void *model_param);

typedef struct
{
	strbuf name;
	model_func *func;
	size_t npar;
	size_t ndim;
	strbuf plot_fmt;
} fit_model;

/** access **/
void fit_model_get_name(strbuf name, const fit_model *model);
void fit_model_get_plot_fmt(strbuf plot_fmt, const fit_model *model);
double fit_model_eval(const fit_model *model, mat *x,mat *p,\
					  void *model_param);

/* fit data structure */
typedef struct
{
	size_t ndata;
	size_t ndim;
	mat *x;
	mat *x_varinv;
	mat *data;
	mat *data_varinv;
	bool is_data_correlated;
	bool is_x_correlated;
	bool have_x_var;
	bool save_chi2pdof;
	bool *to_fit;
	const fit_model *model;
	void *model_param;
	int stage;
	double chi2pdof;
	mat *buf_chi2[4];
} fit_data;

/** allocation **/
fit_data *fit_data_create(const size_t ndata, const size_t ndim);
void fit_data_destroy(fit_data *d);

/** access **/
/*** chi2 value ***/
void fit_data_save_chi2pdof(fit_data *d, bool save);
double fit_data_get_chi2pdof(fit_data *d);

/*** fit points ***/
void fit_data_set_x(fit_data *d, const size_t i, const size_t j,\
					const double x_i);
double fit_data_get_x(const fit_data *d, const size_t i, const size_t j);
mat *fit_data_pt_x(const fit_data *d);
latan_errno fit_data_set_x_var(fit_data *d, mat *var);
bool fit_data_have_x_var(fit_data *d);
void fit_data_fit_all_points(fit_data *d, bool fit);
void fit_data_fit_point(fit_data *d, size_t i, bool fit);
void fit_data_fit_range(fit_data *d, size_t start, size_t end, bool fit);
bool fit_data_is_fit_point(const fit_data *d, size_t i);
size_t fit_data_fit_point_num(const fit_data *d);

/*** data ***/
void fit_data_set_data(fit_data *d, const size_t i, const double data_i);
double fit_data_get_data(const fit_data *d, const size_t i);
size_t fit_data_get_ndata(const fit_data *d);
mat *fit_data_pt_data(const fit_data *d);
latan_errno fit_data_set_data_var(fit_data *d, mat *var);
bool fit_data_is_correlated(const fit_data *d);

/*** fit model ***/
latan_errno fit_data_set_model(fit_data *d, const fit_model *model);
const fit_model *fit_data_pt_model(fit_data *d);
void fit_data_set_model_param(fit_data *d, void *model_param);
double fit_data_model_eval(const fit_data *d, const size_t i,mat *p);

/*** stages ***/
void fit_data_set_stage(fit_data *d, const int stage);
int fit_data_get_stage(const fit_data *d);

/*** dof ***/
int fit_data_get_dof(const fit_data *d);

/* chi2 functions, have min_func type */
double chi2(mat *var, void *d);

/* fit functions */
latan_errno data_fit(mat *p, fit_data *d);
latan_errno rs_data_fit(rs_sample *p, rs_sample *data, fit_data *d);
latan_errno rs_x_data_fit(rs_sample *p, rs_sample *x, rs_sample *data,\
						  fit_data *d);

__END_DECLS

#endif