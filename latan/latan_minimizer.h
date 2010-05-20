#ifndef LATAN_MINIMIZER_H_
#define LATAN_MINIMIZER_H_

#include <latan/latan_globals.h>

__BEGIN_DECLS

/* prototype of function to minimize */
typedef double min_func(const mat var, void* param);

/* the minimizer */
latan_errno minimize(mat var, double* f_min, min_func* f, void* param);

/* fit model structure */
typedef double model_func(const double x, const mat func_param,\
						  void* model_param);

typedef struct
{
	stringbuf name;
	model_func* func;
	size_t npar;
	stringbuf plot_fmt;
} fit_model;

/** access **/
void fit_model_get_name(stringbuf name, const fit_model* model);
void fit_model_get_plot_fmt(stringbuf plot_fmt, const fit_model* model);
double fit_model_eval(const fit_model* model, const double x,\
					  const mat func_param, void* model_param);

/** some useful models **/
/*** constant: y(x) = p0 ***/
double fm_const_func(const double x, const mat func_param, void* nothing);
extern const fit_model fm_const;
/*** linear: y(x) = p0 + p1*x ***/
double fm_lin_func(const double x, const mat func_param, void* nothing);
extern const fit_model fm_lin;
/*** exponential decay: y(x) = p0*exp(-p1*x) ***/
double fm_expdec_func(const double x, const mat func_param, void* nothing);
extern const fit_model fm_expdec;
/*** hyperbolic cosine: y(x) = p0*cosh(p1*x) ***/
double fm_cosh_func(const double x, const mat func_param, void* nothing);
extern const fit_model fm_cosh;

/* fit data structure */
typedef struct
{
	size_t ndata;
	mat x;
	mat data;
	mat var_inveigval;
	mat var_eigvec;
	bool is_correlated;
	bool save_chi2pdof;
	bool* to_fit;
	const fit_model* model;
	void* model_param;
	int stage;
	double chi2pdof;
}* fit_data;

/** allocation **/
fit_data fit_data_create(const size_t ndata);
void fit_data_destroy(fit_data d);

/** access **/
void fit_data_save_chi2pdof(fit_data d, bool save);
double fit_data_get_chi2pdof(fit_data d);
void fit_data_set_x(fit_data d, const size_t i, const double x_i);
double fit_data_get_x(const fit_data d, const size_t i);
void fit_data_fit_all_points(fit_data d, bool fit);
void fit_data_fit_point(fit_data d, size_t i, bool fit);
void fit_data_fit_range(fit_data d, size_t start, size_t end, bool fit);
bool fit_data_is_fit_point(fit_data d, size_t i);
size_t fit_data_fit_point_num(fit_data d);
mat fit_data_pt_x(const fit_data d);
void fit_data_set_data(fit_data d, const size_t i, const double data_i);
double fit_data_get_data(const fit_data d, const size_t i);
mat fit_data_pt_data(const fit_data d);
latan_errno fit_data_set_var(fit_data d, const mat var);
void fit_data_set_model(fit_data d, const fit_model* model);
const fit_model* fit_data_pt_model(fit_data d);
void fit_data_set_model_param(fit_data d, void* model_param);
double fit_data_model_eval(const fit_data d, const size_t i,\
						   const mat func_param);
void fit_data_set_stage(fit_data d, const int stage);
int fit_data_get_stage(const fit_data d);
int fit_data_get_dof(const fit_data d);
bool fit_data_is_correlated(const fit_data d);

/* chi2 functions, have min_func type */
double chi2(const mat var, void* d);

/* fit functions */
latan_errno data_fit(mat fit_param, fit_data d);

__END_DECLS

#endif