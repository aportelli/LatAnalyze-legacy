#ifndef LATAN_MINIMIZER_H_
#define LATAN_MINIMIZER_H_

#include <latan/latan_globals.h>

__BEGIN_DECLS

/* prototype of function to minimize */
typedef double min_func(const mat var, void* param);

/* the minimizer */
latan_errno minimize(mat var, double* f_min, min_func* f, void* param);

/* fit data structure */
typedef double model_func(const double x, const mat fit_param,\
						  void* model_param);

typedef struct
{
	mat x;
	mat data;
	mat var_inveigval;
	mat var_eigvec;
	bool is_correlated;
	model_func* model;
	void* model_param;
	int stage;
}* fit_data;

/** allocation **/
fit_data fit_data_create(mat x, mat data, mat var, model_func* model,\
						 void* model_param);
void fit_data_destroy(fit_data d);

/* chi2 functions, got min_func type */
double chi2(const mat var, void* d);

__END_DECLS

#endif