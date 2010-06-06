#ifndef LATAN_MODELS_H_
#define LATAN_MODELS_H_

#include <latan/latan_globals.h>
#include <latan/latan_minimizer.h>

/*** constant: y(x) = p0 ***/
double fm_const_func(const mat x, const mat func_param, void* nothing);
extern const fit_model fm_const;

/*** linear: y(x) = p0 + p1*x ***/
double fm_lin_func(const mat x, const mat func_param, void* nothing);
extern const fit_model fm_lin;

/*** exponential decay: y(x) = p0*exp(-p1*x) ***/
double fm_expdec_func(const mat x, const mat func_param, void* nothing);
extern const fit_model fm_expdec;

/*** hyperbolic cosine: y(x) = p0*cosh(p1*x) ***/
double fm_cosh_func(const mat x, const mat func_param, void* nothing);
extern const fit_model fm_cosh;

/*** 2D polynomial models ***/
double fm_polyn_2d_00_func(const mat X, const mat func_param,void* nothing);
extern const fit_model fm_polyn_2d_00;
double fm_polyn_2d_01_func(const mat X, const mat func_param,void* nothing);
extern const fit_model fm_polyn_2d_01;
double fm_polyn_2d_02_func(const mat X, const mat func_param,void* nothing);
extern const fit_model fm_polyn_2d_02;
double fm_polyn_2d_10_func(const mat X, const mat func_param,void* nothing);
extern const fit_model fm_polyn_2d_10;
double fm_polyn_2d_11_func(const mat X, const mat func_param,void* nothing);
extern const fit_model fm_polyn_2d_11;
double fm_polyn_2d_12_func(const mat X, const mat func_param,void* nothing);
extern const fit_model fm_polyn_2d_12;
double fm_polyn_2d_20_func(const mat X, const mat func_param,void* nothing);
extern const fit_model fm_polyn_2d_20;
double fm_polyn_2d_21_func(const mat X, const mat func_param,void* nothing);
extern const fit_model fm_polyn_2d_21;
double fm_polyn_2d_22_func(const mat X, const mat func_param,void* nothing);
extern const fit_model fm_polyn_2d_22;
extern const fit_model* fm_polyn_2d[3][3];

#endif
