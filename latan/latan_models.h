/* latan_models.h, part of LatAnalyze library
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

#ifndef LATAN_MODELS_H_
#define LATAN_MODELS_H_

#include <latan/latan_globals.h>
#include <latan/latan_fit.h>

/* 1D models */
/** 1D polynomial models **/
#define FM_POLYN_1D_MAXDEG 4
double fm_polyn_1d_0_func(const mat *X, const mat *p, void *nothing);
extern const fit_model fm_polyn_1d_0;
double fm_polyn_1d_1_func(const mat *X, const mat *p, void *nothing);
extern const fit_model fm_polyn_1d_1;
double fm_polyn_1d_2_func(const mat *X, const mat *p, void *nothing);
extern const fit_model fm_polyn_1d_2;
double fm_polyn_1d_3_func(const mat *X, const mat *p, void *nothing);
extern const fit_model fm_polyn_1d_3;
double fm_polyn_1d_4_func(const mat *X, const mat *p, void *nothing);
extern const fit_model fm_polyn_1d_4;
#define fm_lin fm_polyn_1d_1
#define fm_const fm_polyn_1d_0
extern const fit_model *fm_polyn_1d[FM_POLYN_1D_MAXDEG+1];

/** exponential decay **/
double fm_expdec_func(const mat *X, const mat *p, void *nothing);
extern const fit_model fm_expdec;

/** hyperbolic cosine **/
double fm_cosh_func(const mat *X, const mat *p, void *nothing);
extern const fit_model fm_cosh;

/* 2D models */
/** 2D polynomial models **/
#define FM_POLYN_2D_MAXDEG 2
double fm_polyn_2d_00_func(const mat *X, const mat *p,void *nothing);
extern const fit_model fm_polyn_2d_00;
double fm_polyn_2d_01_func(const mat *X, const mat *p,void *nothing);
extern const fit_model fm_polyn_2d_01;
double fm_polyn_2d_02_func(const mat *X, const mat *p,void *nothing);
extern const fit_model fm_polyn_2d_02;
double fm_polyn_2d_10_func(const mat *X, const mat *p,void *nothing);
extern const fit_model fm_polyn_2d_10;
double fm_polyn_2d_11_func(const mat *X, const mat *p,void *nothing);
extern const fit_model fm_polyn_2d_11;
double fm_polyn_2d_12_func(const mat *X, const mat *p,void *nothing);
extern const fit_model fm_polyn_2d_12;
double fm_polyn_2d_20_func(const mat *X, const mat *p,void *nothing);
extern const fit_model fm_polyn_2d_20;
double fm_polyn_2d_21_func(const mat *X, const mat *p,void *nothing);
extern const fit_model fm_polyn_2d_21;
double fm_polyn_2d_22_func(const mat *X, const mat *p,void *nothing);
extern const fit_model fm_polyn_2d_22;
extern const fit_model *fm_polyn_2d[FM_POLYN_2D_MAXDEG+1][FM_POLYN_2D_MAXDEG+1];

#endif
