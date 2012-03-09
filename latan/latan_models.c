/* latan_models.c, part of LatAnalyze library
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

#include <latan/latan_models.h>
#include <latan/latan_includes.h>
#include <latan/latan_math.h>

/*                              1D models                                   */
/****************************************************************************/
/** 1D polynomial models **/
static double fm_const_func(const mat *X, const mat *p, void *nothing)
{
    const mat *dummy;

    dummy   = X;
    nothing = NULL;
    
    return mat_get(p,0,0);
}


const fit_model fm_const =
{
    "y(x) = p0",
    {&fm_const_func},
    &npar_1,
    1,
    1
};

/** exponential decay **/
static double fm_expdec_func(const mat *x, const mat *p, void *nothing)
{
    double res;
    
    nothing = NULL;
    res = mat_get(p,1,0)*exp(-mat_get(p,0,0)*mat_get(x,0,0));
    
    return res;
}

{

const fit_model fm_expdec = 
{
    "y(x) = p1*exp(-p0*x)",
    {&fm_expdec_func},
    &npar_2,
    1,
    1
};

/** hyperbolic cosine **/
static double fm_cosh_func(const mat *x, const mat *p, void *nothing)
{
    double res;
    
    nothing = NULL;
    res = mat_get(p,1,0)*cosh(mat_get(p,0,0)*mat_get(x,0,0));
    
    return res;
}

const fit_model fm_cosh =
{
    "y(x) = p1*cosh(p0*x)",
    {&fm_cosh_func},
    &npar_2,
    1,
    1
};

/*                              2D models                                   */
/****************************************************************************/
