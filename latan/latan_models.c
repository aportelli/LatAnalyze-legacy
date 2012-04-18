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
static double fm_const_func(const mat *X __dumb, const mat *p,\
                            void *nothing __dumb)
{
    return mat_get(p,0,0);
}


fit_model fm_const =
{
    "y(x) = p0",
    {&fm_const_func},
    &npar_1,
    1,
    1
};

/** exponential decay **/
static double fm_expdec_func(const mat *x, const mat *p, void *nothing __dumb)
{
    double res,m1,A1,t;

    t   = mat_get(x,0,0);
    m1  = mat_get(p,0,0);
    A1  = mat_get(p,1,0);
    res = exp(-m1*t+A1);
    
    return res;
}

fit_model fm_expdec = 
{
    "y(x) = exp(-p0*x+p1)",
    {&fm_expdec_func},
    &npar_2,
    1,
    1
};

static double fm_dbl_expdec_func(const mat *x, const mat *p,\
                                 void *nothing __dumb)
{
    double res,m1,m2,A1,A2,t;
    
    t   = mat_get(x,0,0);
    m1  = mat_get(p,0,0);
    A1  = mat_get(p,1,0);
    m2  = m1*exp(SQ(mat_get(p,2,0)));
    A2  = mat_get(p,3,0);
    res = exp(-m1*t+A1)+exp(-m2*t+A2);
    
    return res;
}

fit_model fm_dbl_expdec = 
{
    "y(x) = exp(-p0*x+p1)+exp(-p0*exp(p2^2)*x+p3)",
    {&fm_dbl_expdec_func},
    &npar_4,
    1,
    1
};

/** hyperbolic cosine **/
static double fm_cosh_func(const mat *x, const mat *p, void *vnt)
{
    double res,m1,A1,t;
    size_t nt;
    
    t   = mat_get(x,0,0);
    m1  = mat_get(p,0,0);
    A1  = mat_get(p,1,0);
    
    if (vnt)
    {
        nt = *((size_t *)(vnt));
    }
    else
    {
        nt = 0;
    }
    res = exp(A1)*cosh(m1*(t-DRATIO(nt,2)));
    
    return res;
}

fit_model fm_cosh =
{
    "y(x) = exp(p1)*cosh(p0*(x-nt/2))",
    {&fm_cosh_func},
    &npar_2,
    1,
    1
};

static double fm_dbl_cosh_func(const mat *x, const mat *p, void *vnt)
{
    double res,m1,A1,m2,A2,t;
    size_t nt;
    
    t   = mat_get(x,0,0);
    m1  = mat_get(p,0,0);
    A1  = mat_get(p,1,0);
    m2  = m1*exp(SQ(mat_get(p,2,0)));
    A2  = mat_get(p,3,0);
    
    if (vnt)
    {
        nt = *((size_t *)(vnt));
    }
    else
    {
        nt = 0;
    }
    res = exp(A1)*cosh(m1*(t-DRATIO(nt,2)))+exp(A2)*cosh(m2*(t-DRATIO(nt,2)));
    
    return res;
}

fit_model fm_dbl_cosh =
{
    "y(x) = exp(p1)*cosh(p0*(x-nt/2))+exp(p3)*cosh(p0*exp(p2^2)*(x-nt/2))",
    {&fm_dbl_cosh_func},
    &npar_4,
    1,
    1
};

/*                              2D models                                   */
/****************************************************************************/
