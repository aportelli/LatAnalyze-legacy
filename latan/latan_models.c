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

static double fm_expdec_ex_func(const mat *x, const mat *p,\
                                 void *nothing __dumb)
{
    double res,m1,m2,A1,A2,t;
    
    t   = mat_get(x,0,0);
    m1  = mat_get(p,0,0);
    m2  = mat_get(p,1,0);
    A1  = mat_get(p,2,0);
    A2  = mat_get(p,3,0);
    res = exp(-m1*t+A1)+exp(-m2*t+A2);
    
    return res;
}

fit_model fm_expdec_ex = 
{
    "y(x) = exp(-p0*x+p2) + exp(-p1*x+p3))",
    {&fm_expdec_ex_func},
    &npar_4,
    1,
    1
};

static double fm_expdec_splitsum_func0(const mat *x, const mat *p,\
                                       void *nothing __dumb)
{
    double res,m,dm,A0,t;
    
    t   = mat_get(x,0,0);
    m   = mat_get(p,0,0);
    dm  = mat_get(p,1,0);
    A0  = mat_get(p,2,0);
    res = exp(-(m+0.5*dm)*t+A0);
    
    return res;
}

static double fm_expdec_splitsum_func1(const mat *x, const mat *p,\
                                       void *nothing __dumb)
{
    double res,m,dm,A1,t;
    
    t   = mat_get(x,0,0);
    m   = mat_get(p,0,0);
    dm  = mat_get(p,1,0);
    A1  = mat_get(p,3,0);
    res = exp(-(m-0.5*dm)*t+A1);
    
    return res;
}

fit_model fm_expdec_splitsum = 
{
    "y0(x) = exp(-(p0-p1)*x+p2), y1(x) = exp(-(p0+p1)*x+p3)",
    {&fm_expdec_splitsum_func0,&fm_expdec_splitsum_func1},
    &npar_4,
    1,
    2
};

static double fm_expdec_ex_splitsum_func0(const mat *x, const mat *p,\
                                          void *nothing __dumb)
{
    double res,m,dm,A1_0,A2_0,E2_0,t;
    
    t    = mat_get(x,0,0);
    m    = mat_get(p,0,0);
    dm   = mat_get(p,1,0);
    E2_0 = mat_get(p,2,0);
    A1_0 = mat_get(p,4,0);
    A2_0 = mat_get(p,6,0);
    
    res = exp(-(m+0.5*dm)*t+A1_0) + exp(-E2_0*t+A2_0);
    
    return res;
}

static double fm_expdec_ex_splitsum_func1(const mat *x, const mat *p,\
                                          void *nothing __dumb)
{
    double res,m,dm,A1_1,A2_1,E2_1,t;
    
    t    = mat_get(x,0,0);
    m    = mat_get(p,0,0);
    dm   = mat_get(p,1,0);
    E2_1 = mat_get(p,3,0);
    A1_1 = mat_get(p,5,0);
    A2_1 = mat_get(p,7,0);
    
    res = exp(-(m-0.5*dm)*t+A1_1) + exp(-E2_1*t+A2_1);
    
    return res;
}

fit_model fm_expdec_ex_splitsum = 
{
    "y0(x) = exp(-(p0-p1)*x+p3)+exp(-p2*x+p5), y1(x) = exp(-(p0+p1)*x+p4)+exp(-p2*x+p6)",
    {&fm_expdec_ex_splitsum_func0,&fm_expdec_ex_splitsum_func1},
    &npar_8,
    1,
    2
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

static double fm_cosh_ex_func(const mat *x, const mat *p, void *vnt)
{
    double res,m1,A1,m2,A2,t;
    size_t nt;
    
    t   = mat_get(x,0,0);
    m1  = mat_get(p,0,0);
    m2  = mat_get(p,1,0);
    A1  = mat_get(p,2,0);
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

fit_model fm_cosh_ex =
{
    "y(x) = exp(p2)*cosh(p0*(x-nt/2)) + exp(p3)*cosh(p1*(x-nt/2))",
    {&fm_cosh_ex_func},
    &npar_4,
    1,
    1
};

static double fm_cosh_splitsum_func0(const mat *x, const mat *p,\
                                      void *vnt)
{
    double res,m,dm,A0,t;
    size_t nt;
    
    t   = mat_get(x,0,0);
    m   = mat_get(p,0,0);
    dm  = mat_get(p,1,0);
    A0  = mat_get(p,2,0);
    
    if (vnt)
    {
        nt = *((size_t *)(vnt));
    }
    else
    {
        nt = 0;
    }
    
    res = exp(A0)*cosh((m+0.5*dm)*(t-DRATIO(nt,2)));
    
    return res;
}

static double fm_cosh_splitsum_func1(const mat *x, const mat *p,\
                                     void *vnt)
{
    double res,m,dm,A1,t;
    size_t nt;
    
    t   = mat_get(x,0,0);
    m   = mat_get(p,0,0);
    dm  = mat_get(p,1,0);
    A1  = mat_get(p,3,0);
    
    if (vnt)
    {
        nt = *((size_t *)(vnt));
    }
    else
    {
        nt = 0;
    }
    
    res = exp(A1)*cosh((m-0.5*dm)*(t-DRATIO(nt,2)));
    
    return res;
}

fit_model fm_cosh_splitsum = 
{
    "y0(x) = exp(p2)*cosh((p0-p1)*(x-nt/2)), y1(x) = exp(p3)*cosh((p0+p1)*(x-nt/2))",
    {&fm_cosh_splitsum_func0,&fm_cosh_splitsum_func1},
    &npar_4,
    1,
    2
};

static double fm_cosh_ex_splitsum_func0(const mat *x, const mat *p,\
                                        void *vnt)
{
    double res,m,dm,A1_0,A2_0,E2_0,t;
    size_t nt;
    
    t    = mat_get(x,0,0);
    m    = mat_get(p,0,0);
    dm   = mat_get(p,1,0);
    E2_0 = mat_get(p,2,0);
    A1_0 = mat_get(p,4,0);
    A2_0 = mat_get(p,6,0);
    
    if (vnt)
    {
        nt = *((size_t *)(vnt));
    }
    else
    {
        nt = 0;
    }
    
    res = exp(A1_0)*cosh((m+0.5*dm)*(t-DRATIO(nt,2)))\
          +exp(A2_0)*cosh(E2_0*(t-DRATIO(nt,2)));
    
    return res;
}

static double fm_cosh_ex_splitsum_func1(const mat *x, const mat *p,\
                                        void *vnt)
{
    double res,m,dm,A1_1,A2_1,E2_1,t;
    size_t nt;
    
    t    = mat_get(x,0,0);
    m    = mat_get(p,0,0);
    dm   = mat_get(p,1,0);
    E2_1 = mat_get(p,3,0);
    A1_1 = mat_get(p,5,0);
    A2_1 = mat_get(p,7,0);
    
    if (vnt)
    {
        nt = *((size_t *)(vnt));
    }
    else
    {
        nt = 0;
    }
    
    res = exp(A1_1)*cosh((m-0.5*dm)*(t-DRATIO(nt,2)))\
          +exp(A2_1)*cosh(E2_1*(t-DRATIO(nt,2)));
    
    return res;
}

fit_model fm_cosh_ex_splitsum = 
{
    "y0(x) = exp(p3)*cosh((p0-p1)*(x-nt/2))+exp(p5)*cosh(p2*(x-nt/2)), y1(x) = exp(p4)*cosh((p0+p1)*(x-nt/2))+exp(p6)*cosh(p2*(x-nt/2))",
    {&fm_cosh_ex_splitsum_func0,&fm_cosh_ex_splitsum_func1},
    &npar_8,
    1,
    2
};
