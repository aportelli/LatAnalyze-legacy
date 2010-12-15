#include <latan/latan_models.h>
#include <latan/latan_includes.h>
#include <latan/latan_math.h>

/*                              1D models                                   */
/****************************************************************************/
/** 1D polynomial models **/
double fm_polyn_1d_0_func(mat *X, mat *p, void *nothing)
{
    mat *dummy;
    
    nothing = NULL;
    dummy   = X;
    
    return mat_get(p,0,0);
}

const fit_model fm_polyn_1d_0 =
{
    "y(x) = p0",
    {fm_polyn_1d_0_func},
    &npar_1,
    1,
    "%e"
};

double fm_polyn_1d_1_func(mat *X, mat *p, void *nothing)
{
    double x,res;
    
    nothing = NULL;
    x = mat_get(X,0,0);
    
    res = mat_get(p,0,0) + mat_get(p,1,0)*x;
    
    return res;
}

const fit_model fm_polyn_1d_1 =
{
    "y(x) = p0 + p1*x",
    {fm_polyn_1d_1_func},
    &npar_2,
    1,
    "%e+%e*x"
};

double fm_polyn_1d_2_func(mat *X, mat *p, void *nothing)
{
    double x,res;
    
    nothing = NULL;
    x = mat_get(X,0,0);
    
    res = mat_get(p,0,0)             \
          + mat_get(p,1,0)*x         \
          + mat_get(p,2,0)*SQ(x);
    
    return res;
}

const fit_model fm_polyn_1d_2 =
{
    "y(x) = p0 + p1*x + p2*x^2",
    {fm_polyn_1d_2_func},
    &npar_3,
    1,
    "%e+%e*x+%e*x**2"
};

double fm_polyn_1d_3_func(mat *X, mat *p, void *nothing)
{
    double x,res;
    
    nothing = NULL;
    x = mat_get(X,0,0);
    
    res = mat_get(p,0,0)             \
          + mat_get(p,1,0)*x         \
          + mat_get(p,2,0)*SQ(x)     \
          + mat_get(p,3,0)*x*x*x;
    return res;
}

const fit_model fm_polyn_1d_3 =
{
    "y(x) = p0 + p1*x + p2*x^2 + p3*x^3",
    {fm_polyn_1d_3_func},
    &npar_4,
    1,
    "%e+%e*x+%e*x**2+%e*x**3"
};

double fm_polyn_1d_4_func(mat *X, mat *p, void *nothing)
{
    double x,res;
    
    nothing = NULL;
    x = mat_get(X,0,0);
    
    res = mat_get(p,0,0)             \
          + mat_get(p,1,0)*x         \
          + mat_get(p,2,0)*SQ(x)     \
          + mat_get(p,3,0)*x*x*x     \
          + mat_get(p,4,0)*SQ(SQ(x));
    return res;
}

const fit_model fm_polyn_1d_4 =
{
    "y(x) = p0 + p1*x + p2*x^2 + p3*x^3 + p4*x^4",
    {fm_polyn_1d_4_func},
    &npar_5,
    1,
    "%e+%e*x+%e*x**2+%e*x**3+%e*x**4"
};

const fit_model *fm_polyn_1d[FM_POLYN_1D_MAXDEG+1] = 
{
    &fm_polyn_1d_0,\
    &fm_polyn_1d_1,\
    &fm_polyn_1d_2,\
    &fm_polyn_1d_3,\
    &fm_polyn_1d_4
};

/** exponential decay **/
double fm_expdec_func(mat *x, mat *p, void *nothing)
{
    double res;
    
    nothing = NULL;
    res = mat_get(p,1,0)*exp(-mat_get(p,0,0)*mat_get(x,0,0));
    
    return res;
}

const fit_model fm_expdec = 
{
    "y(x) = p1*exp(-p0*x)",
    {fm_expdec_func},
    &npar_2,
    1,
    "exp(-%e*x)*%e"
};

/** hyperbolic cosine **/
double fm_cosh_func(mat *x, mat *p, void *nothing)
{
    double res;
    
    nothing = NULL;
    res = mat_get(p,1,0)*cosh(mat_get(p,0,0)*mat_get(x,0,0));
    
    return res;
}

const fit_model fm_cosh =
{
    "y(x) = p1*cosh(p0*x)",
    {fm_cosh_func},
    &npar_2,
    1,
    "cosh(%e*x)*%e"
};

/*                              2D models                                   */
/****************************************************************************/
/** 2D polynomial models **/
double fm_polyn_2d_00_func(mat *X, mat *p, void *nothing)
{
    mat *dummy;
    
    nothing = NULL;
    dummy   = X;
    
    return mat_get(p,0,0);
}

const fit_model fm_polyn_2d_00 =
{
    "z(x,y) = p0",
    {fm_polyn_2d_00_func},
    &npar_1,
    2,
    "%e"
};

double fm_polyn_2d_01_func(mat *X, mat *p, void *nothing)
{
    double res,x,y;
    
    nothing = NULL;
    x = mat_get(X,0,0);
    y = mat_get(X,1,0);
    
    res = mat_get(p,0,0)             \
          + mat_get(p,1,0)*y;
    
    return res;
}

const fit_model fm_polyn_2d_01 =
{
    "z(x,y) = p0+p1*y",
    {fm_polyn_2d_01_func},
    &npar_2,
    2,
    "%e+%e*y"
};

double fm_polyn_2d_02_func(mat *X, mat *p, void *nothing)
{
    double res,x,y;
    
    nothing = NULL;
    x = mat_get(X,0,0);
    y = mat_get(X,1,0);
    
    res = mat_get(p,0,0)             \
          + mat_get(p,1,0)*y         \
          + mat_get(p,2,0)*SQ(y);
    
    return res;
}

const fit_model fm_polyn_2d_02 =
{
    "z(x,y) = p0 + p1*y + p2*y^2",
    {fm_polyn_2d_02_func},
    &npar_3,
    2,
    "%e+%e*y+%e*y**2"
};

double fm_polyn_2d_10_func(mat *X, mat *p, void *nothing)
{
    double res,x,y;
    
    nothing = NULL;
    x = mat_get(X,0,0);
    y = mat_get(X,1,0);
    
    res = mat_get(p,0,0)             \
          + mat_get(p,1,0)*x;
    
    return res;
}

const fit_model fm_polyn_2d_10 =
{
    "z(x,y) = p0 + p1*x",
    {fm_polyn_2d_10_func},
    &npar_2,
    2,
    "%e+%e*x"
};

double fm_polyn_2d_11_func(mat *X, mat *p, void *nothing)
{
    double res,x,y;
    
    nothing = NULL;
    x = mat_get(X,0,0);
    y = mat_get(X,1,0);
    
    res = mat_get(p,0,0)             \
          + mat_get(p,1,0)*x         \
          + mat_get(p,2,0)*y         \
          + mat_get(p,3,0)*x*y;
    
    return res;
}

const fit_model fm_polyn_2d_11 =
{
    "z(x,y) = p0 + p1*x + p2*y + p3*x*y",
    {fm_polyn_2d_11_func},
    &npar_4,
    2,
    "%e+%e*x+%e*y+%e*x*y"
};

double fm_polyn_2d_12_func(mat *X, mat *p, void *nothing)
{
    double res,x,y;
    
    nothing = NULL;
    x = mat_get(X,0,0);
    y = mat_get(X,1,0);
    
    res = mat_get(p,0,0)             \
          + mat_get(p,1,0)*x         \
          + mat_get(p,2,0)*y         \
          + mat_get(p,3,0)*x*y       \
          + mat_get(p,4,0)*SQ(y);
    
    return res;
}

const fit_model fm_polyn_2d_12 =
{
    "z(x,y) = p0 + p1*x + p2*y + p3*x*y + p4*y^2",
    {fm_polyn_2d_12_func},
    &npar_5,
    2,
    "%e+%e*x+%e*y+%e*x*y+%e*y**2"
};

double fm_polyn_2d_20_func(mat *X, mat *p, void *nothing)
{
    double res,x,y;
    
    nothing = NULL;
    x = mat_get(X,0,0);
    y = mat_get(X,1,0);
    
    res = mat_get(p,0,0)             \
          + mat_get(p,1,0)*x         \
          + mat_get(p,2,0)*SQ(x);
    
    return res;
}

const fit_model fm_polyn_2d_20 =
{
    "z(x,y) = p0 + p1*x + p2*x^2",
    {fm_polyn_2d_20_func},
    &npar_3,
    2,
    "%e+%e*x+%e*x**2"
};

double fm_polyn_2d_21_func(mat *X, mat *p, void *nothing)
{
    double res,x,y;
    
    nothing = NULL; 
    x = mat_get(X,0,0);
    y = mat_get(X,1,0);
    
    res = mat_get(p,0,0)             \
          + mat_get(p,1,0)*x         \
          + mat_get(p,2,0)*y         \
          + mat_get(p,3,0)*x*y       \
          + mat_get(p,4,0)*SQ(x);
    
    return res;
}

const fit_model fm_polyn_2d_21 =
{
    "z(x,y) = p0 + p1*x + p2*y + p3*x*y + p4*x^2",
    {fm_polyn_2d_21_func},
    &npar_5,
    2,
    "%e+%e*x+%e*y+%e*x*y+%e*x**2"
};

double fm_polyn_2d_22_func(mat *X, mat *p, void *nothing)
{
    double res,x,y;
    
    nothing = NULL;
    x = mat_get(X,0,0);
    y = mat_get(X,1,0);
    
    res = mat_get(p,0,0)             \
          + mat_get(p,1,0)*x         \
          + mat_get(p,2,0)*y         \
          + mat_get(p,3,0)*x*y       \
          + mat_get(p,4,0)*SQ(x)     \
          + mat_get(p,5,0)*SQ(y);
    
    return res;
}

const fit_model fm_polyn_2d_22 =
{
    "z(x,y) = p0 + p1*x + p2*y + p3*x*y + p4*x^2 + p5*y^2",
    {fm_polyn_2d_22_func},
    &npar_6,
    2,
    "%e+%e*x+%e*y+%e*x*y+%e*x**2+%e*y**2"
};

const fit_model *fm_polyn_2d[FM_POLYN_2D_MAXDEG+1][FM_POLYN_2D_MAXDEG+1] =
{
    {&fm_polyn_2d_00, &fm_polyn_2d_01, &fm_polyn_2d_02},\
    {&fm_polyn_2d_10, &fm_polyn_2d_11, &fm_polyn_2d_12},\
    {&fm_polyn_2d_20, &fm_polyn_2d_21, &fm_polyn_2d_22}
};
