#include <latan/latan_models.h>
#include <latan/latan_includes.h>
#include <latan/latan_math.h>

/*** constant: y(x) = p0 ***/
double fm_const_func(const mat x, const mat func_param, void* nothing)
{
	mat dummy;
	
	nothing = NULL;
	dummy = x;
	
	return mat_get(func_param,0,0);
}

const fit_model fm_const =
{
	"y(x) = p0",
	&fm_const_func,
	1,
	1,
	"%e"
};

double fm_lin_func(const mat x, const mat func_param, void* nothing)
{
	double res;
	
	nothing = NULL;
	res = mat_get(func_param,0,0) + mat_get(func_param,1,0)*mat_get(x,0,0);
	
	return res;
}

const fit_model fm_lin =
{
	"y(x) = p0 + p1*x",
	&fm_lin_func,
	2,
	1,
	"%e+%e*x"
};

/*** exponential decay: y(x) = p1*exp(-p0*x) ***/
double fm_expdec_func(const mat x, const mat func_param, void* nothing)
{
	double res;
	
	nothing = NULL;
	res = mat_get(func_param,1,0)*exp(-mat_get(func_param,0,0)*mat_get(x,0,0));
	
	return res;
}

const fit_model fm_expdec = 
{
	"y(x) = p1*exp(-p0*x)",
	&fm_expdec_func,
	2,
	1,
	"exp(-%e*x)*%e"
};

/*** hyperbolic cosine: y(x) = p1*cosh(p0*x) ***/
double fm_cosh_func(const mat x, const mat func_param, void* nothing)
{
	double res;
	
	nothing = NULL;
	res = mat_get(func_param,1,0)*cosh(mat_get(func_param,0,0)*mat_get(x,0,0));
	
	return res;
}

const fit_model fm_cosh =
{
	"y(x) = p1*cosh(p0*x)",
	&fm_cosh_func,
	2,
	1,
	"cosh(%e*x)*%e"
};

double fm_polyn_2d_00_func(const mat X, const mat func_param, void* nothing)
{
	double res,x,y;
	
	nothing = NULL;
	x = mat_get(X,0,0);
	y = mat_get(X,0,1);
	
	res = mat_get(func_param,0,0);
	
	return res;
}

const fit_model fm_polyn_2d_00 =
{
	"z(x,y) = p_0",
	&fm_polyn_2d_00_func,
	1,
	2,
	"%e"
};

double fm_polyn_2d_01_func(const mat X, const mat func_param, void* nothing)
{
	double res,x,y;
	
	nothing = NULL;
	x = mat_get(X,0,0);
	y = mat_get(X,0,1);
	
	res = mat_get(func_param,0,0)             \
	+ mat_get(func_param,1,0)*y;
	
	return res;
}

const fit_model fm_polyn_2d_01 =
{
	"z(x,y) = p_0+p_1*y",
	&fm_polyn_2d_01_func,
	2,
	2,
	"%e+%e*y"
};

double fm_polyn_2d_02_func(const mat X, const mat func_param, void* nothing)
{
	double res,x,y;
	
	nothing = NULL;
	x = mat_get(X,0,0);
	y = mat_get(X,0,1);
	
	res = mat_get(func_param,0,0)             \
	      + mat_get(func_param,1,0)*y         \
	      + mat_get(func_param,2,0)*SQ(y);
	
	return res;
}

const fit_model fm_polyn_2d_02 =
{
	"z(x,y) = p_0+p_1*y+p_2*y^2",
	&fm_polyn_2d_02_func,
	3,
	2,
	"%e+%e*y+%e*y**2"
};

double fm_polyn_2d_10_func(const mat X, const mat func_param, void* nothing)
{
	double res,x,y;
	
	nothing = NULL;
	x = mat_get(X,0,0);
	y = mat_get(X,0,1);
	
	res = mat_get(func_param,0,0)             \
	      + mat_get(func_param,1,0)*x;
	
	return res;
}

const fit_model fm_polyn_2d_10 =
{
	"z(x,y) = p_0+p_1*x",
	&fm_polyn_2d_10_func,
	2,
	2,
	"%e+%e*x"
};

double fm_polyn_2d_11_func(const mat X, const mat func_param, void* nothing)
{
	double res,x,y;
	
	nothing = NULL;
	x = mat_get(X,0,0);
	y = mat_get(X,0,1);
	
	res = mat_get(func_param,0,0)             \
	      + mat_get(func_param,1,0)*x         \
	      + mat_get(func_param,2,0)*y         \
	      + mat_get(func_param,3,0)*x*y;
	
	return res;
}

const fit_model fm_polyn_2d_11 =
{
	"z(x,y) = p_0+p_1*x+p_2*y+p_3*x*y",
	&fm_polyn_2d_11_func,
	4,
	2,
	"%e+%e*x+%e*y+%e*x*y"
};

double fm_polyn_2d_12_func(const mat X, const mat func_param, void* nothing)
{
	double res,x,y;
	
	nothing = NULL;
	x = mat_get(X,0,0);
	y = mat_get(X,0,1);
	
	res = mat_get(func_param,0,0)             \
	      + mat_get(func_param,1,0)*x         \
	      + mat_get(func_param,2,0)*y         \
	      + mat_get(func_param,3,0)*x*y       \
	      + mat_get(func_param,4,0)*SQ(y);
	
	return res;
}

const fit_model fm_polyn_2d_12 =
{
	"z(x,y) = p_0+p_1*x+p_2*y+p_3*x*y+p_4*y^2",
	&fm_polyn_2d_12_func,
	5,
	2,
	"%e+%e*x+%e*y+%e*x*y+%e*y**2"
};

double fm_polyn_2d_20_func(const mat X, const mat func_param, void* nothing)
{
	double res,x,y;
	
	nothing = NULL;
	x = mat_get(X,0,0);
	y = mat_get(X,0,1);
	
	res = mat_get(func_param,0,0)             \
	      + mat_get(func_param,1,0)*x         \
	      + mat_get(func_param,2,0)*SQ(x);
	
	return res;
}

const fit_model fm_polyn_2d_20 =
{
	"z(x,y) = p_0+p_1*x+p_2*x^2",
	&fm_polyn_2d_20_func,
	3,
	2,
	"%e+%e*x+%e*x**2"
};

double fm_polyn_2d_21_func(const mat X, const mat func_param, void* nothing)
{
	double res,x,y;
	
	nothing = NULL;	
	x = mat_get(X,0,0);
	y = mat_get(X,0,1);
	
	res = mat_get(func_param,0,0)             \
	      + mat_get(func_param,1,0)*x         \
	      + mat_get(func_param,2,0)*y         \
	      + mat_get(func_param,3,0)*x*y       \
	      + mat_get(func_param,4,0)*SQ(x);
	
	return res;
}

const fit_model fm_polyn_2d_21 =
{
	"z(x,y) = p_0+p_1*x+p_2*y+p_3*x*y+p_4*x^2",
	&fm_polyn_2d_21_func,
	5,
	2,
	"%e+%e*x+%e*y+%e*x*y+%e*x**2"
};

double fm_polyn_2d_22_func(const mat X, const mat func_param, void* nothing)
{
	double res,x,y;
	
	nothing = NULL;
	x = mat_get(X,0,0);
	y = mat_get(X,0,1);
	
	res = mat_get(func_param,0,0)             \
	      + mat_get(func_param,1,0)*x         \
	      + mat_get(func_param,2,0)*y         \
	      + mat_get(func_param,3,0)*x*y       \
	      + mat_get(func_param,4,0)*SQ(x)     \
	      + mat_get(func_param,5,0)*SQ(y);
	
	return res;
}

const fit_model fm_polyn_2d_22 =
{
	"z(x,y) = p_0+p_1*x+p_2*y+p_3*x*y+p_4*x^2+p_5*y_2",
	&fm_polyn_2d_22_func,
	6,
	2,
	"%e+%e*x+%e*y+%e*x*y+%e*x**2+%e*y**2"
};

const fit_model* fm_polyn_2d[3][3] =
{
	{&fm_polyn_2d_00, &fm_polyn_2d_01, &fm_polyn_2d_02},\
	{&fm_polyn_2d_10, &fm_polyn_2d_11, &fm_polyn_2d_12},\
	{&fm_polyn_2d_20, &fm_polyn_2d_21, &fm_polyn_2d_22}
};
