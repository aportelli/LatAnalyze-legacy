#include <latan/latan_min_minuit2.h>
#include <latan/latan_includes.h>

#ifdef HAVE_MINUIT2

#include <iostream>
#include <vector>
#include <Minuit2/FCNBase.h>
#include <Minuit2/CombinedMinimizer.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnPrint.h>

#ifndef INIT_RERROR
#define INIT_RERROR 1.0e-2
#endif
#ifndef STRATEGY
#define STRATEGY 2
#endif
#ifndef FIT_TOL
#define FIT_TOL 1.0e-6
#endif

// WARNING : std namespace already contain a stringbuf type
using namespace ROOT;
using namespace Minuit2;

class Minuit2MinFunc: public FCNBase
{
public:
	Minuit2MinFunc(min_func *init_f, void *init_param);
	~Minuit2MinFunc(void);
	
	virtual double operator()(const std::vector<double>& v_var) const;
	virtual double Up(void) const;
	
private:
	min_func *f;
	void *param;
};

Minuit2MinFunc::Minuit2MinFunc(min_func *init_f, void *init_param)
{
	f = init_f;
	param = init_param;
}

Minuit2MinFunc::~Minuit2MinFunc(void)
{
}

double Minuit2MinFunc::operator()(const std::vector<double>& v_x) const
{
	mat x;
	size_t x_size;
	size_t i;
	double res;
	
	x_size = v_x.size();
	
	x = mat_create(x_size,1);
	
	for (i=0;i<x_size;i++)
	{
		mat_set(x,i,0,v_x[i]);
	}
	res = f(x,param);
	
	mat_destroy(x);
	
	return res;
}

double Minuit2MinFunc::Up(void) const
{
	return 1.0;
}
#endif

latan_errno minimize_minuit2(mat x, double *f_min, min_func *f, void *param)
{
#ifdef HAVE_MINUIT2
	latan_errno status;
	size_t x_size;
	size_t i;
	double init_x_i,x_i;
	unsigned int max_iteration;
	std::vector<double> v_init_x;
	std::vector<double> v_init_err;
	Minuit2MinFunc minuit2_f(f,param);
	VariableMetricMinimizer minimizer_migrad;
	SimplexMinimizer minimizer_simplex;
	ModularFunctionMinimizer *minimizer;
	
	status = LATAN_SUCCESS;
	x_size = nrow(x);
	max_iteration = minimizer_get_max_iteration();
	
	for (i=0;i<x_size;i++)
	{
		init_x_i = mat_get(x,i,0);
		v_init_x.push_back(init_x_i);
		v_init_err.push_back(init_x_i*INIT_RERROR);
	}
	switch (minimizer_get_alg())
	{
		case MIN_MIGRAD:
			minimizer = &minimizer_migrad;
			break;
		case MIN_SIMPLEX:
			minimizer = &minimizer_simplex;
			break;
		default:
			LATAN_ERROR("invalid MINUIT minimization algorithm flag",\
						LATAN_EINVAL);
	}
	FunctionMinimum minuit2_min = minimizer->Minimize(minuit2_f,v_init_x, \
													  v_init_err,STRATEGY,  \
													  max_iteration,FIT_TOL);
	if (!minuit2_min.IsValid())
	{
		LATAN_WARNING("MINUIT library reported that minimization result is not valid",\
					  LATAN_FAILURE);
		status = LATAN_FAILURE;
	}
	for (i=0;i<x_size;i++)
	{
		x_i = minuit2_min.UserParameters().Parameter((unsigned int)i).Value();
		mat_set(x,i,0,x_i);
	}
	*f_min = minuit2_min.Fval();

	latan_printf(DEBUG,"MINUIT minimizer call :\n");
	if (latan_get_verb() == DEBUG)
	{
		std::cout << "--------------------------------------------------------";
		std::cout << minuit2_min;
		std::cout << "--------------------------------------------------------";
		std::cout << std::endl;
	}
	
	return status;
#else
	x = NULL;
	f_min = NULL;
	f = NULL;
	param = NULL;
	LATAN_ERROR("MINUIT library support was not compiled",LATAN_FAILURE);
#endif
}