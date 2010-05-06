#include <latan/latan_min_minuit2.h>
#include <latan/latan_includes.h>

#ifdef HAVE_MINUIT2

#include <iostream>
#include <vector>
#include <Minuit2/FCNBase.h>
#include <Minuit2/CombinedMinimizer.h>
#include <Minuit2/FunctionMinimum.h>

#define INIT_RERROR 0.25

// WARNING : std namespace already contain a stringbuf type
using namespace ROOT;
using namespace Minuit2;

class Minuit2MinFunc: public FCNBase
{
public:
	Minuit2MinFunc(min_func* init_f, void* init_param);
	~Minuit2MinFunc(void);
	
	virtual double operator()(const std::vector<double>& v_var) const;
	virtual double Up(void) const;
	
private:
	min_func* f;
	void* param;
};

Minuit2MinFunc::Minuit2MinFunc(min_func* init_f, void* init_param)
{
	f = init_f;
	param = init_param;
}

Minuit2MinFunc::~Minuit2MinFunc(void)
{
}

double Minuit2MinFunc::operator()(const std::vector<double>& v_var) const
{
	mat var;
	size_t var_size;
	size_t i;
	double res;
	
	var_size = v_var.size();
	
	var = mat_create(var_size,1);
	
	for (i=0;i<var_size;i++)
	{
		mat_set(var,i,0,v_var[i]);
	}
	res = f(var,param);
	
	mat_destroy(var);
	
	return res;
}

double Minuit2MinFunc::Up(void) const
{
	return 1.0;
}
#endif

latan_errno minimize_minuit2(mat var, min_func* f, void* param)
{
#ifdef HAVE_MINUIT2
	latan_errno status;
	size_t var_size;
	size_t i;
	double init_var_i,var_i;
	std::vector<double> v_init_var;
	std::vector<double> v_init_err;
	Minuit2MinFunc minuit2_f(f,param);
	CombinedMinimizer minimizer;
	
	status = LATAN_SUCCESS;
	var_size = nrow(var);
	
	for (i=0;i<var_size;i++)
	{
		init_var_i = mat_get(var,i,0);
		v_init_var.push_back(init_var_i);
		v_init_err.push_back(init_var_i*INIT_RERROR);
	}
	FunctionMinimum minuit2_min = minimizer.Minimize(minuit2_f,v_init_var,\
													 v_init_err);
	for (i=0;i<var_size;i++)
	{
		var_i = minuit2_min.UserParameters().Parameter((unsigned int)i).Value();
		mat_set(var,i,0,var_i);
	}

	return status;
#else
	mat dummy;
	
	dummy = res;
	dummy = init_var;
	f = NULL;
	param = NULL;
	LATAN_ERROR("MINUIT library support was not compiled",LATAN_FAILURE);
#endif
}