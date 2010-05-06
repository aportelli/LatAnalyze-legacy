#include <latan/latan_minimizer.h>
#include <latan/latan_includes.h>
#include <latan/latan_min_minuit2.h>

latan_errno minimize(mat var, min_func* f, void* param)
{
	latan_errno status;
	
	switch (latan_get_minimize_lib())
	{
		case GSL:
			LATAN_ERROR("GSL support is not implemented yet",LATAN_FAILURE);
			break;
		case MINUIT:
			status = minimize_minuit2(var,f,param);
			break;
		default:
			LATAN_ERROR("minimizing library flag invalid",LATAN_EINVAL);
			break;
	}
	return status;
}