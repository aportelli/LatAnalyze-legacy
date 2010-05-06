#include <latan/latan_mass.h>
#include <latan/latan_includes.h>
#include <latan/latan_math.h>

latan_errno effmass_PCAC(mat res, const mat mprop_AP, const mat mprop_PP)
{
	latan_errno status;
	size_t t;
	
	if (!mat_issamedim(mprop_AP,mprop_PP))
	{
		LATAN_ERROR("AP and PP propagators have different dimensions",\
					LATAN_EBADLEN);
	}
	
	status = finite_diff(res,mprop_AP);
	for (t=0;t<nrow(res);t++)
	{
		mat_set(res,t,0,mat_get(res,t,0)/(2.0*mat_get(mprop_PP,t+1,0)));
	}
	mat_eqabs(res);
	
	return status;
}