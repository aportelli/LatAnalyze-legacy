#include <latan/latan_minimizer.h>
#include <latan/latan_includes.h>
#include <latan/latan_min_minuit2.h>
#include <latan/latan_io.h>

/*							the minimizer									*/
/****************************************************************************/
latan_errno minimize(mat var, double* f_min, min_func* f, void* param)
{
	latan_errno status;
	
	switch (latan_get_minimize_lib())
	{
		case GSL:
			LATAN_ERROR("GSL support is not implemented yet",LATAN_FAILURE);
			break;
		case MINUIT:
			status = minimize_minuit2(var,f_min,f,param);
			break;
		default:
			LATAN_ERROR("minimizing library flag invalid",LATAN_EINVAL);
			break;
	}
	return status;
}

/*							fit data structure								*/
/****************************************************************************/
/** allocation **/
fit_data fit_data_create(mat x, mat data, mat var, model_func* model,\
						 void* model_param)
{
	fit_data d;
	size_t ndata;
	
	MALLOC_ERRVAL(d,fit_data,1,NULL);
	
	d->is_correlated = mat_issquare(var);
	if (nrow(var) != nrow(data))
	{
		LATAN_ERROR_NULL("correlation matrix have not the same number of row than data vector",LATAN_EBADLEN);
	}
	if (nrow(var) != nrow(x))
	{
		LATAN_ERROR_NULL("correlation matrix have not the same number of row than point vector",LATAN_EBADLEN);
	}
	ndata = nrow(data);
	d->x = mat_create(ndata,1);
	d->data = mat_create(ndata,1);
	d->var_inveigval = mat_create(ndata,1);
	d->var_eigvec = mat_create(ndata,ndata);
	mat_cp(d->x,x);
	mat_cp(d->data,data);
	if (d->is_correlated)
	{
		LATAN_ERROR_NULL("correlated fit is not implemented yet",LATAN_FAILURE);
	}
	else
	{
		mat_id(d->var_eigvec);
		mat_cst(d->var_inveigval,1.0);
		mat_eqdivp(d->var_inveigval,var);
	}
	d->model = model;
	d->model_param = model_param;
	d->stage = 0;
	
	return d;
}

void fit_data_destroy(fit_data d)
{
	mat_destroy(d->x);
	mat_destroy(d->data);
	mat_destroy(d->var_inveigval);
	mat_destroy(d->var_eigvec);
	FREE(d);
}

/*							chi2 functions									*/
/****************************************************************************/
double chi2(const mat var, void* d)
{
	fit_data dt;
	size_t ndata;
	size_t i;
	mat X;
	mat res;
	
	dt = (fit_data)d;
	ndata = nrow(dt->data);
	
	X = mat_create(ndata,1);
	res = mat_create(1,1);
	
	for (i=0;i<ndata;i++)
	{
		mat_set(X,i,0,dt->model(mat_get(dt->x,i,0),var,dt->model_param));
	}
	mat_eqsub(X,dt->data);
	mat_mul_nn(X,dt->var_eigvec,X);
	mat_eqdivp(X,dt->var_inveigval);
	mat_mul_tn(res,X,X);
	
	return mat_get(res,0,0);
}
