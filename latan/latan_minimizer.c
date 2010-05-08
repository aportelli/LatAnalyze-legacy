#include <latan/latan_minimizer.h>
#include <latan/latan_includes.h>
#include <latan/latan_min_minuit2.h>
#include <latan/latan_io.h>
#include <latan/latan_math.h>

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

/*							fit model structure								*/
/****************************************************************************/
/** access **/
void fit_model_get_name(stringbuf name, const fit_model* model)
{
	strcpy(name,model->name);
}

void fit_model_get_plot_fmt(stringbuf plot_fmt, const fit_model* model)
{
	strcpy(plot_fmt,model->plot_fmt);
}

double fit_model_eval(const fit_model* model, const double x,
					  const mat func_param, void* model_param)
{
	return model->func(x,func_param,model_param);
}

/** some useful models **/
/*** constant: y(x) = p0 ***/
double fm_const_func(const double x, const mat func_param, void* nothing)
{
	double dummy;
	
	nothing = NULL;
	dummy = x;
	
	return mat_get(func_param,0,0);
}

const fit_model fm_const =
{
	"y(x) = p0",
	&fm_const_func,
	1,
	"%e"
};

double fm_lin_func(const double x, const mat func_param, void* nothing)
{
	double res;
	
	nothing = NULL;
	res = mat_get(func_param,0,0) + mat_get(func_param,1,0)*x;
	
	return res;
}

const fit_model fm_lin =
{
	"y(x) = p0 + p1*x",
	&fm_lin_func,
	2,
	"%e+%e*x"
};

/*** exponential decay: y(x) = p0*exp(-p1*x) ***/
double fm_expdec_func(const double x, const mat func_param, void* nothing)
{
	double res;
	
	nothing = NULL;
	res = mat_get(func_param,0,0)*exp(-mat_get(func_param,1,0)*x);
	
	return res;
}

const fit_model fm_expdec = 
{
	"y(x) = p0*exp(-p1*x)",
	&fm_expdec_func,
	2,
	"%e*exp(-%e*x)"
};

/*** hyperbolic cosine: y(x) = p0*cosh(p1*x) ***/
double fm_cosh_func(const double x, const mat func_param, void* nothing)
{
	double res;
	
	nothing = NULL;
	res = mat_get(func_param,0,0)*cosh(mat_get(func_param,1,0)*x);
	
	return res;
}

const fit_model fm_cosh =
{
	"y(x) = p0*cosh(p1*x)",
	&fm_cosh_func,
	2,
	"%e*cosh(%e*x)"
};

/*							fit data structure								*/
/****************************************************************************/
/** allocation **/
fit_data fit_data_create(const size_t ndata)
{
	fit_data d;
	size_t i;
	
	MALLOC_ERRVAL(d,fit_data,1,NULL);
	d->x = mat_create(ndata,1);
	d->data = mat_create(ndata,1);
	d->var_inveigval = mat_create(ndata,1);
	d->var_eigvec = mat_create(ndata,ndata);
	MALLOC_ERRVAL(d->to_fit,bool*,ndata,NULL);
	for (i=0;i<ndata;i++)
	{
		d->to_fit[i] = true;
	}
				  
	d->model = NULL;
	d->model_param = NULL;
	d->stage = 0;
	d->ndata = ndata;
	d->save_chi2pdof = true;
	
	return d;
}

void fit_data_destroy(fit_data d)
{
	mat_destroy(d->x);
	mat_destroy(d->data);
	mat_destroy(d->var_inveigval);
	mat_destroy(d->var_eigvec);
	FREE(d->to_fit);
	FREE(d);
}

/** access **/
void fit_data_save_chi2pdof(fit_data d, bool save)
{
	d->save_chi2pdof = save;
}

double fit_data_get_chi2pdof(fit_data d)
{
	return d->chi2pdof;
}

void fit_data_set_x(fit_data d, const size_t i, const double x_i)
{
	mat_set(d->x,i,0,x_i);
}

double fit_data_get_x(const fit_data d, const size_t i)
{
	return mat_get(d->x,i,0);
}

mat fit_data_pt_x(const fit_data d)
{
	return d->x;
}

void fit_data_fit_all_points(fit_data d, bool fit)
{
	size_t i;
	
	for (i=0;i<d->ndata;i++)
	{
		d->to_fit[i] = fit;
	}
}

void fit_data_fit_point(fit_data d, size_t i, bool fit)
{
	if (i>=d->ndata)
	{
		LATAN_ERROR_VOID("index out of range",LATAN_EBADLEN);
	}
	
	d->to_fit[i] = fit;
}

bool fit_data_is_fit_point(fit_data d, size_t i)
{
	if (i>=d->ndata)
	{
		LATAN_ERROR_VAL("index out of range",LATAN_EBADLEN,false);
	}
	
	return d->to_fit[i];
}

size_t fit_data_fit_point_num(fit_data d)
{
	size_t nfitpt;
	size_t i;
	
	nfitpt = 0;
	
	for (i=0;i<d->ndata;i++)
	{
		if (d->to_fit[i])
		{
			nfitpt++;
		}
	}
	
	return nfitpt;
}

void fit_data_set_data(fit_data d, const size_t i, const double data_i)
{
	mat_set(d->data,i,0,data_i);
}

double fit_data_get_data(const fit_data d, const size_t i)
{
	return mat_get(d->data,i,0);
}

mat fit_data_pt_data(fit_data d)
{
	return d->data;
}

latan_errno fit_data_set_var(fit_data d, const mat var)
{
	latan_errno status;
	
	status = LATAN_SUCCESS;
	d->is_correlated = mat_issquare(var);
	
	if (d->is_correlated)
	{
		LATAN_ERROR("correlated fit is not implemented yet",LATAN_FAILURE);
	}
	else
	{
		mat_id(d->var_eigvec);
		mat_cst(d->var_inveigval,1.0);
		LATAN_UPDATE_STATUS(status,mat_eqdivp(d->var_inveigval,var));
	}
	
	return status;
}

void fit_data_set_model(fit_data d, const fit_model* model)
{
	d->model = model;
}

const fit_model* fit_data_pt_model(fit_data d)
{
	return d->model;
}

void fit_data_set_model_param(fit_data d, void* model_param)
{
	d->model_param = model_param;
}

double fit_data_model_eval(const fit_data d, const size_t i,\
						   const mat func_param)
{
	return fit_model_eval(d->model,mat_get(d->x,i,0),func_param,d->model_param);
}

void fit_data_set_stage(fit_data d, const int stage)
{
	d->stage = stage;
}

int fit_data_get_stage(const fit_data d)
{
	return d->stage;
}

int fit_data_get_dof(const fit_data d)
{
	return fit_data_fit_point_num(d) - d->model->npar;
}

bool fit_data_is_correlated(const fit_data d)
{
	return d->is_correlated;
}

/*							chi2 functions									*/
/****************************************************************************/
double chi2(const mat fit_param, void* d)
{
	fit_data dt;
	size_t ndata;
	size_t i;
	mat X;
	mat mres;
	double res;
	
	dt = (fit_data)d;
	ndata = nrow(fit_data_pt_data(d));
	
	X = mat_create(ndata,1);
	mres = mat_create(1,1);
	
	for (i=0;i<ndata;i++)
	{
		if (fit_data_is_fit_point(d,i))
		{
			mat_set(X,i,0,fit_data_model_eval(d,i,fit_param)
					- fit_data_get_data(d,i));
		}
		else
		{
			mat_set(X,i,0,0.0);
		}
	}
	mat_mul_nn(X,dt->var_eigvec,X);
	mat_eqdivp(X,dt->var_inveigval);
	mat_mul_tn(mres,X,X);
	res = mat_get(mres,0,0);
	
	mat_destroy(X);
	mat_destroy(mres);
	
	return res;
}

/*							fit functions									*/
/****************************************************************************/
latan_errno data_fit(mat fit_param, fit_data d)
{
	latan_errno status;
	stringbuf cor_status;
	double chi2_min;
	
	fit_data_is_correlated(d) ? \
	(strcpy(cor_status,"correlated")) : (strcpy(cor_status,"uncorrelated"));
	latan_printf(VERB,"fitting (%s) %u data points with model %s...\n",
				 cor_status,(unsigned int)fit_data_fit_point_num(d),\
				 d->model->name);
	status = minimize(fit_param,&chi2_min,&chi2,d);
	if (d->save_chi2pdof)
	{
		d->chi2pdof = DRATIO(chi2_min,fit_data_get_dof(d));
	}
	
	return status;
}
