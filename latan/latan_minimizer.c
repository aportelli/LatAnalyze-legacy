#include <latan/latan_minimizer.h>
#include <latan/latan_includes.h>
#include <latan/latan_min_minuit2.h>
#include <latan/latan_io.h>
#include <latan/latan_math.h>

/*						minimizer options									*/
/****************************************************************************/

#ifdef HAVE_MINUIT2
#define DEF_LIB MINUIT
#define DEF_ALG MIN_MIGRAD
#else
#define DEF_LIB GSL
#define DEF_ALG GSL_GRAD
#endif

#ifndef DEF_MAX_ITERATION
#define DEF_MAX_ITERATION 200u
#endif

typedef struct
{
	minlib_no lib;
	minalg_no alg;
	unsigned int max_iteration;
} minimizer_env;

static minimizer_env env = 
{
	DEF_LIB,\
	DEF_ALG,\
	DEF_MAX_ITERATION
};

minlib_no minimizer_get_lib(void)
{
	return env.lib;
}

minalg_no minimizer_get_alg(void)
{
	return env.alg;
}

latan_errno minimizer_set_alg(minalg_no alg)
{
	if (alg <= GSL_SIMPLEX_NM)
	{
		env.lib = GSL;
		env.alg = alg;
	}
	else if (alg <= MIN_SIMPLEX)
	{
		env.lib = MINUIT;
		env.alg = alg;
	}
	else
	{
		LATAN_ERROR("minimization algorithm flag invalid",LATAN_EINVAL);
	}
	
	return LATAN_SUCCESS;
}

latan_errno minimizer_get_alg_name(stringbuf name)
{
	switch (env.alg)
	{
		case GSL_GRAD_FR:
			strcpy(name,"Fletcher-Reeves conjugate gradient (GSL)");
			break;
		case GSL_GRAD_PR:
			strcpy(name,"Polak-Ribiere conjugate gradient (GSL)");
			break;
		case GSL_VEC_BFGS:
			strcpy(name,\
				   "improved Broyden-Fletcher-Goldfarb-Shanno vector (GSL)");
			break;
		case GSL_SIMPLEX_NM:
			strcpy(name,"improved Nelder-Mead simplex (GSL)");
			break;
		case MIN_MIGRAD:
			strcpy(name,"variable metric (MINUIT)");
			break;
		case MIN_SIMPLEX:
			strcpy(name,"simplex (MINUIT)");
			break;
		default:
			LATAN_ERROR("minimization algorithm flag invalid",LATAN_EINVAL);
			break;
	}
	
	return LATAN_SUCCESS;
}

unsigned int minimizer_get_max_iteration(void)
{
	return env.max_iteration;
}

void minimizer_set_max_iteration(unsigned int max_iteration)
{
	env.max_iteration = max_iteration;
}

/*							the minimizer									*/
/****************************************************************************/
latan_errno minimize(mat var, double* f_min, min_func* f, void* param)
{
	latan_errno status;
	stringbuf name;
	
	minimizer_get_alg_name(name);
	latan_printf(VERB,"minimizing using %s algorithm...\n",name);
	switch (minimizer_get_lib())
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

double fit_model_eval(const fit_model* model, const mat x,
					  const mat func_param, void* model_param)
{
	return model->func(x,func_param,model_param);
}

/*							fit data structure								*/
/****************************************************************************/
/** allocation **/
fit_data fit_data_create(const size_t ndata, const size_t ndim)
{
	fit_data d;
	size_t i;
	
	MALLOC_ERRVAL(d,fit_data,1,NULL);
	d->x           = mat_create(ndata,ndim);
	d->x_varinv    = mat_create(ndata*ndim,ndata*ndim);
	mat_id(d->x_varinv);
	d->data        = mat_create(ndata,1);
	d->data_varinv = mat_create(ndata,ndata);
	mat_id(d->data_varinv);
	MALLOC_ERRVAL(d->to_fit,bool*,ndata,NULL);
	for (i=0;i<ndata;i++)
	{
		d->to_fit[i] = false;
	}
	d->model              = NULL;
	d->model_param        = NULL;
	d->stage              = 0;
	d->ndata              = ndata;
	d->ndim               = ndim;
	d->is_data_correlated = false;
	d->is_x_correlated    = false;
	d->have_x_var         = false;
	d->save_chi2pdof      = true;
	
	return d;
}

void fit_data_destroy(fit_data d)
{
	mat_destroy(d->x);
	mat_destroy(d->x_varinv);
	mat_destroy(d->data);
	mat_destroy(d->data_varinv);
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

void fit_data_set_x(fit_data d, const size_t i, const size_t j,\
					const double x_ij)
{
	mat_set(d->x,i,j,x_ij);
}

double fit_data_get_x(const fit_data d, const size_t i, const size_t j)
{
	return mat_get(d->x,i,j);
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

void fit_data_fit_range(fit_data d, size_t start, size_t end, bool fit)
{
	size_t i;
	
	if ((start >= d->ndata)||(end >= d->ndata))
	{
		LATAN_ERROR_VOID("index out of range",LATAN_EBADLEN);
	}
	if (start > end)
	{
		LATAN_WARNING("trying to modify an empty fit range",LATAN_EBADLEN);
	}
	
	for (i=start;i<=end;i++)
	{
		d->to_fit[i] = fit;
	}
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

latan_errno fit_data_set_data_var(fit_data d, const mat var)
{
	latan_errno status;
	size_t i;
	double diag_i;
	stringbuf warnmsg;
	
	status = LATAN_SUCCESS;
	d->is_data_correlated = mat_issquare(var);
	if (d->is_data_correlated)
	{
		LATAN_UPDATE_STATUS(status,mat_inv(d->data_varinv,var));
	}
	else
	{
		mat_zero(d->data_varinv);
		for (i=0;i<nrow(var);i++)
		{
			if (mat_get(var,i,0) == 0)
			{
				sprintf(warnmsg,"singular point %lu eliminated from fit",\
						(long unsigned)i);
				LATAN_WARNING(warnmsg,LATAN_EDOM);
				diag_i = 0.0;
			}
			else
			{
				diag_i = 1.0/mat_get(var,i,0);
			}
			mat_set(d->data_varinv,i,i,diag_i);
		}
	}
	
	return status;
}

latan_errno fit_data_set_x_var(fit_data d, const mat var)
{
	latan_errno status;
	size_t i;
	double diag_i;
	stringbuf warnmsg;
	
	status = LATAN_SUCCESS;
	d->have_x_var      = true;
	d->is_x_correlated = mat_issquare(var);
	if (d->is_data_correlated)
	{
		LATAN_UPDATE_STATUS(status,mat_inv(d->x_varinv,var));
	}
	else
	{
		mat_zero(d->data_varinv);
		for (i=0;i<nrow(var);i++)
		{
			if (mat_get(var,i,0) == 0)
			{
				sprintf(warnmsg,"singular point %lu eliminated from fit",\
						(long unsigned)i);
				LATAN_WARNING(warnmsg,LATAN_EDOM);
				diag_i = 0.0;
			}
			else
			{
				diag_i = 1.0/mat_get(var,i,0);
			}
			mat_set(d->x_varinv,i,i,diag_i);
		}
	}
	
	return status;
}

latan_errno fit_data_set_model(fit_data d, const fit_model* model)
{
	if (model->ndim != d->ndim)
	{
		LATAN_ERROR("fit model and fit data dimension mismatch",LATAN_EBADLEN);
	}
	
	d->model = model;
	
	return LATAN_SUCCESS;
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
	mat x_i, x_i_t;
	
	x_i   = mat_create(d->ndim,1);
	x_i_t = mat_create(1,d->ndim);
	
	mat_get_subm(x_i,d->x,i,0,i,d->ndim-1);
	mat_transpose(x_i_t,x_i);
	
	return fit_model_eval(d->model,x_i,func_param,d->model_param);
	
	mat_destroy(x_i);
	mat_destroy(x_i_t);
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
	return d->is_data_correlated;
}

/*							chi2 functions									*/
/****************************************************************************/
double chi2(const mat fit_param, void* d)
{
	fit_data dt;
	size_t ndata;
	size_t i;
	mat X, sigX;
	mat mres;
	double res;
	
	dt = (fit_data)d;
	ndata = nrow(fit_data_pt_data(d));
	
	X    = mat_create(ndata,1);
	sigX = mat_create(ndata,1);
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
	mat_mul_nn(sigX,dt->data_varinv,X);
	mat_mul_tn(mres,X,sigX);
	res = mat_get(mres,0,0);
	
	mat_destroy(X);
	mat_destroy(sigX);
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

latan_errno rs_data_fit(rs_sample fit_param, rs_sample data, fit_data d)
{
	latan_errno status;
	mat datavar;
	size_t i;
	int verb_backup;
	double chi2pdof_backup;
	
	if (rs_sample_get_method(fit_param) != GENERIC)
	{
		LATAN_WARNING("resampled sample fit result do not have GENERIC method",\
					  LATAN_EINVAL);
	}
	
	verb_backup = latan_get_verb();
	status = LATAN_SUCCESS;
	
	datavar = mat_create(nrow(rs_sample_pt_cent_val(data)),1);
	
	LATAN_UPDATE_STATUS(status,rs_sample_varp(datavar,data));
	LATAN_UPDATE_STATUS(status,fit_data_set_data_var(d,datavar));
	mat_cp(fit_data_pt_data(d),rs_sample_pt_cent_val(data));
	LATAN_UPDATE_STATUS(status,data_fit(rs_sample_pt_cent_val(fit_param),d));
	chi2pdof_backup = fit_data_get_chi2pdof(d);
	latan_printf(DEBUG,"central value chi^2/dof = %e\n",chi2pdof_backup);
	if (verb_backup != DEBUG)
	{
		latan_set_verb(QUIET);
	}
	for (i=0;i<rs_sample_get_nsample(data);i++)
	{
		mat_cp(rs_sample_pt_sample(fit_param,i),\
			   rs_sample_pt_cent_val(fit_param));
		mat_cp(fit_data_pt_data(d),rs_sample_pt_sample(data,i));
		LATAN_UPDATE_STATUS(status,data_fit(rs_sample_pt_sample(fit_param,i),d));
		latan_printf(DEBUG,"sample %lu chi^2/dof = %e\n",(long unsigned)i,\
					 fit_data_get_chi2pdof(d));
	}
	d->chi2pdof = chi2pdof_backup;
	latan_set_verb(verb_backup);
	
	mat_destroy(datavar);
	
	return status;
}

latan_errno rs_x_data_fit(rs_sample fit_param, rs_sample x, rs_sample data,\
						  fit_data d)
{
	latan_errno status;
	mat datavar;
	size_t i;
	int verb_backup;
	double chi2pdof_backup;
	
	if (rs_sample_get_method(fit_param) != GENERIC)
	{
		LATAN_WARNING("resampled sample fit result do not have GENERIC method",\
					  LATAN_EINVAL);
	}
	
	verb_backup = latan_get_verb();
	status = LATAN_SUCCESS;
	
	datavar = mat_create(nrow(rs_sample_pt_cent_val(data)),1);
	
	LATAN_UPDATE_STATUS(status,rs_sample_varp(datavar,data));
	LATAN_UPDATE_STATUS(status,fit_data_set_data_var(d,datavar));
	mat_cp(fit_data_pt_data(d),rs_sample_pt_cent_val(data));
	mat_cp(fit_data_pt_x(d),rs_sample_pt_cent_val(x));
	LATAN_UPDATE_STATUS(status,data_fit(rs_sample_pt_cent_val(fit_param),d));
	chi2pdof_backup = fit_data_get_chi2pdof(d);
	latan_printf(DEBUG,"central value chi^2/dof = %e\n",chi2pdof_backup);
	if (verb_backup != DEBUG)
	{
		latan_set_verb(QUIET);
	}
	for (i=0;i<rs_sample_get_nsample(data);i++)
	{
		mat_cp(rs_sample_pt_sample(fit_param,i),\
			   rs_sample_pt_cent_val(fit_param));
		mat_cp(fit_data_pt_data(d),rs_sample_pt_sample(data,i));
		mat_cp(fit_data_pt_x(d),rs_sample_pt_sample(x,i));
		LATAN_UPDATE_STATUS(status,data_fit(rs_sample_pt_sample(fit_param,i),d));
		latan_printf(DEBUG,"sample %lu chi^2/dof = %e\n",(long unsigned)i,\
					 fit_data_get_chi2pdof(d));
	}
	d->chi2pdof = chi2pdof_backup;
	latan_set_verb(verb_backup);
	
	mat_destroy(datavar);
	
	return status;
}
