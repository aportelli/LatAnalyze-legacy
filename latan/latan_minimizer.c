#include <latan/latan_minimizer.h>
#include <latan/latan_includes.h>
#include <latan/latan_min_minuit2.h>
#include <latan/latan_io.h>
#include <latan/latan_math.h>
#include <latan/latan_mass.h>

#ifndef PLAT_TOL
#define PLAT_TOL 0.25
#endif
#ifndef NSIGMA
#define NSIGMA 1.0
#endif

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

double fit_model_eval(const fit_model* model, const mat x,
					  const mat func_param, void* model_param)
{
	return model->func(x,func_param,model_param);
}

/** some useful models **/
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

/*							fit data structure								*/
/****************************************************************************/
/** allocation **/
fit_data fit_data_create(const size_t ndata, const size_t ndim)
{
	fit_data d;
	size_t i;
	
	MALLOC_ERRVAL(d,fit_data,1,NULL);
	d->x = mat_create(ndata,ndim);
	d->data = mat_create(ndata,1);
	d->var_inveigval = mat_create(ndata,1);
	d->var_eigvec = mat_create(ndata,ndata);
	MALLOC_ERRVAL(d->to_fit,bool*,ndata,NULL);
	for (i=0;i<ndata;i++)
	{
		d->to_fit[i] = false;
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

void fit_data_set_x(fit_data d, const size_t i, const size_t j,\
					const double x_i)
{
	mat_set(d->x,i,j,x_i);
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

latan_errno fit_data_set_var(fit_data d, const mat var)
{
	latan_errno status;
	size_t i;
	double dieg_i;
	
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
		for (i=0;i<nrow(var);i++)
		{
			dieg_i = (mat_get(var,i,0) == 0) ? (0.0) : (1.0/mat_get(var,i,0));
			mat_set(d->var_inveigval,i,0,dieg_i);
		}
	}
	
	return status;
}

latan_errno fit_data_set_model(fit_data d, const fit_model* model)
{
	if (model->ndim != ncol(d->x))
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
	mat x_i;
	
	x_i = mat_create(1,ncol(d->x));
	
	mat_cp_subm(x_i,d->x,i,0,i,ncol(x_i)-1);
	return fit_model_eval(d->model,x_i,func_param,d->model_param);
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

/** tuning **/
latan_errno fit_data_mass_fit_tune(fit_data d, mat fit_init, const mat prop,\
								   const mat em, const mat sigem,			\
								   const int parity)
{
	plat* em_plat;
	size_t nplat,nt,ntmax;
	size_t p,t;
	double shift,mem,pref,logslope;
	stringbuf ranges,buf;
	const fit_model* model;
	
	nt    = nrow(em) + 2;
	ntmax = (parity == EVEN) ? nt/2 : nt-2;
	
	/* setting fit model */
	switch (parity)
	{
		case EVEN:
			model = &fm_expdec;
			break;
		case ODD:
			model = &fm_cosh;
			break;
		default:
			LATAN_ERROR("wrong parity flag",LATAN_EINVAL);
			break;
	}
	fit_data_set_model(d,model);
	
	/* setting datas */
	switch (parity)
	{
		case EVEN:
			shift = 0.0;
			break;
		case ODD:
			shift = -DRATIO(nt,2.0);
			break;
		default:
			LATAN_ERROR("wrong parity flag",LATAN_EINVAL);
			break;
	}
	for (t=0;t<nt;t++)
	{
		fit_data_set_x(d,t,0,(double)(t)+shift);
	}
	
	/* searching mass plateaux */
	latan_printf(VERB,"searching mass plateaux in range [1,%lu]...\n",
				 (long unsigned)ntmax);
	em_plat = search_plat(&nplat,em,sigem,ntmax-1,NSIGMA,PLAT_TOL);
	
	/* setting points to fit */
	fit_data_fit_all_points(d,false);
	strcpy(ranges,"");
	for (p=0;p<nplat;p++)
	{
		fit_data_fit_range(d,em_plat[p].start+1,em_plat[p].end+1,true);
		sprintf(buf,"[%u,%u] ",(unsigned int)em_plat[p].start+1,\
				(unsigned int)em_plat[p].end+1);
		strcat(ranges,buf);
	}
	latan_printf(VERB,"fit ranges set to : %s\n",ranges);
	
	/* setting initial fit parameters */
	latan_printf(VERB,"searching initial parameter values...\n");
	mem = 0.0;
	for (p=0;p<nplat;p++)
	{
		mem += em_plat[p].mean;
	}
	mem /= (double)(nplat);
	switch (parity)
	{
		case EVEN:
			logslope = log(mat_get(prop,nt/8,0))*(1.0+nt/8) \
			           - log(mat_get(prop,nt/8+1,0))*nt/8;
			pref = exp(logslope);
			break;
		case ODD:
			pref = mat_get(prop,nt/2,0);
			break;
		default:
			LATAN_ERROR("wrong parity flag",LATAN_EINVAL);
			break;
	}
	latan_printf(VERB,"prefactor = %e mass = %e\n",pref,mem);
	mat_set(fit_init,1,0,pref);
	mat_set(fit_init,0,0,mem);
	
	FREE(em_plat);
	
	return LATAN_SUCCESS;
}

/*							chi2 functions									*/
/****************************************************************************/
double chi2(const mat fit_param, void* d)
{
	fit_data dt;
	size_t ndata;
	size_t i;
	mat X,Xonsig;
	mat mres;
	double res;
	
	dt = (fit_data)d;
	ndata = nrow(fit_data_pt_data(d));
	
	X = mat_create(ndata,1);
	Xonsig = mat_create(ndata,1);
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
	mat_mulp(Xonsig,X,dt->var_inveigval);
	mat_mul_tn(mres,Xonsig,X);
	res = mat_get(mres,0,0);
	
	mat_destroy(X);
	mat_destroy(Xonsig);
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

latan_errno rs_sample_fit(rs_sample fit_param, rs_sample data, fit_data d)
{
	latan_errno status;
	mat datavar;
	size_t i;
	bool save_chi2pdof_backup;
	int verb_backup;
	
	if (rs_sample_get_method(fit_param) != GENERIC)
	{
		LATAN_WARNING("resampled sample fit result do not have GENERIC method",\
					  LATAN_EINVAL);
	}
	
	save_chi2pdof_backup = d->save_chi2pdof;
	verb_backup = latan_get_verb();
	status = LATAN_SUCCESS;
	
	datavar = mat_create_from_dim(rs_sample_pt_cent_val(data));
	
	LATAN_UPDATE_STATUS(status,rs_sample_varp(datavar,data));
	LATAN_UPDATE_STATUS(status,fit_data_set_var(d,datavar));
	mat_cp(fit_data_pt_data(d),rs_sample_pt_cent_val(data));
	LATAN_UPDATE_STATUS(status,data_fit(rs_sample_pt_cent_val(fit_param),d));
	fit_data_save_chi2pdof(d,false);
	if (verb_backup != DEBUG)
	{
		latan_set_verb(QUIET);
	}
	for (i=0;i<rs_sample_get_nsample(data);i++)
	{
		mat_cp(fit_data_pt_data(d),rs_sample_pt_sample(data,i));
		LATAN_UPDATE_STATUS(status,data_fit(rs_sample_pt_sample(fit_param,i),d));
	}
	fit_data_save_chi2pdof(d,save_chi2pdof_backup);
	latan_set_verb(verb_backup);
	
	mat_destroy(datavar);
	
	return status;
}
