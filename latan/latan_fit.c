#include <latan/latan_fit.h>
#include <latan/latan_includes.h>
#include <latan/latan_math.h>
#include <latan/latan_minimizer.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*							fit model structure								*/
/****************************************************************************/
/** some useful constant npar_func **/
#define DEFINE_CST_NPAR_FUNC(n)\
size_t npar_##n(int stage, void *nothing)\
{\
	stage   = 0;\
	nothing = NULL;\
	return n;\
}

DEFINE_CST_NPAR_FUNC(1)
DEFINE_CST_NPAR_FUNC(2)
DEFINE_CST_NPAR_FUNC(3)
DEFINE_CST_NPAR_FUNC(4)
DEFINE_CST_NPAR_FUNC(5)
DEFINE_CST_NPAR_FUNC(6)
DEFINE_CST_NPAR_FUNC(7)
DEFINE_CST_NPAR_FUNC(8)
DEFINE_CST_NPAR_FUNC(9)
DEFINE_CST_NPAR_FUNC(10)

/** access **/
void fit_model_get_name(strbuf name, const fit_model *model)
{
	strcpy(name,model->name);
}

void fit_model_get_plot_fmt(strbuf plot_fmt, const fit_model *model)
{
	strcpy(plot_fmt,model->plot_fmt);
}

double fit_model_eval(const fit_model *model, mat *x, mat *p,\
					  const size_t stage, void *model_param)
{
	size_t i;
	double res;
	
	res = 0.0;
	for (i=0;i<=stage;i++)
	{
		res += model->func[i](x,p,model_param);
	}
	return res;
}

/*							fit data structure								*/
/****************************************************************************/
/** allocation **/
fit_data *fit_data_create(const size_t ndata, const size_t ndim)
{
	fit_data *d;
	size_t i;
	
	MALLOC_ERRVAL(d,fit_data *,1,NULL);
	d->x           = mat_create(ndata*ndim,1);
	d->x_varinv    = mat_create(ndata*ndim,ndata*ndim);
	mat_id(d->x_varinv);
	d->data        = mat_create(ndata,1);
	d->data_varinv = mat_create(ndata,ndata);
	mat_id(d->data_varinv);
	d->buf_chi2[0] = mat_create(ndata,1);
	d->buf_chi2[1] = mat_create(ndata,1);
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

void fit_data_destroy(fit_data *d)
{
	mat_destroy(d->x);
	mat_destroy(d->x_varinv);
	mat_destroy(d->data);
	mat_destroy(d->data_varinv);
	mat_destroy(d->buf_chi2[0]);
	mat_destroy(d->buf_chi2[1]);
	FREE(d->to_fit);
	FREE(d);
}

/** access **/
/*** chi2 value ***/
void fit_data_save_chi2pdof(fit_data *d, bool save)
{
	d->save_chi2pdof = save;
}

double fit_data_get_chi2pdof(fit_data *d)
{
	return d->chi2pdof;
}

/*** fit points ***/
void fit_data_set_x(fit_data *d, const size_t i, const size_t j,\
					const double x_ij)
{
	mat_set(d->x,i*d->ndim+j,0,x_ij);
}

double fit_data_get_x(const fit_data *d, const size_t i, const size_t j)
{
	return mat_get(d->x,i*d->ndim+j,0);
}

mat *fit_data_pt_x(const fit_data *d)
{
	return d->x;
}

latan_errno fit_data_set_x_var(fit_data *d, mat *var)
{
	latan_errno status;
	size_t i;
	double diag_i;
	strbuf warnmsg;
	
	status = LATAN_SUCCESS;
	d->have_x_var      = true;
	d->is_x_correlated = mat_is_square(var);
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

bool fit_data_have_x_var(fit_data *d)
{
	return d->have_x_var;
}

void fit_data_fit_all_points(fit_data *d, bool fit)
{
	size_t i;
	
	for (i=0;i<d->ndata;i++)
	{
		d->to_fit[i] = fit;
	}
}

void fit_data_fit_point(fit_data *d, size_t i, bool fit)
{
	if (i>=d->ndata)
	{
		LATAN_ERROR_VOID("index out of range",LATAN_EBADLEN);
	}
	
	d->to_fit[i] = fit;
}

void fit_data_fit_range(fit_data *d, size_t start, size_t end, bool fit)
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

bool fit_data_is_fit_point(const fit_data *d, size_t i)
{
	if (i>=d->ndata)
	{
		LATAN_ERROR_VAL("index out of range",LATAN_EBADLEN,false);
	}
	
	return d->to_fit[i];
}

size_t fit_data_fit_point_num(const fit_data *d)
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

/*** data ***/
void fit_data_set_data(fit_data *d, const size_t i, const double data_i)
{
	mat_set(d->data,i,0,data_i);
}

double fit_data_get_data(const fit_data *d, const size_t i)
{
	return mat_get(d->data,i,0);
}

size_t fit_data_get_ndata(const fit_data *d)
{
	return d->ndata;
}

mat *fit_data_pt_data(const fit_data *d)
{
	return d->data;
}

latan_errno fit_data_set_data_var(fit_data *d, mat *var)
{
	latan_errno status;
	size_t i;
	double diag_i;
	strbuf warnmsg;
	
	status = LATAN_SUCCESS;
	d->is_data_correlated = mat_is_square(var);
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

bool fit_data_is_correlated(const fit_data *d)
{
	return d->is_data_correlated;
}

/*** model ***/
latan_errno fit_data_set_model(fit_data *d, const fit_model *model)
{
	if (model->ndim != d->ndim)
	{
		LATAN_ERROR("fit model and fit data dimension mismatch",LATAN_EBADLEN);
	}
	
	d->model = model;
	
	return LATAN_SUCCESS;
}

const fit_model *fit_data_pt_model(fit_data *d)
{
	return d->model;
}

void fit_data_set_model_param(fit_data *d, void *model_param)
{
	d->model_param = model_param;
}

/* CRITICALLY CALLED FUNCTION
 * fit_data_model_eval is heavily called during one call of minimize function,
 * optimization is done using MAT_PT_SUBM macro
 */
double fit_data_model_eval(const fit_data *d, const size_t i,\
						   mat *p)
{
	mat x_i;
	
	x_i.mem_flag = CPU_LAST;
	
	MAT_PT_SUBM(&x_i,d->x,i*d->ndim,0,i*d->ndim+d->ndim-1,0);
	
	return fit_model_eval(d->model,&x_i,p,d->stage,d->model_param);
}

/*** stages ***/
void fit_data_set_stage(fit_data *d, const int stage)
{
	d->stage = stage;
}

int fit_data_get_stage(const fit_data *d)
{
	return d->stage;
}

/*** dof ***/
int fit_data_get_dof(const fit_data *d)
{
	return fit_data_fit_point_num(d) - d->model->npar(d->stage,d->model_param);
}

/*							chi2 function									*/
/****************************************************************************/
/* CRITICALLY CALLED FUNCTION
 * chi2 is heavily called during one call of minimize function,
 * optimization is done using vector/matrix views from GSL, matrix buffers
 * in fit_data structure, and ddot BLAS operation
 */
double chi2(mat *p, void *d)
{
	fit_data *dt;
	size_t ndata;
	size_t i;
	mat *X,*sigX;
	double res;
	
	dt         = (fit_data *)d;
	ndata      = fit_data_get_ndata(dt);

#ifdef _OPENMP
	#pragma omp critical
#endif
	{
		X          = dt->buf_chi2[0];
		sigX       = dt->buf_chi2[1];
		
		mat_zero(X);
		mat_zero(sigX);
		for (i=0;i<ndata;i++)
		{
			if (fit_data_is_fit_point(d,i))
			{
				mat_set(X,i,0,fit_data_model_eval(d,i,p)
						- fit_data_get_data(d,i));
			}
			else
			{
				mat_set(X,i,0,0.0);
			}
		}
		mat_mul_nn(sigX,dt->data_varinv,X);
		mat_dvmul(&res,sigX,X);
	}
	
	return res;
}

/*							fit functions									*/
/****************************************************************************/
latan_errno data_fit(mat *p, fit_data *d)
{
	latan_errno status;
	strbuf cor_status;
	double chi2_min;
	
	fit_data_is_correlated(d) ? \
	(strcpy(cor_status,"correlated")) : (strcpy(cor_status,"uncorrelated"));
	latan_printf(VERB,"fitting (%s) %u data points with model %s...\n",
				 cor_status,(unsigned int)fit_data_fit_point_num(d),\
				 d->model->name);
	status = minimize(p,&chi2_min,&chi2,d);
	if (d->save_chi2pdof)
	{
		d->chi2pdof = DRATIO(chi2_min,fit_data_get_dof(d));
	}
	
	return status;
}

latan_errno rs_data_fit(rs_sample *p, rs_sample *data, fit_data *d)
{
	latan_errno status;
	mat *datavar;
	size_t i;
	int verb_backup;
	double chi2pdof_backup;
	
	if (rs_sample_get_method(p) != GENERIC)
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
	LATAN_UPDATE_STATUS(status,data_fit(rs_sample_pt_cent_val(p),d));
	chi2pdof_backup = fit_data_get_chi2pdof(d);
	latan_printf(DEBUG,"central value chi^2/dof = %e\n",chi2pdof_backup);
	if (verb_backup != DEBUG)
	{
		latan_set_verb(QUIET);
	}
	for (i=0;i<rs_sample_get_nsample(data);i++)
	{
		mat_cp(rs_sample_pt_sample(p,i),\
			   rs_sample_pt_cent_val(p));
		mat_cp(fit_data_pt_data(d),rs_sample_pt_sample(data,i));
		LATAN_UPDATE_STATUS(status,data_fit(rs_sample_pt_sample(p,i),d));
		latan_printf(DEBUG,"sample %lu chi^2/dof = %e\n",(long unsigned)i,\
					 fit_data_get_chi2pdof(d));
	}
	d->chi2pdof = chi2pdof_backup;
	latan_set_verb(verb_backup);
	
	mat_destroy(datavar);
	
	return status;
}

latan_errno rs_x_data_fit(rs_sample *p, rs_sample *x, rs_sample *data,\
						  fit_data *d)
{
	latan_errno status;
	mat *datavar;
	size_t i;
	int verb_backup;
	double chi2pdof_backup;
	
	if (rs_sample_get_method(p) != GENERIC)
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
	LATAN_UPDATE_STATUS(status,data_fit(rs_sample_pt_cent_val(p),d));
	chi2pdof_backup = fit_data_get_chi2pdof(d);
	latan_printf(DEBUG,"central value chi^2/dof = %e\n",chi2pdof_backup);
	if (verb_backup != DEBUG)
	{
		latan_set_verb(QUIET);
	}
	for (i=0;i<rs_sample_get_nsample(data);i++)
	{
		mat_cp(rs_sample_pt_sample(p,i),\
			   rs_sample_pt_cent_val(p));
		mat_cp(fit_data_pt_data(d),rs_sample_pt_sample(data,i));
		mat_cp(fit_data_pt_x(d),rs_sample_pt_sample(x,i));
		LATAN_UPDATE_STATUS(status,data_fit(rs_sample_pt_sample(p,i),d));
		latan_printf(DEBUG,"sample %lu chi^2/dof = %e\n",(long unsigned)i,\
					 fit_data_get_chi2pdof(d));
	}
	d->chi2pdof = chi2pdof_backup;
	latan_set_verb(verb_backup);
	
	mat_destroy(datavar);
	
	return status;
}
