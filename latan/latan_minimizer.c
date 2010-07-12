#include <latan/latan_minimizer.h>
#include <latan/latan_includes.h>
#include <latan/latan_min_gsl.h>
#include <latan/latan_min_minuit2.h>

/*						minimizer algorithms								*/
/****************************************************************************/

static const stringbuf minalg_id[NMINALG] =
{
	"GSL_GRAD_FR"   ,\
	"GSL_GRAD_PR"   ,\
	"GSL_VEC_BFGS"  ,\
	"GSL_SIMPLEX_NM",\
	"MIN_MIGRAD"    ,\
	"MIN_SIMPLEX"   
};

minalg_no minalg_no_get(const stringbuf m_id)
{
	minalg_no i;
	
	for (i=0;i<NMINALG;i++)
	{
		if (strcmp(m_id,minalg_id[i]) == 0)
		{
			return i;
		}
	}
	LATAN_ERROR("wrong minimizer name",LATAN_FAILURE);
}

latan_errno minalg_id_get(stringbuf m_id, const minalg_no n)
{
	if (n >= NMINALG)
	{
		LATAN_ERROR("wrong minimizer flag",LATAN_EINVAL);
	}
	
	strcpy(m_id,minalg_id[n]);
	
	return LATAN_SUCCESS;
}

/*						minimizer options									*/
/****************************************************************************/

#ifdef HAVE_MINUIT2
#define DEF_LIB MINUIT
#define DEF_ALG MIN_MIGRAD
#else
#define DEF_LIB GSL
#define DEF_ALG GSL_GRAD_FR
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
latan_errno minimize(mat *x, double *f_min, min_func *f, void *param)
{
	latan_errno status;
	stringbuf name;
	
	minimizer_get_alg_name(name);
	latan_printf(VERB,"minimizing using %s algorithm...\n",name);
	switch (minimizer_get_lib())
	{
		case GSL:
			status = minimize_gsl(x,f_min,f,param);
			break;
		case MINUIT:
			status = minimize_minuit2(x,f_min,f,param);
			break;
		default:
			LATAN_ERROR("minimizing library flag invalid",LATAN_EINVAL);
			break;
	}
	return status;
}

