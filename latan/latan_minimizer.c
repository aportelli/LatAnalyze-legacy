/* latan_minimizer.c, part of LatAnalyze library
 *
 * Copyright (C) 2010, 2011, 2012 Antonin Portelli
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <latan/latan_minimizer.h>
#include <latan/latan_includes.h>
#include <latan/latan_min_gsl.h>
#ifdef HAVE_MINUIT2
#include <latan/latan_min_minuit2.h>
#endif

/*                      minimizer algorithms                                */
/****************************************************************************/

static const strbuf minalg_id[NMINALG] =
{
    "GSL_GRAD_FR"   ,\
    "GSL_GRAD_PR"   ,\
    "GSL_VEC_BFGS"  ,\
    "GSL_SIMPLEX_NM",\
    "MIN_MIGRAD"    ,\
    "MIN_SIMPLEX"   
};

minalg_no minalg_no_get(const strbuf m_id)
{
    unsigned int i;
    
    for (i=0;i<NMINALG;i++)
    {
        if (strbufcmp(m_id,minalg_id[i]) == 0)
        {
            return (minalg_no)i;
        }
    }
    LATAN_ERROR_VAL("wrong minimizer name",LATAN_EINVAL,(minalg_no)0);
}

latan_errno minalg_id_get(strbuf m_id, const minalg_no n)
{
    if ((unsigned int)n >= NMINALG)
    {
        LATAN_ERROR("wrong minimizer flag",LATAN_EINVAL);
    }
    
    strbufcpy(m_id,minalg_id[n]);
    
    return LATAN_SUCCESS;
}

/*                      minimizer options                                   */
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

latan_errno minimizer_get_alg_name(strbuf name)
{
    switch (env.alg)
    {
        case GSL_GRAD_FR:
            strbufcpy(name,"Fletcher-Reeves conjugate gradient (GSL)");
            break;
        case GSL_GRAD_PR:
            strbufcpy(name,"Polak-Ribiere conjugate gradient (GSL)");
            break;
        case GSL_VEC_BFGS:
            strbufcpy(name,\
                   "improved Broyden-Fletcher-Goldfarb-Shanno vector (GSL)");
            break;
        case GSL_SIMPLEX_NM:
            strbufcpy(name,"improved Nelder-Mead simplex (GSL)");
            break;
        case MIN_MIGRAD:
            strbufcpy(name,"variable metric (MINUIT)");
            break;
        case MIN_SIMPLEX:
            strbufcpy(name,"simplex (MINUIT)");
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

/*                          the minimizer                                   */
/****************************************************************************/
latan_errno minimize(mat *x, const mat *x_limit, double *f_min, min_func *f,\
                     void *param)
{
    latan_errno status;
    strbuf name;
    
    minimizer_get_alg_name(name);
    latan_printf(VERB,"minimizing using %s algorithm...\n",name);
    switch (minimizer_get_lib())
    {
        case GSL:
            status = minimize_gsl(x,x_limit,f_min,f,param);
            break;
        case MINUIT:
#ifdef HAVE_MINUIT2
            status = minimize_minuit2(x,x_limit,f_min,f,param);
#else
            LATAN_ERROR("MINUIT support was not compiled",LATAN_EINVAL);
#endif
            break;
        default:
            LATAN_ERROR("minimizing library flag invalid",LATAN_EINVAL);
            break;
    }
    return status;
}

