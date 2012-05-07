/* latan_mass.c, part of LatAnalyze library
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

#include <latan/latan_mass.h>
#include <latan/latan_includes.h>
#include <latan/latan_math.h>
#include <latan/latan_models.h>
#include <latan/latan_nunits.h>

/*                         effective mass functions                         */
/****************************************************************************/
#ifdef HAVE_ACOSH
extern double acosh(double x); /* acosh is not ANSI compliant */
#endif

latan_errno get_effmass_size(size_t dim[2], const mat *mprop,   \
                             const size_t nstate, const int type)
{
    dim[1] = nstate;
    switch (type) 
    {
        case EM_LOG:
            dim[0] = nrow(mprop) - 2*nstate + 1;
            break;
        case EM_ACOSH:
            dim[0] = nrow(mprop) - 4*nstate + 2;
            break;
        default:
            dim[0] = 0;
            LATAN_ERROR("wrong effective mass type flag",LATAN_EINVAL);
            break;
    }
    
    return LATAN_SUCCESS;
}

/* effective mass for 1 or 2 states ( cf. http://arxiv.org/abs/0903.2314 ) */
latan_errno effmass(mat *res, mat *t, const mat *mprop, const size_t nstate,\
                    const int type)
{
    size_t dim[2],t_0;
    size_t i,j,k,t_i;
    unsigned int bin;
    double *y,*em,lim,sign;
    double (*inv_func)(double);
    
    get_effmass_size(dim,mprop,nstate,type);
    if (res)
    {
        if ((nrow(res) != dim[0])||(ncol(res) != dim[1]))
        {
            LATAN_ERROR("effective mass matrix has wrong dimensions",\
                        LATAN_EBADLEN);
        }
    }
    if (t)
    {
        if (nrow(t) != dim[0])
        {
            LATAN_ERROR("time vector has wrong dimensions",\
                        LATAN_EBADLEN);
        }
    }
    
    MALLOC(y,double *,2*nstate);
    MALLOC(em,double *,nstate);
    
    switch (type)
    {
        case EM_LOG:
            t_0      = 0;
            inv_func = &log;
            lim      = 0.0;
            sign     = -1.0;
            break;
        case EM_ACOSH:
            t_0      = 2*nstate-1;
            inv_func = &acosh;
            lim      = 1.0;
            sign     = 1.0;
            break;
        default:
            t_0      = 0;
            inv_func = NULL;
            lim      = 0.0;
            sign     = 0.0;
            LATAN_ERROR("wrong type flag",LATAN_EINVAL);
            break;
    }
    for (i=0;i<dim[0];i++)
    {
        t_i = t_0 + i;
        if (t)
        {
            mat_set(t,i,0,(double)(t_i));
        }
        if (res)
        {
            for (j=0;j<2*nstate;j++) 
            {
                switch (type)
                {
                    case EM_LOG:
                        y[j] = mat_get(mprop,t_i+j,0);
                        break;
                    case EM_ACOSH:
                        y[j] = 0.0;
                        for (k=0;k<=j;k++)
                        {
                            bin   = binomial((unsigned int)(j),\
                                             (unsigned int)(k));
                            y[j] += ((double)(bin))*mat_get(mprop,t_i+j-2*k,0);
                        }
                        y[j] *= pow(0.5,j);
                        break;
                    default:
                        LATAN_ERROR("wrong type flag",LATAN_EINVAL);
                        break;
                }
            }
            switch (nstate) 
            {
                case 1:
                    em[0] = y[1]/y[0];
                    break;
                case 2:
                    {
                        double a,b,c,d,x1,x2;
                        a     = y[0]*y[2] - SQ(y[1]);
                        b     = y[1]*y[2] - y[0]*y[3];
                        c     = y[1]*y[3] - SQ(y[2]);
                        d     = SQ(b)-4.0*a*c;
                        if (d > 0.0)
                        {
                            x1 = (-b-sqrt(d))/(2.0*a);
                            x2 = (-b+sqrt(d))/(2.0*a);
                            x1 = (x1 > lim) ? x1 : latan_nan();
                            x2 = (x2 > lim) ? x2 : latan_nan();
                            if (latan_isnan(x1)&&!latan_isnan(x2))
                            {
                                em[0] = x2;
                                em[1] = x1;
                            }
                            else if (latan_isnan(x2)&&!latan_isnan(x1))
                            {
                                em[0] = x1;
                                em[1] = x2;
                            }
                            else if (!latan_isnan(x1)&&!latan_isnan(x2))
                            {
                                em[0] = MIN(x1,x2);
                                em[1] = MAX(x1,x2);
                            }
                            else
                            {
                                em[0] = latan_nan();
                                em[1] = latan_nan();
                            }
                        }
                        else
                        {
                            em[0] = latan_nan();
                            em[1] = latan_nan();
                        }
                    }
                    break;
                default:
                    LATAN_ERROR("only 1 or 2 states effective mass implemented",\
                                LATAN_EINVAL);
                    break;
            }
            for (j=0;j<nstate;j++)
            {
                mat_set(res,i,j,sign*inv_func(em[j]));
            }
        }
    }
    
    FREE(y);
    FREE(em);
    
    return LATAN_SUCCESS;
}

#define CENT_VAL(s)  rs_sample_pt_cent_val(s)
#define ELEMENT(s,i) rs_sample_pt_sample(s,i)

latan_errno rs_sample_effmass(rs_sample *s_res, mat *t,                     \
                              const rs_sample *s_mprop, const size_t nstate,\
                              const int type)
{
    size_t i;
    size_t nsample;
    latan_errno status;
    
    if (rs_sample_get_nsample(s_res) != rs_sample_get_nsample(s_mprop))
    {
        LATAN_ERROR("operation between samples with different number of elements",\
                    LATAN_EINVAL);
    }
    
    nsample = rs_sample_get_nsample(s_res);
    status  = LATAN_SUCCESS;
    
    USTAT(effmass(CENT_VAL(s_res),t,CENT_VAL(s_mprop),nstate,type));
    for (i=0;i<nsample;i++)
    {
        USTAT(effmass(ELEMENT(s_res,i),NULL,ELEMENT(s_mprop,i),nstate,type));
    }
    
    return status;
}

#undef CENT_VAL
#undef ELEMENT

latan_errno effmass_PCAC(mat *res, const mat *mprop_AP, const mat *mprop_PP)
{
    latan_errno status;
    size_t t;
    
    if (!mat_is_samedim(mprop_AP,mprop_PP))
    {
        LATAN_ERROR("AP and PP propagators have different dimensions",\
                    LATAN_EBADLEN);
    }
    
    status = finite_diff(res,mprop_AP);
    for (t=0;t<nrow(res);t++)
    {
        mat_set(res,t,0,mat_get(res,t,0)/(2.0*mat_get(mprop_PP,t+1,0)));
    }
    
    return status;
}

/*            interface to read experimental masses in latan_nunits.h       */
/****************************************************************************/
#define IF_NAME(pname)\
if (strbufcmp(name,#pname) == 0)\
{\
    mass[0] = NU_M_##pname;\
    mass[1] = NU_M_##pname##_ERR;\
}
#define ELIF_NAME(pname)\
else IF_NAME(pname)

latan_errno get_mass(double mass[2], const strbuf name)
{
    IF_NAME(pi_p)
    ELIF_NAME(pi_0)
    ELIF_NAME(pi_m)
    ELIF_NAME(pi_iso)
    ELIF_NAME(pi)
    ELIF_NAME(pi_p_miso)
    ELIF_NAME(pi_0_miso)
    ELIF_NAME(K_p)
    ELIF_NAME(K_0)
    ELIF_NAME(K_m)
    ELIF_NAME(K)
    ELIF_NAME(K_iso)
    ELIF_NAME(Kchi)
    ELIF_NAME(rho_p)
    ELIF_NAME(rho_0)
    ELIF_NAME(rho_m)
    ELIF_NAME(rho)
    ELIF_NAME(Kst_p)
    ELIF_NAME(Kst_0)
    ELIF_NAME(Kst_m)
    ELIF_NAME(Kst)
    ELIF_NAME(p)
    ELIF_NAME(n)
    ELIF_NAME(N)
    ELIF_NAME(Sigma_p)
    ELIF_NAME(Sigma_0)
    ELIF_NAME(Sigma_m)
    ELIF_NAME(Sigma)
    ELIF_NAME(Xi_0)
    ELIF_NAME(Xi_m)
    ELIF_NAME(Xi)
    ELIF_NAME(Delta_pp)
    ELIF_NAME(Delta_p)
    ELIF_NAME(Delta_0)
    ELIF_NAME(Delta_m)
    ELIF_NAME(Delta)
    ELIF_NAME(Sigmast_p)
    ELIF_NAME(Sigmast_0)
    ELIF_NAME(Sigmast_m)
    ELIF_NAME(Sigmast)
    ELIF_NAME(Xist_0)
    ELIF_NAME(Xist_m)
    ELIF_NAME(Xist)
    ELIF_NAME(Omega_m)
    ELIF_NAME(Omega)
    else
    {
        LATAN_ERROR("particle name unknown",LATAN_EINVAL);
    }
    
    return LATAN_SUCCESS;
}

#undef IF_NAME
#undef ELIF_NAME
