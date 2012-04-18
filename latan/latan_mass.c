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

latan_errno effmass(mat *res, const mat *mprop, const int parity)
{
    size_t i;
    double em;
    
    if (nrow(res) != nrow(mprop) - 2)
    {
        LATAN_ERROR("effective mass matrix have wrong dimensions",\
                    LATAN_EBADLEN);
    }
    
    /* fabs/MAX(.,1.0) : you don't want NaN at half-time */
    switch (parity)
    {
        case EM_LOG:
            for (i=0;i<nrow(res);i++)
            {
                em = fabs(log(fabs(mat_get(mprop,i+1,0)/mat_get(mprop,i+2,0))));
                mat_set(res,i,0,em);
            }
            break;
        case EM_ACOSH:
            for (i=0;i<nrow(res);i++)
            {
                em  = mat_get(mprop,i,0) + mat_get(mprop,i+2,0);
                em /= 2.0*mat_get(mprop,i+1,0);
                em  = MAX(em,cosh(1.0e-3));
                em  = acosh(em);
                mat_set(res,i,0,em);
            }
            break;
        default:
            LATAN_ERROR("wrong parity flag",LATAN_EINVAL);
            break;
    }
    
    return LATAN_SUCCESS;
}

#define CENT_VAL(s)  rs_sample_pt_cent_val(s)
#define ELEMENT(s,i) rs_sample_pt_sample(s,i)

latan_errno rs_sample_effmass(rs_sample *s_res, const rs_sample *s_mprop,\
                              const int parity)
{
    size_t i;
    size_t nsample;
    latan_errno status;
    
    if (rs_sample_get_nsample(s_res) != rs_sample_get_nsample(s_mprop))
    {
        LATAN_ERROR("operation between samples with different numbers of elements",\
                    LATAN_EINVAL);
    }
    
    nsample = rs_sample_get_nsample(s_res);
    status  = LATAN_SUCCESS;
    
    USTAT(effmass(CENT_VAL(s_res),CENT_VAL(s_mprop),parity));
    for (i=0;i<nsample;i++)
    {
        USTAT(effmass(ELEMENT(s_res,i),ELEMENT(s_mprop,i),parity));
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
