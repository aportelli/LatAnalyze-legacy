/* latan_mass.c, part of LatAnalyze library
 *
 * Copyright (C) 2010, 2011 Antonin Portelli
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
        case EVEN:
            for (i=0;i<nrow(res);i++)
            {
                em = fabs(log(fabs(mat_get(mprop,i+1,0)/mat_get(mprop,i+2,0))));
                mat_set(res,i,0,em);
            }
            break;
        case ODD:
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
if (strcmp(name,#pname) == 0)\
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
    ELIF_NAME(K_iso)
    ELIF_NAME(rho_p)
    ELIF_NAME(rho_0)
    ELIF_NAME(rho_m)
    ELIF_NAME(Kst_p)
    ELIF_NAME(Kst_0)
    ELIF_NAME(Kst_m)
    ELIF_NAME(p)
    ELIF_NAME(n)
    ELIF_NAME(Sigma_p)
    ELIF_NAME(Sigma_0)
    ELIF_NAME(Sigma_m)
    ELIF_NAME(Xi_0)
    ELIF_NAME(Xi_m)
    ELIF_NAME(Delta_pp)
    ELIF_NAME(Delta_p)
    ELIF_NAME(Delta_0)
    ELIF_NAME(Delta_m)
    ELIF_NAME(Sigmast_p)
    ELIF_NAME(Sigmast_0)
    ELIF_NAME(Sigmast_m)
    ELIF_NAME(Xist_0)
    ELIF_NAME(Xist_m)
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

/*                     mass fit parameter tuning functions                  */
/****************************************************************************/
#ifndef PLAT_TOL
#define PLAT_TOL 0.40
#endif
#ifndef NSIGMA
#define NSIGMA 1.0
#endif

plat *search_plat(size_t *nplat, mat *data, mat *sigdata,\
                  const size_t ntmax, const double nsig, const double tol)
{
    bool in,toend;
    size_t maxlength,length,cplat,bplat;
    size_t i,j;
    double sup1,sup2,inf1,inf2,re1,re2;
    plat *plat_ar;
    
    in        = false;
    toend     = false;
    j         = 0;
    maxlength = 0;
    *nplat    = 0;
    cplat     = 0;
    bplat     = 0;
    plat_ar   = NULL;
    
    if (nrow(data) != nrow(sigdata))
    {
        LATAN_ERROR_NULL("data and error vectors must have the same number of rows",\
                         LATAN_EBADLEN);
    }
    latan_printf(VERB,"Plateau searching at with %.2f sigma(s) with %.2f%% of maximum relative error\n",\
                 nsig,tol*100.0);
    for (i=0;i<ntmax;i++)
    {
        sup1 = mat_get(data,i,0)+nsig*mat_get(sigdata,i,0);
        sup2 = mat_get(data,i+1,0)+nsig*mat_get(sigdata,i+1,0);
        inf1 = mat_get(data,i,0)-nsig*mat_get(sigdata,i,0);
        inf2 = mat_get(data,i+1,0)-nsig*mat_get(sigdata,i+1,0);
        re1  = fabs(mat_get(sigdata,i,0)/mat_get(data,i,0));
        re2  = fabs(mat_get(sigdata,i+1,0)/mat_get(data,i+1,0));
        if ((sup1>inf2)&&(sup2>inf1)&&(re1<tol)&&(re2<tol))
        {
            if (!in)
            {
                (*nplat)++;
                cplat = *nplat - 1;
                REALLOC_ERRVAL(plat_ar,plat_ar,plat *,*nplat,NULL);
                plat_ar[cplat].mean  = 0.0;
                plat_ar[cplat].sig   = 0.0;
                plat_ar[cplat].start = i;
                latan_printf(VERB,"Plateau %i\n",*nplat);
                in = true;
            }
            plat_ar[cplat].mean += mat_get(data,i,0);
            plat_ar[cplat].sig  += SQ(mat_get(data,i,0));
            latan_printf(DEBUG,                                                \
                         "%i ]%.10e,%.10e[ val = %e aerr = %e rerr = %.2f%%\n",\
                         i,inf1,sup1,mat_get(data,i,0),mat_get(sigdata,i,0),   \
                         re1*100.0);
            if (i == ntmax-1)
            {
                j = i + 1;
                toend = true;
            }
        }
        else if (in)
        {
            j = i;
            toend = true;
        }
        if (toend)
        {
            plat_ar[cplat].mean += mat_get(data,j,0);
            plat_ar[cplat].sig  += SQ(mat_get(data,j,0));
            plat_ar[cplat].end   = j;
            latan_printf(DEBUG,
                         "%i ]%.10e,%.10e[ val = %e aerr = %e rerr = %.2f%%\n",\
                         j,mat_get(data,j,0)-nsig*mat_get(sigdata,j,0),        \
                         mat_get(data,j,0)+nsig*mat_get(sigdata,j,0),          \
                         mat_get(data,j,0),mat_get(sigdata,j,0),               \
                         fabs(mat_get(sigdata,j,0)/mat_get(data,j,0))*100.0);
            length = plat_ar[cplat].end-plat_ar[cplat].start+1;
            plat_ar[cplat].mean /= (double)(length);
            plat_ar[cplat].sig  /= (double)(length-1);
            plat_ar[cplat].sig  -= ((double)(length)/(double)(length-1))\
                                   *SQ(plat_ar[cplat].mean);
            plat_ar[cplat].sig   = sqrt(plat_ar[cplat].sig);
            latan_printf(VERB,"length = %i mean = %.10e stddev = %.10e\n",\
                         length,plat_ar[cplat].mean,plat_ar[cplat].sig);
            in = false;
            toend = false;
            if (length > maxlength)
            {
                maxlength = length;
                bplat = cplat;
            }
            else if (length == maxlength)
            {
                if (plat_ar[cplat].sig < plat_ar[bplat].sig)
                {
                    bplat = cplat;
                }
            }   
        }
    }
    latan_printf(VERB,"Plateau %i is the best one.\n",bplat+1);
    
    return plat_ar;
}

latan_errno fit_data_mass_fit_tune(fit_data *d, mat *fit_init, mat *prop,\
                                   mat *em, mat *sigem,                  \
                                   const int parity)
{
    plat *em_plat;
    size_t nplat,nt,ntmax;
    size_t p,t;
    double shift,mem,pref;
    strbuf ranges,buf;
    const fit_model *model;
    
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
    fit_data_set_model(d,model,NULL);
    
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
    strbufcpy(ranges,"");
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
    mem = mat_get(em,nt/8-1,0);
    switch (parity)
    {
        case EVEN:
            pref = fabs(mat_get(prop,nt/8,0)*exp((int)(nt)*mem/8));
            break;
        case ODD:
            pref = fabs(mat_get(prop,nt/2,0));
            break;
        default:
            LATAN_ERROR("wrong parity flag",LATAN_EINVAL);
            break;
    }
    latan_printf(VERB,"prefactor = %e mass = %e\n",pref,mem);
    mat_set(fit_init,0,0,mem);
    mat_set(fit_init,1,0,pref);
    FREE(em_plat);
    
    return LATAN_SUCCESS;
}
