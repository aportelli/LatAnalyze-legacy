/* latan_hadron.c, part of LatAnalyze library
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

#include <latan/latan_hadron.h>
#include <latan/latan_includes.h>

/* channels */
static strbuf channel_id[] = 
{
    "SS"    ,\
    "ViVi"  ,\
    "PP"    ,\
    "PA4"   ,\
    "A4P"   ,\
    "AiAi"  ,\
    "NUCL"  ,\
    "LAMBDA",\
    "DELTA"  \
};

channel_no channel_get_no_from_label(const strbuf label)
{
    channel_no res;
    
    if (strcmp(label,"SS") == 0)
    {
        res = ch_SS;
    }
    else if (strcmp(label,"VV") == 0)
    {
        res = ch_VV;
    }
    else if (strcmp(label,"PP") == 0)
    {
        res = ch_PP;
    }
    else if (strcmp(label,"PA") == 0)
    {
        res = ch_PA;
    }
    else if (strcmp(label,"AP") == 0)
    {
        res = ch_AP;
    }
    else if (strcmp(label,"AA") == 0)
    {
        res = ch_AA;
    }
    else if (strcmp(label,"N") == 0)
    {
        res = ch_N;
    }
    else if (strcmp(label,"Lambda") == 0)
    {
        res = ch_Lambda;
    }
    else if (strcmp(label,"Delta") == 0)
    {
        res = ch_Delta;
    }
    else
    {
        LATAN_ERROR("wrong channel name",LATAN_FAILURE);
    }
    
    return res;
}

channel_no channel_get_no_from_id(const strbuf id)
{
    channel_no res;

    if (strcmp(id,channel_id[ch_SS]) == 0)
    {
        res = ch_SS;
    }
    else if (strcmp(id,channel_id[ch_VV]) == 0)
    {
        res = ch_VV;
    }
    else if (strcmp(id,channel_id[ch_PP]) == 0)
    {
        res = ch_PP;
    }
    else if (strcmp(id,channel_id[ch_PA]) == 0)
    {
        res = ch_PA;
    }
    else if (strcmp(id,channel_id[ch_AP]) == 0)
    {
        res = ch_AP;
    }
    else if (strcmp(id,channel_id[ch_AA]) == 0)
    {
        res = ch_AA;
    }
    else if (strcmp(id,channel_id[ch_N]) == 0)
    {
        res = ch_N;
    }
    else if (strcmp(id,channel_id[ch_Lambda]) == 0)
    {
        res = ch_Lambda;
    }
    else if (strcmp(id,channel_id[ch_Delta]) == 0)
    {
        res = ch_Delta;
    }
    else
    {
        LATAN_ERROR("unknown channel ID",LATAN_FAILURE);
    }

    return res;
}

void channel_id_set(const channel_no i, const strbuf new_id)
{
    strbufcpy(channel_id[i],new_id);
}

void channel_id_get(strbuf str, const channel_no i)
{
    strbufcpy(str,channel_id[i]);
}

/* quarks */
static strbuf quark_id[] =
{
    "0" ,\
    "0" ,\
    "1" ,\
    "2"
};

void quark_id_set(const quark_no i, const strbuf new_id)
{
    strbufcpy(quark_id[i],new_id);
}

void quark_id_get(strbuf str, const quark_no i)
{
    strbufcpy(str,quark_id[i]);
}

/* sources/sinks */
static strbuf ss_id[] =
{
    "P",
    "W",
    "G"
};

void ss_id_set(const ss_no i, const strbuf new_id)
{
    strbufcpy(ss_id[i],new_id);
}

void ss_id_get(strbuf str, const ss_no i)
{
    strbufcpy(str,ss_id[i]);
}

/* hadron */
/** allocation **/
hadron *hadron_create(void)
{
    hadron *h;
    
    MALLOC_ERRVAL(h,hadron *,1,NULL);
    
    return h;
}

void hadron_destroy(hadron *h)
{
    FREE(h);
}

/** access **/
void hadron_set_2q_nomix(hadron *h, const strbuf name, const int parity,  \
                         const channel_no channel, const quark_no q1,     \
                         const quark_no q2)
{
    strbufcpy(h->name,name);
    h->parity        = parity;
    h->chmix         = NOMIX;
    h->stmix         = NOMIX;
    h->channel[0]    = channel;
    h->quarkst[0][0] = q1;
    h->quarkst[0][1] = q2;
}

void hadron_set_2q_2stmean(hadron *h, const strbuf name, const int parity,\
                           const channel_no channel, const quark_no q11,  \
                           const quark_no q12, const quark_no q21,        \
                           const quark_no q22)
{
    strbufcpy(h->name,name);
    h->parity        = parity;
    h->chmix         = NOMIX;
    h->stmix         = MEAN;
    h->channel[0]    = channel;
    h->quarkst[0][0] = q11;
    h->quarkst[0][1] = q12;
    h->quarkst[1][0] = q21;
    h->quarkst[1][1] = q22;
}

void hadron_get_name(strbuf str, const hadron *h)
{
    strbufcpy(str,h->name);
}

/* spectrum */
/** allocation **/
spectrum *spectrum_create(const size_t nparticle)
{
    spectrum *s;
    size_t i;
    
    MALLOC_ERRVAL(s,spectrum *,1,NULL);
    
    s->nparticle = nparticle;
    
    MALLOC_ERRVAL(s->particle,hadron **,s->nparticle,NULL);
    for (i=0;i<s->nparticle;i++)
    {
        s->particle[i] = hadron_create();
    }
    
    return s;
}

void spectrum_destroy(spectrum *s)
{
    size_t i;
    
    for (i=0;i<s->nparticle;i++)
    {
        hadron_destroy(s->particle[i]);
    }
    FREE(s->particle);
    s->nparticle = 0;
    FREE(s);
}

/** QCD spectrum **/
#define SPECT_QCD_SIZE 12
enum
{
    h_pi        = 0,
    h_K         = 1,
    h_rho       = 2,
    h_Kst       = 3,
    h_N         = 4,
    h_Lambda    = 5,
    h_Sigma     = 6,
    h_Xi        = 7,
    h_Delta     = 8,
    h_Sigmast   = 9,
    h_Xist      = 10,
    h_Omega     = 11
};

spectrum *spectrum_create_qcd(void)
{
    spectrum *s_qcd;
    
    s_qcd = spectrum_create(SPECT_QCD_SIZE);
    
    /* mesons */
    hadron_set_2q_nomix(s_qcd->particle[h_pi],"pi",ODD,ch_PP,qu_l,qu_l);
    hadron_set_2q_nomix(s_qcd->particle[h_K],"K",ODD,ch_PP,qu_l,qu_s);
    hadron_set_2q_nomix(s_qcd->particle[h_rho],"rho",ODD,ch_VV,qu_l,qu_l);
    hadron_set_2q_nomix(s_qcd->particle[h_Kst],"K*",ODD,ch_VV,qu_l,qu_s);
    
    /* baryons */
    hadron_set_2q_nomix(s_qcd->particle[h_N],"N",EVEN,ch_N,qu_l,qu_l);
    hadron_set_2q_nomix(s_qcd->particle[h_Lambda],"Lambda",EVEN,ch_Lambda,\
                        qu_s,qu_l);
    hadron_set_2q_nomix(s_qcd->particle[h_Sigma],"Sigma",EVEN,ch_N,qu_s,qu_l);
    hadron_set_2q_nomix(s_qcd->particle[h_Xi],"Xi",EVEN,ch_N,qu_l,qu_s);
    hadron_set_2q_nomix(s_qcd->particle[h_Delta],"Delta",EVEN,ch_Delta,qu_l,\
                        qu_l);
    hadron_set_2q_nomix(s_qcd->particle[h_Sigmast],"Sigma*",EVEN,ch_Delta,qu_s,\
                        qu_l);
    hadron_set_2q_nomix(s_qcd->particle[h_Xist],"Xi*",EVEN,ch_Delta,qu_l,qu_s);
    hadron_set_2q_nomix(s_qcd->particle[h_Omega],"Omega",EVEN,ch_Delta,qu_s,\
                        qu_s);
    
    return s_qcd;
}

/** QCD+QED spectrum **/
#define SPECT_QCDQED_SIZE 29
enum
{
    h_pi_0      = 0,
    h_pi_p      = 1,
    h_K_0       = 2,
    h_K_p       = 3,
    h_rho_0     = 4,
    h_rho_p     = 5,
    h_Kst_0     = 6,
    h_Kst_p     = 7,
    h_n         = 8,
    h_p         = 9,
    h_Sigma_m   = 10,
    h_Sigma_p   = 11,
    h_Xi_m      = 12,
    h_Xi_0      = 13,
    h_Delta_m   = 14,
    h_Delta_0   = 15,
    h_Delta_p   = 16,
    h_Delta_pp  = 17,
    h_Sigmast_m = 18,
    h_Sigmast_p = 19,
    h_Xist_m    = 20,
    h_Xist_0    = 21,
    h_Omega_m   = 22,
    h_PP_uu     = 23,
    h_PP_dd     = 24,
    h_PP_ss     = 25,
    h_VV_uu     = 26,
    h_VV_dd     = 27,
    h_VV_ss     = 28
};

spectrum *spectrum_create_qcdqed(void)
{
    spectrum *s_qcdqed;
    
    s_qcdqed = spectrum_create(SPECT_QCDQED_SIZE);
    
    /* mesons */
    hadron_set_2q_2stmean(s_qcdqed->particle[h_pi_0],"pi0",ODD,ch_PP,\
                          qu_u,qu_u,qu_d,qu_d);
    hadron_set_2q_nomix(s_qcdqed->particle[h_pi_p],"pi+",ODD,ch_PP,qu_u,qu_d);
    hadron_set_2q_nomix(s_qcdqed->particle[h_K_0],"K0",ODD,ch_PP,qu_d,qu_s);
    hadron_set_2q_nomix(s_qcdqed->particle[h_K_p],"K+",ODD,ch_PP,qu_u,qu_s);
    hadron_set_2q_2stmean(s_qcdqed->particle[h_rho_0],"rho0",ODD,ch_VV,\
                          qu_u,qu_u,qu_d,qu_d);
    hadron_set_2q_nomix(s_qcdqed->particle[h_rho_p],"rho+",ODD,ch_VV,qu_u,qu_d);
    hadron_set_2q_nomix(s_qcdqed->particle[h_Kst_0],"K*0",ODD,ch_VV,qu_d,qu_s);
    hadron_set_2q_nomix(s_qcdqed->particle[h_Kst_p],"K*+",ODD,ch_VV,qu_u,qu_s);
    hadron_set_2q_nomix(s_qcdqed->particle[h_PP_uu],"PP_uu",ODD,ch_PP,qu_u,\
                        qu_u);
    hadron_set_2q_nomix(s_qcdqed->particle[h_PP_dd],"PP_dd",ODD,ch_PP,qu_d,\
                        qu_d);
    hadron_set_2q_nomix(s_qcdqed->particle[h_PP_ss],"PP_ss",ODD,ch_PP,qu_s,\
                        qu_s);
    hadron_set_2q_nomix(s_qcdqed->particle[h_VV_uu],"VV_uu",ODD,ch_VV,qu_u,\
                        qu_u);
    hadron_set_2q_nomix(s_qcdqed->particle[h_VV_dd],"VV_dd",ODD,ch_VV,qu_d,\
                        qu_d);
    hadron_set_2q_nomix(s_qcdqed->particle[h_VV_ss],"VV_ss",ODD,ch_VV,qu_s,\
                        qu_s);
    
    /* baryons */
    hadron_set_2q_nomix(s_qcdqed->particle[h_n],"n",EVEN,ch_N,qu_u,qu_d);
    hadron_set_2q_nomix(s_qcdqed->particle[h_p],"p",EVEN,ch_N,qu_d,qu_u);
    hadron_set_2q_nomix(s_qcdqed->particle[h_Sigma_m],"Sigma-",EVEN,ch_N,qu_s,\
                        qu_d);
    hadron_set_2q_nomix(s_qcdqed->particle[h_Sigma_p],"Sigma+",EVEN,ch_N,qu_s,\
                        qu_u);
    hadron_set_2q_nomix(s_qcdqed->particle[h_Xi_m],"Xi-",EVEN,ch_N,qu_d,qu_s);
    hadron_set_2q_nomix(s_qcdqed->particle[h_Xi_0],"Xi0",EVEN,ch_N,qu_u,qu_s);
    hadron_set_2q_nomix(s_qcdqed->particle[h_Delta_m],"Delta-",EVEN,ch_Delta,\
                  qu_d,qu_d);
    hadron_set_2q_nomix(s_qcdqed->particle[h_Delta_0],"Delta0",EVEN,ch_Delta,\
                  qu_u,qu_d);
    hadron_set_2q_nomix(s_qcdqed->particle[h_Delta_p],"Delta+",EVEN,ch_Delta,\
                  qu_d,qu_u);
    hadron_set_2q_nomix(s_qcdqed->particle[h_Delta_pp],"Delta++",EVEN,ch_Delta,\
                  qu_u,qu_u);
    hadron_set_2q_nomix(s_qcdqed->particle[h_Sigmast_m],"Sigma*-",EVEN,\
                        ch_Delta,qu_s,qu_d);
    hadron_set_2q_nomix(s_qcdqed->particle[h_Sigmast_p],"Sigma*+",EVEN,\
                        ch_Delta,qu_s,qu_u);
    hadron_set_2q_nomix(s_qcdqed->particle[h_Xist_m],"Xi*-",EVEN,ch_Delta,qu_d,\
                        qu_s);
    hadron_set_2q_nomix(s_qcdqed->particle[h_Xist_0],"Xi*0",EVEN,ch_Delta,qu_u,\
                        qu_s);
    hadron_set_2q_nomix(s_qcdqed->particle[h_Omega_m],"Omega-",EVEN,ch_Delta,\
                        qu_s,qu_s);
    return s_qcdqed;
}

/** access **/
hadron *spectrum_get(const spectrum *s, const strbuf part_name)
{
    size_t i;
    hadron *h;
    bool found;
    strbuf errmsg;
    
    found = false;
    
    for (i=0;i<s->nparticle;i++)
    {
        if (strcmp(s->particle[i]->name,part_name) == 0)
        {
            h = s->particle[i];
            found = true;
            break;
        }
    }
    if (!found)
    {
        sprintf(errmsg,"hadron %s not found",part_name);
        LATAN_ERROR_NULL(errmsg,LATAN_EINVAL);
    }
    
    return h;
}


