#include <latan/latan_hadron.h>
#include <latan/latan_includes.h>

/* channels */
static stringbuf channel_id[] = 
{
	"SS"	,\
	"ViVi"	,\
	"PP"	,\
	"PA4"	,\
	"A4P"	,\
	"A4A4"	,\
	"NUCL"	,\
	"LAMBDA",\
	"DELTA"
};

channel_no channel_no_get(const stringbuf label)
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

void channel_id_set(const channel_no i, const stringbuf new_id)
{
	strcpy(channel_id[i],new_id);
}

void channel_id_get(stringbuf str, const channel_no i)
{
	strcpy(str,channel_id[i]);
}

/* quarks */
static stringbuf quark_id[] =
{
	"0"	,\
	"0"	,\
	"1"	,\
	"2"
};

quark_no quark_no_get(const char c)
{
	quark_no res;
	
	switch (c)
	{
		case 'l':
			res = qu_l;
			break;
		case 'u':
			res = qu_u;
			break;
		case 'd':
			res = qu_d;
			break;
		case 's':
			res = qu_s;
			break;
		default:
			LATAN_ERROR("wrong quark name",LATAN_FAILURE);
			break;
	}
	
	return res;
}

void quark_id_set(const quark_no i, const stringbuf new_id)
{
	strcpy(quark_id[i],new_id);
}

void quark_id_get(stringbuf str, const quark_no i)
{
	strcpy(str,quark_id[i]);
}

void diquark_id_get(stringbuf str, const quark_no i1, const quark_no i2)
{
	sprintf(str,"%s%s",quark_id[i1],quark_id[i2]);
}

void triquark_id_get(stringbuf str, const quark_no i1, const quark_no i2,\
					 const quark_no i3)
{
	sprintf(str,"%s%s%s",quark_id[i1],quark_id[i2],quark_id[i3]);
}

/* sources/sinks */
static stringbuf ss_id[] =
{
	"0",
	"1",
	"2"
};

ss_no ss_no_get(const char c)
{
	ss_no res;
	
	switch (c)
	{
		case 'P':
			res = ss_P;
			break;
		case 'W':
			res = ss_W;
			break;
		case 'G':
			res = ss_G;
			break;
		default:
			LATAN_ERROR("wrong source/sink name",LATAN_FAILURE);
			break;
	}
	
	return res;
}

void ss_id_set(const ss_no i, const stringbuf new_id)
{
	strcpy(ss_id[i],new_id);
}

void ss_id_get(stringbuf str, const ss_no i)
{
	strcpy(str,ss_id[i]);
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
void hadron_set_2q_nomix(hadron *h, const stringbuf name, const int parity,	\
						 const channel_no channel, const quark_no q1,		\
						 const quark_no q2)
{
	strcpy(h->name,name);
	h->parity		= parity;
	h->chmix		= NOMIX;
	h->stmix		= NOMIX;
	channel_id_get(h->channel[0],channel);
	diquark_id_get(h->quarkst[0],q1,q2);
}

void hadron_set_2q_2stmean(hadron *h, const stringbuf name, const int parity,\
						   const channel_no channel, const quark_no q11,	\
						   const quark_no q12, const quark_no q21,			\
						   const quark_no q22)
{
	strcpy(h->name,name);
	h->parity		= parity;
	h->chmix		= NOMIX;
	h->stmix		= MEAN;
	channel_id_get(h->channel[0],channel);
	diquark_id_get(h->quarkst[0],q11,q12);
	diquark_id_get(h->quarkst[1],q21,q22);
}

void hadron_get_name(stringbuf str, const hadron *h)
{
	strcpy(str,h->name);
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
	h_pi		= 0,
	h_K			= 1,
	h_rho		= 2,
	h_Kst		= 3,
	h_N			= 4,
	h_Lambda	= 5,
	h_Sigma		= 6,
	h_Xi		= 7,
	h_Delta		= 8,
	h_Sigmast	= 9,
	h_Xist		= 10,
	h_Omega		= 11
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
#define SPECT_QCDQED_SIZE 26
enum
{
	h_pi_0		= 0,
	h_pi_p		= 1,
	h_K_0		= 2,
	h_K_p		= 3,
	h_rho_0		= 4,
	h_rho_p		= 5,
	h_Kst_0		= 6,
	h_Kst_p		= 7,
	h_n			= 8,
	h_p			= 9,
	h_Sigma_m	= 10,
	h_Sigma_p	= 11,
	h_Xi_m		= 12,
	h_Xi_0		= 13,
	h_Delta_m	= 14,
	h_Delta_0	= 15,
	h_Delta_p	= 16,
	h_Delta_pp	= 17,
	h_Sigmast_m	= 18,
	h_Sigmast_p = 19,
	h_Xist_m	= 20,
	h_Xist_0	= 21,
	h_PP_uu		= 22,
	h_PP_dd		= 23,
	h_VV_uu		= 24,
	h_VV_dd		= 25
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
	hadron_set_2q_nomix(s_qcdqed->particle[h_VV_uu],"VV_uu",ODD,ch_VV,qu_u,\
						qu_u);
	hadron_set_2q_nomix(s_qcdqed->particle[h_VV_dd],"VV_dd",ODD,ch_VV,qu_d,\
						qu_d);
	
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
	return s_qcdqed;
}

/** access **/
hadron *spectrum_get(const spectrum *s, const stringbuf part_name)
{
	size_t i;
	hadron *h;
	bool found;
	stringbuf errmsg;
	
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


