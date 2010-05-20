#include <latan/latan_mass.h>
#include <latan/latan_includes.h>
#include <latan/latan_math.h>

latan_errno effmass(mat res, const mat mprop, const int parity)
{
	size_t i;
	double em;
	
	if (nrow(res) != nrow(mprop) - 2)
	{
		LATAN_ERROR("effective mass matrix have wrong dimensions",\
					LATAN_EBADLEN);
	}
	
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
				em  = MAX(em,1.0);
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

latan_errno effmass_PCAC(mat res, const mat mprop_AP, const mat mprop_PP)
{
	latan_errno status;
	size_t t;
	
	if (!mat_issamedim(mprop_AP,mprop_PP))
	{
		LATAN_ERROR("AP and PP propagators have different dimensions",\
					LATAN_EBADLEN);
	}
	
	status = finite_diff(res,mprop_AP);
	for (t=0;t<nrow(res);t++)
	{
		mat_set(res,t,0,mat_get(res,t,0)/(2.0*mat_get(mprop_PP,t+1,0)));
	}
	mat_eqabs(res);
	
	return status;
}

plat* search_plat(size_t *nplat, const mat data, const mat sigdata,\
				  const size_t ntmax, const double nsig, const double tol)
{
	bool in,toend;
	size_t maxlength,length,cplat,bplat;
	size_t i,j;
	double sup1,sup2,inf1,inf2,re1,re2;
	plat *plat_ar;
	
	in        = false;
	j         = 0;
	length    = 0;
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
				REALLOC_ERRVAL(plat_ar,plat_ar,plat*,*nplat,NULL);
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
