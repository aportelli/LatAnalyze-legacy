#include <latan/hadron.h>
#include <latan/includes.h>
#include <latan/io.h>
#include <latan/statistics.h>

static stringbuf propmark = "PROP";
static stringbuf propidfmt = "%s_%s_%d_%d";

int hadron_getnt(const hadron h, const int source, const int sink,\
				 const stringbuf manfname)
{
	stringbuf ffname,fullpropid;
	int nt;
	
	get_firstfname(ffname,manfname);
	sprintf(fullpropid,propidfmt,h->channel[0],h->quarkst[0],source,sink);
	nt = mat_load_nrow(propmark,fullpropid,ffname);
	
	return nt;
}

latan_errno hadron_prop(mat* prop, const hadron h, const int source,\
						const int sink, const stringbuf manfname)
{
	int i,p,s;
	size_t j;
	int ndat, chmix, stmix;
	size_t nt;
	double mean;
	mat* dat[MAXPROP][MAXQUARKST];
	stringbuf fullpropid;
	latan_errno status;
	
	nt		= nrow(prop[0]);
	ndat	= get_nfile(manfname);
	chmix	= ((h->chmix) == NOMIX) ? 1 : MAXPROP;
	stmix	= ((h->stmix) == NOMIX) ? 1 : MAXQUARKST;
	status	= LATAN_SUCCESS;
	
	if (ndat == -1)
	{
		LATAN_ERROR("error while reading manifest file",LATAN_ESYSTEM);
	}
	
	for (p=0;p<chmix;p++)
	{
		for (s=0;s<stmix;s++)
		{
			dat[p][s] = mat_create_ar((size_t)(ndat),(size_t)(nt),1);
			sprintf(fullpropid,propidfmt,h->channel[p],h->quarkst[s],source,\
					sink);
			LATAN_UPDATE_STATUS(status,mat_load_ar(dat[p][s],propmark,\
												   fullpropid,manfname));
			for (i=0;i<ndat;i++)
			{
				LATAN_UPDATE_STATUS(status,mat_eqabs(dat[p][s][i]));
			}
		}
	}
	
	for (i=0;i<ndat;i++)
	{
		mat_zero(prop[i]);
		for (p=0;p<chmix;p++)
		{
			for (s=0;s<stmix;s++)
			{
				LATAN_UPDATE_STATUS(status,mat_eqadd(prop[i],dat[p][s][i]));
			}
		}
		if ((h->chmix) == MEAN)
		{
			LATAN_UPDATE_STATUS(status,mat_eqmuls(prop[i],DRATIO(1,MAXPROP)));
		}
		if ((h->stmix) == MEAN)
		{
			LATAN_UPDATE_STATUS(status,mat_eqmuls(prop[i],\
												  DRATIO(1,MAXQUARKST)));
		}
	}
	
	if ((h->parity) == ODD)
	{
		for (i=0;i<ndat;i++)
		{
			for (j=1;j<nt/2;j++)
			{
				mean = 0.5*(mat_get(prop[i],(size_t)(j),0)	\
							+ mat_get(prop[i],(size_t)(nt-j),0));
				mat_set(prop[i],(size_t)(j),0,mean);
				mat_set(prop[i],(size_t)(nt-j),0,mean);
			}
		}
	}
	
	for (p=0;p<chmix;p++)
	{
		for (s=0;s<stmix;s++)
		{
			mat_destroy_ar(dat[p][s],(size_t)(ndat));
		}
	}
	
	return status;
}