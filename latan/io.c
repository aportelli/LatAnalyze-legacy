#include <latan/io.h>
#include <latan/includes.h>

/*							general I/O										*/
/****************************************************************************/

int get_nfile(const stringbuf manifestfname)
{
	stringbuf buf1, buf2;
	int nfile;
	FILE* manifest = NULL;
	
	nfile = 0;
	
	FOPEN_ERRVAL(manifest,manifestfname,"r",LATAN_FAILURE);
	while (!feof(manifest))
	{
		if ((fgets(buf1,STRING_LENGTH,manifest))&&(sscanf(buf1,"%s\n",buf2)>0))
		{
			nfile++;
		}
	}
	fclose(manifest);
	
	return nfile;
}

latan_errno get_firstfname(stringbuf fname, const stringbuf manifestfname)
{
	stringbuf buf;
	FILE* manifest = NULL;
	
	FOPEN(manifest,manifestfname,"r");
	fgets(buf,STRING_LENGTH,manifest);
	fclose(manifest);
	sscanf(buf,"%s\n",fname);
	
	return LATAN_SUCCESS;
}

/*								mat I/O										*/
/****************************************************************************/

void mat_dump(FILE* stream, const mat m)
{
	size_t i,j;
	
	for (i=0;i<nrow(m);i++)
	{
		for (j=0;j<ncol(m)-1;j++)
		{
			fprintf(stream,"%.10e ",mat_get(m,i,j));
		}
		fprintf(stream,"%.10e\n",mat_get(m,i,ncol(m)-1));
	}
}

int mat_load_nrow(const stringbuf mark, const stringbuf matid,\
				  const stringbuf inputfname)
{
	stringbuf buf1, buf2, startfmt, end;
	double dumb;
	int dat_nrow;
	bool in, found;
	FILE* inputf = NULL;
	
	dat_nrow = 0;
	in = false;
	found = false;
	
	sprintf(startfmt,"START_%s %%s",mark);
	sprintf(end,"END_%s",mark);
	FOPEN_ERRVAL(inputf,inputfname,"r",LATAN_FAILURE);
	while (!feof(inputf))
	{
		fgets(buf1,STRING_LENGTH,inputf);
		if (!in)
		{
			if (sscanf(buf1,startfmt,buf2)>0)
			{
				if (strcmp(buf2,matid)==0)
				{
					in = true;
					if (!found)
					{
						found = true;
					}
				}
			}	
		}
		else
		{
			if (sscanf(buf1,"%lf",&dumb)>0)
			{
				dat_nrow++;
			}
			else if (sscanf(buf1,"%s",buf2)>0)
			{
				if (strcmp(buf2,end)==0)
				{
					break;
				}
				else 
				{
					stringbuf errmsg;
					
					sprintf(errmsg,"unexpected expression %s found before END_%s",\
							buf2,mark);
					LATAN_ERROR_VAL(errmsg,LATAN_ELATSYN,LATAN_FAILURE);
				}
			}
		}
	}
	fclose(inputf);
	if (!found)
	{
		stringbuf errmsg ;
		
		sprintf(errmsg,"no matrix %s %s in file %s",mark,matid,inputfname);
		LATAN_ERROR_VAL(errmsg,LATAN_ELATSYN,LATAN_FAILURE);
	}
	
	return dat_nrow;
}

latan_errno mat_load(mat m, const stringbuf mark, const stringbuf matid,\
					 const stringbuf inputfname)
{
	stringbuf buf1, buf2, startfmt, datfmt, tmp, end;
	double buf;
	size_t row, col;
	bool in, broken, found;
	FILE* inputf = NULL;
	
	row = 0;
	in = false;
	broken = false;
	found = false;
	
	sprintf(startfmt,"START_%s %%s",mark);
	sprintf(end,"END_%s",mark);
	FOPEN(inputf,inputfname,"r");
	while (!feof(inputf))
	{
		fgets(buf1,STRING_LENGTH,inputf);
		if (!in)
		{
			if (sscanf(buf1,startfmt,buf2)>0)
			{
				if (strcmp(buf2,matid)==0)
				{
					in = true;
					if (!found)
					{
						found = true;
					}
				}
			}	
		}
		else
		{
			strcpy(datfmt,"%lf");
			for (col=0;col<ncol(m);col++)
			{
				if (sscanf(buf1,datfmt,&buf)>0)
				{
					mat_set(m,row,col,buf);
					strcpy(tmp,datfmt);
					strcpy(datfmt,"%*lf ");
					strcat(datfmt,tmp);
				}
				else if (sscanf(buf1,"%s",buf2)>0)
				{
					if (strcmp(buf2,end)==0)
					{
						broken = true;
						break;
					}
					else 
					{
						stringbuf errmsg;
						
						sprintf(errmsg,"error while reading matrix %s %s (%u,%u) in file %s (trying to read a %ux%u matrix)",\
								mark,matid,(unsigned)(row),(unsigned)(col),	\
								inputfname,(unsigned)(nrow(m)),				\
								(unsigned)(ncol(m)));
						LATAN_ERROR(errmsg,LATAN_ELATSYN);
					}
				}
			}
			row++;
			if (broken)
			{
				break;
			}
		}
	}
	fclose(inputf);
	if (!found)
	{
		stringbuf errmsg;
		
		sprintf(errmsg,"no matrix %s %s in file %s",mark,matid,inputfname);
		LATAN_ERROR(errmsg,LATAN_ELATSYN);
	}
	
	return LATAN_SUCCESS;
}

latan_errno mat_load_ar(mat* m, const stringbuf mark, const stringbuf matid,\
						const stringbuf manifestfname)
{
	stringbuf buf, fname;
	int nfile, status;
	int i;
	FILE* manifest = NULL;
	
	status = LATAN_SUCCESS;
	nfile = get_nfile(manifestfname);
	
	FOPEN(manifest,manifestfname,"r");
	for (i=0;i<nfile;i++)
	{
		do
		{
			fgets(buf,STRING_LENGTH,manifest);
		} while (sscanf(buf,"%s\n",fname)<=0);
		LATAN_UPDATE_STATUS(status,mat_load(m[i],mark,matid,fname));
	}
	fclose(manifest);
	
	return status;
}

latan_errno mat_save_plotdat(const mat dat, const double xstart,\
							 const double xstep, const stringbuf fname)
{
	FILE* f = NULL;
	size_t i;
	
	FOPEN(f,fname,"w");
	
	for (i=0;i<nrow(dat);i++)
	{
		fprintf(f,"%.10e %.10e\n",(double)(xstart+(int)(i)*xstep),\
				mat_get(dat,i,0));
	}
	fclose(f);
	
	return LATAN_SUCCESS;
}

latan_errno mat_save_plotdaterr(const mat dat, const mat sig,				\
								const double xstart, const double xstep,	\
								const stringbuf fname)
{
	FILE* f = NULL;
	size_t i;
	
	FOPEN(f,fname,"w");
	for (i=0;i<nrow(dat);i++)
	{
		fprintf(f,"%.10e %.10e %.10e\n",(double)(xstart+(int)(i)*xstep),\
				mat_get(dat,i,0),mat_get(sig,i,0));
	}
	fclose(f);
	
	return LATAN_SUCCESS;
}

/*						random generator state I/O							*/
/****************************************************************************/

latan_errno randgen_save_state(const stringbuf prefname,\
							   const randgen_state state)
{
	stringbuf fname;
	FILE *f;
	int i;
	
	sprintf(fname,"%s.rand",prefname);
	FOPEN(f,fname,"w");
	for (i=0;i<RLXG_STATE_SIZE;i++)
	{
		fprintf(f,"%d ",state[i]);
	}
	fprintf(f,"\n");
	fclose(f);
	
	return LATAN_SUCCESS;
}

latan_errno randgen_load_state(randgen_state state, const stringbuf prefname)
{
	stringbuf fname,errmsg;
	FILE *f;
	int i;
	
	sprintf(fname,"%s.rand",prefname);
	FOPEN(f,fname,"r");
	for (i=0;i<RLXG_STATE_SIZE;i++)
	{
		if(fscanf(f,"%d ",state+i)<0)
		{
			sprintf(errmsg,"error while reading generator state component %d in file %s",\
					i,fname);
			LATAN_ERROR(errmsg,LATAN_ELATSYN);
		}
	}
	fclose(f);
	
	return LATAN_SUCCESS;
}
