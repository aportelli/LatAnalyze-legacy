#include <latan/io.h>
#include <latan/includes.h>

int get_nrow(const stringbuf mark, const stringbuf matid,\
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
					LATAN_ERROR(errmsg,LATAN_ELATSYN);
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

int get_firstfname(stringbuf fname, const stringbuf manifestfname)
{
	stringbuf buf;
	FILE* manifest = NULL;
	
	FOPEN(manifest,manifestfname,"r");
	fgets(buf,STRING_LENGTH,manifest);
	fclose(manifest);
	sscanf(buf,"%s\n",fname);
	
	return LATAN_SUCCESS;
}

int get_mat(mat m, const stringbuf mark, const stringbuf matid,\
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
						
						sprintf(errmsg,"error reading matrix %s %s (%u,%u) in file %s (trying to read a %ux%u matrix)",\
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

int get_matn(mat* m, const stringbuf mark, const stringbuf matid,\
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
		LATAN_UPDATE_STATUS(status,get_mat(m[i],mark,matid,fname));
	}
	fclose(manifest);
	
	return status;
}

int save_plotdat(const mat dat, const double xstart, const double xstep,\
				 const stringbuf fname)
{
	FILE* f = NULL;
	size_t i;
	
	FOPEN(f,fname,"r");
	for (i=0;i<nrow(dat);i++)
	{
		fprintf(f,"%.10e %.10e\n",(double)(xstart+(int)(i)*xstep),\
				mat_get(dat,i,0));
	}
	fclose(f);
	
	return LATAN_SUCCESS;
}

int save_plotdaterr(const mat dat, const mat sig, const double xstart,\
					const double xstep, const stringbuf fname)
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