#include <latan/plot.h>
#include <latan/includes.h>

/*								internal code								*/
/****************************************************************************/
static char* gnuplot_get_program_path(const char* pname);
static void gnuplot_cmd(FILE *ctrl, const char* cmd, ...);
static void plot_add_tmpf(plot p, const stringbuf tmpfname);

/** part of Nicolas Devillard gnuplot interface **/
/*** Maximal size of a gnuplot command ***/
#define GP_CMD_SIZE     	2048
/*** Maximal size of a name in the PATH ***/
#define PATH_MAXNAMESZ       4096

#define GNUPLOT_OPEN(ctrl)\
{\
	if (getenv("DISPLAY") == NULL)\
	{\
		LATAN_ERROR("cannot find DISPLAY variable: is it set ?",LATAN_ESYSTEM);\
	}\
	if (gnuplot_get_program_path("gnuplot") == NULL)\
	{\
		LATAN_ERROR("cannot find gnuplot in your PATH",LATAN_ESYSTEM);\
	}\
	ctrl = popen("gnuplot -persist","w");\
	if (ctrl == NULL)\
	{\
		LATAN_ERROR("error starting gnuplot",LATAN_ESYSTEM);\
	}\
}

#define GNUPLOT_CLOSE(ctrl)\
{\
	if (pclose(ctrl) == -1)\
	{\
		LATAN_ERROR("problem closing communication to gnuplot",LATAN_ESYSTEM);\
	}\
}

static char* gnuplot_get_program_path(const char *pname)
{
    int         i, j, lg;
    char*		path;
    static char buf[PATH_MAXNAMESZ];
	
    /* Trivial case: try in CWD */
    sprintf(buf, "./%s", pname) ;
    if (access(buf, X_OK) == 0) {
        sprintf(buf, ".");
        return buf ;
    }
    /* Try out in all paths given in the PATH variable */
    buf[0] = 0;
    path = getenv("PATH") ;
    if (path!=NULL) {
        for (i=0; path[i]; ) {
            for (j=i ; (path[j]) && (path[j]!=':') ; j++);
            lg = j - i;
            strncpy(buf, path + i,(size_t)(lg));
            if (lg == 0) buf[lg++] = '.';
            buf[lg++] = '/';
            strcpy(buf + lg, pname);
            if (access(buf, X_OK) == 0) {
                /* Found it! */
                break ;
            }
            buf[0] = 0;
            i = j;
            if (path[i] == ':') i++ ;
        }
    } else {
		fprintf(stderr, "PATH variable not set\n");
	}
    /* If the buffer is still empty, the command was not found */
    if (buf[0] == 0) return NULL ;
    /* Otherwise truncate the command name to yield path only */
    lg = (int)(strlen(buf) - 1);
    while (buf[lg]!='/') {
        buf[lg] = 0;
        lg--;
    }
    buf[lg] = 0;
    return buf ;
}

static void gnuplot_cmd(FILE* ctrl, const char* cmd, ...)
{
    va_list ap;
    char local_cmd[GP_CMD_SIZE];
	
    va_start(ap,cmd);
    vsprintf(local_cmd,cmd,ap);
    va_end(ap);
    strcat(local_cmd,"\n");
    fputs(local_cmd,ctrl) ;
    fflush(ctrl);
}

/** temporary file management **/
static void plot_add_tmpf(plot p, const stringbuf tmpfname)
{
	(p->ntmpf)++;
	REALLOC_NOERRET(p->tmpfname,p->tmpfname,stringbuf*,p->ntmpf);
	strcpy(p->tmpfname[p->ntmpf-1],tmpfname);
}

/*								allocation									*/
/****************************************************************************/
plot plot_create(void)
{
	plot p;
	
	MALLOC_ERRVAL(p,plot,1,NULL);
	
	p->nplot = 0;
	p->plotbuf = NULL;
	p->ntmpf = 0;
	p->tmpfname = NULL;
	strcpy(p->title,"");
	strcpy(p->output,DEFTERM);
	p->scale = AUTO;
	p->log = NOLOG;
	p->xmin = 0.0;
	p->xmax = 0.0;
	p->ymin = 0.0;
	p->ymax = 0.0;
	
	return p;
}

void plot_destroy(plot p)
{
	size_t i;
	
	FREE(p->plotbuf);
	for (i=0;i<p->ntmpf;i++)
	{
		remove(p->tmpfname[i]);
	}
	FREE(p->tmpfname);
	FREE(p);
}

/*								plot configuration 							*/
/****************************************************************************/
void plot_set_scale_auto(plot p)
{
	p->scale = AUTO;
}

void plot_set_scale_manual(plot p, const double xmin, const double xmax,\
						   const double ymin, const double ymax)
{
	p->scale = MANUAL;
	p->xmin = xmin;
	p->xmax = xmax;
	p->ymin = ymin;
	p->ymax = ymax;
}

void plot_set_scale_xmanual(plot p, const double xmin, const double xmax)
{
	p->scale = XMANUAL;
	p->xmin = xmin;
	p->xmax = xmax;
}

void plot_set_scale_ymanual(plot p, const double ymin, const double ymax)
{
	p->scale = YMANUAL;
	p->ymin = ymin;
	p->ymax = ymax;
}

void plot_set_scale_lin(plot p)
{
	p->log = NOLOG;
}

void plot_set_scale_xlog(plot p)
{
	p->log = XLOG;
}

void plot_set_scale_ylog(plot p)
{
	p->log = YLOG;
}

void plot_set_scale_xylog(plot p)
{
	p->log = XYLOG;
}

void plot_set_title(plot p, const stringbuf title)
{
	strcpy(p->title,title);
}

/*								plot functions								*/
/****************************************************************************/
void plot_add_plot(plot p, const stringbuf cmd)
{
	(p->nplot)++;
	REALLOC_NOERRET(p->plotbuf,p->plotbuf,stringbuf*,p->nplot);
	strcpy(p->plotbuf[p->nplot-1],cmd);
}

void plot_add_dat(plot p, const mat dat, const double start, const double step,\
				  const stringbuf title)
{
	FILE* tmpf;
	stringbuf tmpfname, plotcmd;
	size_t i;
	
	strcpy(tmpfname,"latan_plot_tmp_XXXXXX");
	mkstemp(tmpfname);
	FOPEN_NOERRET(tmpf,tmpfname,"w");
	for (i=0;i<nrow(dat);i++)
	{
		fprintf(tmpf,"%.10e %.10e\n",start+(double)(i)*step,mat_get(dat,i,0));
	}
	fclose(tmpf);
	plot_add_tmpf(p,tmpfname);
	sprintf(plotcmd,"\"%s\" title \"%s\"",tmpfname,title);
	plot_add_plot(p,plotcmd);
}

void plot_add_daterr(plot p, const mat dat, const mat err, const double start,\
					const double step, const stringbuf title)
{
	FILE* tmpf;
	stringbuf tmpfname, plotcmd;
	size_t i;
	
	strcpy(tmpfname,"latan_plot_tmp_XXXXXX");
	mkstemp(tmpfname);
	FOPEN_NOERRET(tmpf,tmpfname,"w");
	for (i=0;i<nrow(dat);i++)
	{
		fprintf(tmpf,"%.10e %.10e %.10e\n",start+(double)(i)*step,\
				mat_get(dat,i,0),mat_get(err,i,0));
	}
	fclose(tmpf);
	plot_add_tmpf(p,tmpfname);
	sprintf(plotcmd,"\"%s\" title \"%s\"",tmpfname,title);
	plot_add_plot(p,plotcmd);
}

void plot_add_hline(plot p, const double y, const stringbuf style,\
					const stringbuf color)
{
	stringbuf plotcmd;

	sprintf(plotcmd,"%.10e lt %s lc %s notitle",y,style,color);
	plot_add_plot(p,plotcmd);
}

void plot_add_hlineerr(plot p, const double y, const double err,		\
					   const stringbuf style, const stringbuf color1,	\
					   const stringbuf color2)
{
	stringbuf plotcmd;
	
	sprintf(plotcmd,"%.10e lc %s with filledcurve y1=%.10e notitle",\
			y+err,color2,y-err);
	plot_add_plot(p,plotcmd);
	plot_add_hline(p,y,style,color1);
}

/*								plot parsing								*/
/****************************************************************************/
int plot_parse(FILE* outstr, const plot p)
{
	stringbuf begin, end;
	size_t i;
	
	gnuplot_cmd(outstr,"set term %s",p->output);
	switch (p->scale)
	{
		case AUTO:
			break;
		case MANUAL:
			gnuplot_cmd(outstr,"set xrange [%f:%f]",p->xmin,p->xmax);
			gnuplot_cmd(outstr,"set yrange [%f:%f]",p->ymin,p->ymax);
			break;
		case XMANUAL:
			gnuplot_cmd(outstr,"set xrange [%f:%f]",p->xmin,p->xmax);
			break;
		case YMANUAL:
			gnuplot_cmd(outstr,"set yrange [%f:%f]",p->ymin,p->ymax);
			break;
		default:
			LATAN_ERROR("plot scale mode unknow",LATAN_EINVAL);
			break;
	}
	switch (p->log)
	{
		case NOLOG:
			gnuplot_cmd(outstr,"set nolog xy");
			break;
		case XLOG:
			gnuplot_cmd(outstr,"set log x");
			break;
		case YLOG:
			gnuplot_cmd(outstr,"set log y");
			break;
		case XYLOG:
			gnuplot_cmd(outstr,"set log xy");
			break;
		default:
			LATAN_ERROR("plot log mode unknow",LATAN_EINVAL);
			break;
	}
	gnuplot_cmd(outstr,"set title \"%s\"",p->title);
	for (i=0;i<p->nplot;i++)
	{
		if (i == 0)
		{
			strcpy(begin,"plot ");
		}
		else 
		{
			strcpy(begin,"");
		}
		if (i == p->nplot-1)
		{
			strcpy(end,"");
		}
		else
		{
			strcpy(end,",\\");
		}
		gnuplot_cmd(outstr,"%s%s%s",begin,p->plotbuf[i],end);
	}
	
	return LATAN_SUCCESS;
}

int plot_disp(const plot p)
{
	int status;
	FILE* ctrl;
	
	GNUPLOT_OPEN(ctrl);
	
	status = plot_parse(ctrl,p);
	
	GNUPLOT_CLOSE(ctrl);
						
	return status;
}