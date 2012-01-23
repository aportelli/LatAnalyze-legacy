/* latan_plot.c, part of LatAnalyze library
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

#define _POSIX_C_SOURCE 199209L

#include <latan/latan_plot.h>
#include <latan/latan_includes.h>

static char * gnuplot_get_program_path(const char *pname);
static void gnuplot_cmd(FILE *ctrl, const char *cmd, ...);
static void plot_add_tmpf(plot *p, const strbuf tmpfname);

static size_t ntmpf = 0;

/*                              internal code                               */
/****************************************************************************/

#ifndef GNUPLOT_CMD
#define GNUPLOT_CMD "gnuplot"
#endif
#ifndef GNUPLOT_CMD_ARGS
#define GNUPLOT_CMD_ARGS "-persist"
#endif

/** part of Nicolas Devillard gnuplot interface **/
/*** Maximal size of a gnuplot command ***/
#define GP_CMD_SIZE 2048
/*** Maximal size of a name in the PATH ***/
#define PATH_MAXNAMESZ 4096

#define GNUPLOT_OPEN(ctrl)\
{\
    strbuf cmd;\
    if (getenv("DISPLAY") == NULL)\
    {\
        LATAN_ERROR("cannot find DISPLAY variable: is it set ?",LATAN_ESYSTEM);\
    }\
    if (gnuplot_get_program_path(GNUPLOT_CMD) == NULL)\
    {\
        LATAN_ERROR("cannot find gnuplot in your PATH",LATAN_ESYSTEM);\
    }\
    sprintf(cmd,"%s %s",GNUPLOT_CMD,GNUPLOT_CMD_ARGS);\
    ctrl = popen(cmd,"w");\
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

static char * gnuplot_get_program_path(const char *pname)
{
    int         i, j, lg;
    char *      path;
    static char buf[PATH_MAXNAMESZ];
    
    /* Trivial case: try in CWD */
    sprintf(buf,"./%s", pname) ;
    if (access(buf,X_OK) == 0)
    {
        sprintf(buf,".");
        return buf ;
    }
    /* Try out in all paths given in the PATH variable */
    buf[0] = 0;
    path = getenv("PATH") ;
    if (path != NULL)
    {
        for (i=0;path[i];)
        {
            for (j=i;(path[j])&&(path[j]!=':');j++);
            lg = j - i;
            strncpy(buf,path + i,(size_t)(lg));
            if (lg == 0)
            {
                buf[lg++] = '.';
            }
            buf[lg++] = '/';
            strbufcpy(buf + lg, pname);
            if (access(buf, X_OK) == 0)
            {
                /* Found it! */
                break ;
            }
            buf[0] = 0;
            i = j;
            if (path[i] == ':') i++ ;
        }
    } 
    else
    {
        fprintf(stderr, "PATH variable not set\n");
    }
    /* If the buffer is still empty, the command was not found */
    if (buf[0] == 0)
    {
        return NULL;
    }
    /* Otherwise truncate the command name to yield path only */
    lg = (int)(strlen(buf) - 1);
    while (buf[lg]!='/')
    {
        buf[lg] = 0;
        lg--;
    }
    buf[lg] = 0;
    
    return buf ;
}

static void gnuplot_cmd(FILE* ctrl, const char *cmd, ...)
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
static void plot_add_tmpf(plot *p, const strbuf tmpfname)
{
    (ntmpf)++;
    REALLOC_NOERRET(p->tmpfname,p->tmpfname,strbuf*,ntmpf);
    strbufcpy(p->tmpfname[ntmpf-1],tmpfname);
}

/*                              allocation                                  */
/****************************************************************************/
plot *plot_create(void)
{
    plot *p;
    
    MALLOC_ERRVAL(p,plot *,1,NULL);
    
    p->nplot = 0;
    p->plotbuf = NULL;
    p->tmpfname = NULL;
    strbufcpy(p->title,"");
    strbufcpy(p->term,DEFTERM);
    strbufcpy(p->output,"");
    p->scale = AUTO;
    p->log = NOLOG;
    p->xmin = 0.0;
    p->xmax = 0.0;
    strbufcpy(p->xlabel,"");
    p->ymin = 0.0;
    p->ymax = 0.0;
    strbufcpy(p->ylabel,"");
    
    return p;
}

void plot_destroy(plot *p)
{
    size_t i;
    
    FREE(p->plotbuf);
    for (i=0;i<ntmpf;i++)
    {
        remove(p->tmpfname[i]);
    }
    FREE(p->tmpfname);
    FREE(p);
}

/*                              plot configuration                          */
/****************************************************************************/
void plot_set_scale_auto(plot *p)
{
    p->scale = AUTO;
}

void plot_set_scale_manual(plot *p, const double xmin, const double xmax,\
                           const double ymin, const double ymax)
{
    p->scale = MANUAL;
    p->xmin = xmin;
    p->xmax = xmax;
    p->ymin = ymin;
    p->ymax = ymax;
}

void plot_set_scale_xmanual(plot *p, const double xmin, const double xmax)
{
    p->scale = XMANUAL;
    p->xmin = xmin;
    p->xmax = xmax;
}

void plot_set_scale_ymanual(plot *p, const double ymin, const double ymax)
{
    p->scale = YMANUAL;
    p->ymin = ymin;
    p->ymax = ymax;
}

void plot_set_scale_lin(plot *p)
{
    p->log = NOLOG;
}

void plot_set_scale_xlog(plot *p)
{
    p->log = XLOG;
}

void plot_set_scale_ylog(plot *p)
{
    p->log = YLOG;
}

void plot_set_scale_xylog(plot *p)
{
    p->log = XYLOG;
}

void plot_set_title(plot *p, const strbuf title)
{
    strbufcpy(p->title,title);
}

void plot_set_xlabel(plot *p, const strbuf xlabel)
{
    strbufcpy(p->xlabel,xlabel);
}

void plot_set_ylabel(plot *p, const strbuf ylabel)
{
    strbufcpy(p->ylabel,ylabel);
}

void plot_set_term(plot *p, const strbuf term)
{
    strbufcpy(p->term,term);
}

void plot_set_output(plot *p, const strbuf output)
{
    strbufcpy(p->output,output);
}

/*                              plot functions                              */
/****************************************************************************/
void plot_add_plot(plot *p, const strbuf cmd)
{
    (p->nplot)++;
    REALLOC_NOERRET(p->plotbuf,p->plotbuf,strbuf*,p->nplot);
    strbufcpy(p->plotbuf[p->nplot-1],cmd);
}

void plot_add_dat(plot *p, const mat *x, const mat *dat, const strbuf title,\
                  const strbuf color)
{
    FILE* tmpf;
    strbuf tmpfname, plotcmd, colorcmd;
    size_t i;
    
    sprintf(tmpfname,"latan_plot_tmp_%lu",(long unsigned)ntmpf);
    FOPEN_NOERRET(tmpf,tmpfname,"w");
    for (i=0;i<nrow(dat);i++)
    {
        fprintf(tmpf,"%.10e %.10e\n",mat_get(x,i,0),mat_get(dat,i,0));
    }
    fclose(tmpf);
    plot_add_tmpf(p,tmpfname);
    if (strlen(color) == 0)
    {
        strbufcpy(colorcmd,"");
    }
    else
    {
        sprintf(colorcmd,"lc %s",color);
    }
    sprintf(plotcmd,"'%s' u 1:2 t '%s' lt -1 %s",tmpfname,title,colorcmd);
    plot_add_plot(p,plotcmd);
}

void plot_add_dat_yerr(plot *p, const mat *x, const mat *dat, const mat *yerr,\
                     const strbuf title, const strbuf color)
{
    FILE* tmpf;
    strbuf tmpfname, plotcmd, colorcmd;
    size_t i;
    
    sprintf(tmpfname,"latan_plot_tmp_%lu",(long unsigned)ntmpf);
    FOPEN_NOERRET(tmpf,tmpfname,"w");
    for (i=0;i<nrow(dat);i++)
    {
        fprintf(tmpf,"%.10e %.10e %.10e\n",mat_get(x,i,0),\
                mat_get(dat,i,0),mat_get(yerr,i,0));
    }
    fclose(tmpf);
    plot_add_tmpf(p,tmpfname);
    if (strlen(color) == 0)
    {
        strbufcpy(colorcmd,"");
    }
    else
    {
        sprintf(colorcmd,"lc %s",color);
    }
    sprintf(plotcmd,"'%s' u 1:2:3 w yerr t '%s' lt -1 %s",tmpfname,title,\
            colorcmd);
    plot_add_plot(p,plotcmd);
}

void plot_add_dat_xyerr(plot *p, const mat *x, const mat *dat, const mat *xerr,\
                        const mat *yerr, const strbuf title, const strbuf color)
{
    FILE* tmpf;
    strbuf tmpfname, plotcmd, colorcmd;
    size_t i;
    
    sprintf(tmpfname,"latan_plot_tmp_%lu",(long unsigned)ntmpf);
    FOPEN_NOERRET(tmpf,tmpfname,"w");
    for (i=0;i<nrow(dat);i++)
    {
        fprintf(tmpf,"%.10e %.10e %.10e %.10e\n",mat_get(x,i,0),\
                mat_get(dat,i,0),mat_get(xerr,i,0),mat_get(yerr,i,0));
    }
    fclose(tmpf);
    plot_add_tmpf(p,tmpfname);
    if (strlen(color) == 0)
    {
        strbufcpy(colorcmd,"");
    }
    else
    {
        sprintf(colorcmd,"lc %s",color);
    }
    sprintf(plotcmd,"'%s' u 1:2:3:4 w xyerr t '%s' lt -1 %s",tmpfname,title,\
            colorcmd);
    plot_add_plot(p,plotcmd);
}

void plot_add_hline(plot *p, const double y, const strbuf style,\
                    const strbuf color)
{
    strbuf plotcmd;

    sprintf(plotcmd,"%.10e lt %s lc %s notitle",y,style,color);
    plot_add_plot(p,plotcmd);
}

void plot_add_hlineerr(plot *p, const double y, const double err,       \
                       const strbuf style, const strbuf color1, \
                       const strbuf color2)
{
    strbuf plotcmd;
    
    sprintf(plotcmd,"%.10e lc %s with filledcurve y1=%.10e notitle",\
            y+err,color2,y-err);
    plot_add_plot(p,plotcmd);
    plot_add_hline(p,y,style,color1);
}

void plot_add_fit(plot *p, const fit_data *d, const size_t ky, const mat *x_ex,\
                  const size_t kx, const mat *par, const bool do_sub,          \
                  const unsigned int obj, const strbuf title,                  \
                  const strbuf style, const strbuf pcolor, const strbuf lcolor)
{
    mat *x,*x_err,*y,*y_err,*cor_data;
    bool have_x_err;
    size_t npt,ndata;
    size_t i,j;
    strbuf plotstr,plotcmd,lcolcmd,stcmd;
    
    npt        = fit_data_fit_point_num(d);
    ndata      = fit_data_get_ndata(d);
    j          = 0;
    have_x_err = fit_data_have_x_covar(d,kx);
    
    x          = mat_create(npt,1);
    y          = mat_create(npt,1);
    x_err      = mat_create(npt,1);
    y_err      = mat_create(npt,1);
    cor_data   = mat_create(ndata,1);

    if (par != NULL)
    {
        if (obj & PF_FIT)
        {
            fit_data_plot2dstr(plotstr,d,ky,x_ex,kx,par);
            if (strlen(lcolor) == 0)
            {
                strbufcpy(lcolcmd,"");
            }
            else
            {
                sprintf(lcolcmd,"lc %s",lcolor);
            }
            if (strlen(style) == 0)
            {
                strbufcpy(stcmd,"");
            }
            else
            {
                sprintf(stcmd,"lt %s",style);
            }
            sprintf(plotcmd,"%s notitle %s %s",plotstr,stcmd,lcolcmd);
            plot_add_plot(p,plotcmd);
        }
        if (do_sub)
        {
            fit_partresidual(cor_data,d,ky,x_ex,kx,par);
        }
        else
        {
            mat_cp(cor_data,fit_data_pt_y(d));
        }
    }
    else
    {
        fit_data_get_y_k(cor_data,d,ky);
    }
    if (obj & PF_DATA)
    {
        for (i=0;i<ndata;i++)
        {
            if (fit_data_is_fit_point(d,i))
            {
                mat_set(x,j,0,fit_data_get_x(d,i,kx));
                if (have_x_err)
                {
                    mat_set(x_err,j,0,                                       \
                            sqrt(mat_get(fit_data_pt_x_covar(d,kx,kx),i,i)));
                }
                mat_set(y,j,0,mat_get(cor_data,i,0));
                mat_set(y_err,j,0,                                       \
                        sqrt(mat_get(fit_data_pt_y_covar(d,ky,ky),i,i)));
                j++;
            }
        }
        if (have_x_err)
        {
            plot_add_dat_xyerr(p,x,y,x_err,y_err,title,pcolor);
        }
        else
        {
            plot_add_dat_yerr(p,x,y,y_err,title,pcolor);
        }
    }
    
    mat_destroy(x);
    mat_destroy(x_err);
    mat_destroy(y);
    mat_destroy(y_err);
    mat_destroy(cor_data);
}

/*                              plot parsing                                */
/****************************************************************************/
latan_errno plot_parse(FILE* outstr, const plot *p)
{
    strbuf begin, end;
    size_t i;
    
    gnuplot_cmd(outstr,"set term %s",p->term);
    if (strlen(p->output) != 0)
    {
        gnuplot_cmd(outstr,"set output '%s'",p->output);
    }
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
    gnuplot_cmd(outstr,"set title '%s'",p->title);
    gnuplot_cmd(outstr,"set xlabel '%s'",p->xlabel);
    gnuplot_cmd(outstr,"set ylabel '%s'",p->ylabel);
    for (i=0;i<p->nplot;i++)
    {
        if (i == 0)
        {
            strbufcpy(begin,"plot ");
        }
        else 
        {
            strbufcpy(begin,"");
        }
        if (i == p->nplot-1)
        {
            strbufcpy(end,"");
        }
        else
        {
            strbufcpy(end,",\\");
        }
        gnuplot_cmd(outstr,"%s%s%s",begin,p->plotbuf[i],end);
    }
    
    return LATAN_SUCCESS;
}

latan_errno plot_disp(const plot *p)
{
    latan_errno status;
    FILE* ctrl;
    
    GNUPLOT_OPEN(ctrl);
    
    status = plot_parse(ctrl,p);
    
    GNUPLOT_CLOSE(ctrl);
                        
    return status;
}
