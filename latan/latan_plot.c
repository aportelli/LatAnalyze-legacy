/* latan_plot.c, part of LatAnalyze library
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

#define _POSIX_C_SOURCE 199209L

#include <latan/latan_plot.h>
#include <latan/latan_includes.h>
#include <latan/latan_io.h>
#include <latan/latan_math.h>

enum
{
    AUTO    = 0,
    XMANUAL = 1 << 0,
    YMANUAL = 1 << 1
};

enum
{
    NOLOG = 0,
    XLOG  = 1 << 0,
    YLOG  = 1 << 1
};

static char * gnuplot_get_program_path(const char *pname);
static void gnuplot_cmd(FILE *ctrl, const char *cmd, ...);

static size_t ntmpf = 0;

/*                              internal code                               */
/****************************************************************************/

#ifndef GNUPLOT_CMD
#define GNUPLOT_CMD "gnuplot"
#endif
#ifndef GNUPLOT_CMD_ARGS
#define GNUPLOT_CMD_ARGS "-p"
#endif

/** part of Nicolas Devillard gnuplot interface **/
/*** Maximal size of a gnuplot command ***/
#define GP_CMD_SIZE STRING_LENGTH
/*** Maximal size of a name in the PATH ***/
#define PATH_MAXNAMESZ 4096

#define GNUPLOT_OPEN(ctrl)\
{\
    strbuf cmd;\
    if (getenv("DISPLAY") == NULL)\
    {\
        LATAN_ERROR("cannot find DISPLAY variable: is it set ?",LATAN_ESYSTEM);\
    }\
    gnuplot_get_program_path(GNUPLOT_CMD);\
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
        LATAN_ERROR_VAL("PATH variable not set",LATAN_ESYSTEM,NULL);
    }
    /* If the buffer is still empty, the command was not found */
    if (buf[0] == 0)
    {
        LATAN_ERROR_VAL("cannot find gnuplot in your PATH",LATAN_ESYSTEM,NULL);
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

/*                              allocation                                  */
/****************************************************************************/
plot *plot_create(void)
{
    plot *p;
    
    MALLOC_ERRVAL(p,plot *,1,NULL);
    
    p->nplot    = 0;
    p->plotbuf  = NULL;
    p->nhead    = 0;
    p->headbuf  = NULL;
    p->datfname = NULL;
    strbufcpy(p->title,"");
    strbufcpy(p->term,DEFTERM);
    strbufcpy(p->output,"");
    p->scale = AUTO;
    p->log   = NOLOG;
    p->xmin  = 0.0;
    p->xmax  = 0.0;
    strbufcpy(p->xlabel,"");
    p->ymin  = 0.0;
    p->ymax  = 0.0;
    strbufcpy(p->ylabel,"");
    
    return p;
}

void plot_destroy(plot *p)
{
    size_t i;
    strbuf war_msg;
    
    if (p)
    {
        FREE(p->headbuf);
        FREE(p->plotbuf);
        for (i=0;i<p->nplot;i++)
        {
            if (strlen(p->datfname[i]))
            {
                if (remove(p->datfname[i]))
                {
                    sprintf(war_msg,"impossible to remove plot temporary file %s",\
                            p->datfname[i]);
                    LATAN_WARNING(war_msg,LATAN_ESYSTEM);
                }
            }
        }
        FREE(p->datfname);
        FREE(p);
    }
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
    p->scale = XMANUAL|YMANUAL;
    p->xmin  = MIN(xmin,xmax);
    p->xmax  = MAX(xmin,xmax);
    p->ymin  = MIN(ymin,ymax);
    p->ymax  = MAX(ymin,ymax);
}

void plot_set_scale_xmanual(plot *p, const double xmin, const double xmax)
{
    p->scale |= XMANUAL;
    p->xmin   = MIN(xmin,xmax);
    p->xmax   = MAX(xmin,xmax);
}

void plot_set_scale_ymanual(plot *p, const double ymin, const double ymax)
{
    p->scale |= YMANUAL;
    p->ymin   = MIN(ymin,ymax);
    p->ymax   = MAX(ymin,ymax);
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
    p->log = XLOG|YLOG;
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
void plot_add_plot(plot *p, const strbuf cmd, const strbuf datfname)
{
    strbuf war_msg;
    
    (p->nplot)++;
    REALLOC_NOERRET(p->plotbuf,p->plotbuf,strbuf *,p->nplot);
    REALLOC_NOERRET(p->datfname,p->datfname,strbuf *,p->nplot);
    strbufcpy(p->plotbuf[p->nplot-1],cmd);
    strbufcpy(p->datfname[p->nplot-1],datfname);
    if (strlen(datfname)&&(access(datfname,R_OK) != 0))
    {
        sprintf(war_msg,"cannot access datafile %s",datfname);
        LATAN_WARNING(war_msg,LATAN_ESYSTEM);
    }
}

void plot_add_head(plot *p, const strbuf cmd)
{
    (p->nhead)++;
    REALLOC_NOERRET(p->headbuf,p->headbuf,strbuf*,p->nhead);
    strbufcpy(p->headbuf[p->nhead-1],cmd);
}

enum
{
    NO_ERR = 0,
    X_ERR  = 1 << 0,
    Y_ERR  = 1 << 1
};

void plot_add_datpoint(plot *p, const double x, const double y, \
                       const double xerr,const double yerr,     \
                       const strbuf title, const strbuf color)
{
    strbuf ucmd, echocmd, errcmd, plotcmd, colorcmd;
    unsigned int err_flag;
    
    err_flag = NO_ERR;
    
    err_flag |= (xerr > 0.0) ? X_ERR : NO_ERR;
    err_flag |= (yerr > 0.0) ? Y_ERR : NO_ERR;
    if ((err_flag & X_ERR)&&(err_flag & Y_ERR))
    {
        strbufcpy(ucmd,"1:2:3:4");
        strbufcpy(errcmd,"w xyerr");
        sprintf(echocmd,"'< echo %e %e %e %e'",x,y,xerr,yerr);
    }
    else if (err_flag & X_ERR)
    {
        strbufcpy(ucmd,"1:2:3");
        strbufcpy(errcmd,"w xerr");
        sprintf(echocmd,"'< echo %e %e %e'",x,y,xerr);
    }
    else if (err_flag & Y_ERR)
    {
        strbufcpy(ucmd,"1:2:3");
        strbufcpy(errcmd,"w yerr");
        sprintf(echocmd,"'< echo %e %e %e'",x,y,yerr);
    }
    else
    {
        strbufcpy(ucmd,"1:2");
        strbufcpy(errcmd,"");
        sprintf(echocmd,"'< echo %e %e'",x,y);
    }
    if (strlen(color) == 0)
    {
        strbufcpy(colorcmd,"");
    }
    else
    {
        sprintf(colorcmd,"lc %s",color);
    }
    sprintf(plotcmd,"%s u %s %s t '%s' lt -1 pt 1 %s",echocmd,ucmd,errcmd,\
            title,colorcmd);
    plot_add_plot(p,plotcmd,"");
}

void plot_add_dat(plot *p, const mat *x, const mat *dat, const mat *xerr,\
                  const mat *yerr, const strbuf title, const strbuf color)
{
    FILE* tmpf;
    strbuf tmpfname, ucmd, errcmd, plotcmd, colorcmd;
    unsigned int err_flag;
    size_t i;
    
    err_flag = NO_ERR;
    
    err_flag |= (xerr != NULL) ? X_ERR : NO_ERR;
    err_flag |= (yerr != NULL) ? Y_ERR : NO_ERR;
    sprintf(tmpfname,".latan_tmp_plot_%lu.dat",(long unsigned)ntmpf);
    ntmpf++;
    FOPEN_NOERRET(tmpf,tmpfname,"w");
    if ((err_flag & X_ERR)&&(err_flag & Y_ERR))
    {
        strbufcpy(ucmd,"1:2:3:4");
        strbufcpy(errcmd,"w xyerr");
        for (i=0;i<nrow(dat);i++)
        {
            fprintf(tmpf,"%.10e %.10e %.10e %.10e\n",mat_get(x,i,0),\
                    mat_get(dat,i,0),mat_get(xerr,i,0),             \
                    mat_get(yerr,i,0));
        }
    }
    else if (err_flag & X_ERR)
    {
        strbufcpy(ucmd,"1:2:3");
        strbufcpy(errcmd,"w xerr");
        for (i=0;i<nrow(dat);i++)
        {
            fprintf(tmpf,"%.10e %.10e %.10e\n",mat_get(x,i,0),mat_get(dat,i,0),\
                    mat_get(xerr,i,0));
        }
    }
    else if (err_flag & Y_ERR)
    {
        strbufcpy(ucmd,"1:2:3");
        strbufcpy(errcmd,"w yerr");
        for (i=0;i<nrow(dat);i++)
        {
            fprintf(tmpf,"%.10e %.10e %.10e\n",mat_get(x,i,0),mat_get(dat,i,0),\
                    mat_get(yerr,i,0));
        }
    }
    else
    {
        strbufcpy(ucmd,"1:2");
        strbufcpy(errcmd,"");
        for (i=0;i<nrow(dat);i++)
        {
            fprintf(tmpf,"%.10e %.10e\n",mat_get(x,i,0),mat_get(dat,i,0));
        }
    }
    fclose(tmpf);
    if (strlen(color) == 0)
    {
        strbufcpy(colorcmd,"");
    }
    else
    {
        sprintf(colorcmd,"lc %s",color);
    }
    sprintf(plotcmd,"u %s %s t '%s' lt -1 pt 1 %s",ucmd,errcmd,title,\
            colorcmd);
    plot_add_plot(p,plotcmd,tmpfname);
}

void plot_add_points(plot *p, const mat *x, const mat *y, const strbuf title,\
                     const strbuf color, const strbuf style)
{
    FILE* tmpf;
    strbuf tmpfname, plotcmd, colorcmd;
    size_t i;
    
    sprintf(tmpfname,".latan_tmp_plot_%lu.dat",(long unsigned)ntmpf);
    ntmpf++;
    FOPEN_NOERRET(tmpf,tmpfname,"w");
    for (i=0;i<nrow(y);i++)
    {
        fprintf(tmpf,"%.10e %.10e\n",mat_get(x,i,0),mat_get(y,i,0));
    }
    fclose(tmpf);
    if (strlen(color) == 0)
    {
        strbufcpy(colorcmd,"");
    }
    else
    {
        sprintf(colorcmd,"lc %s",color);
    }
    sprintf(plotcmd,"u 1:2 t '%s' lt -1 %s w %s",title,\
            colorcmd,style);
    plot_add_plot(p,plotcmd,tmpfname);
}

void plot_add_func(plot *p, univar_func *f, void *f_param, const double xmin,\
                   const double xmax, const size_t npt, const strbuf title,  \
                   const strbuf color)
{
    mat *x,*y;
    size_t i;
    double x_i,y_i;
    
    x = mat_create(npt,1);
    y = mat_create(npt,1);
    
    for (i=0;i<npt;i++)
    {
        x_i = xmin + (xmax-xmin)*DRATIO(i,npt-1);
        y_i = f(x_i,f_param);
        mat_set(x,i,0,x_i);
        mat_set(y,i,0,y_i);
    }
    plot_add_line(p,x,y,title,color);
    
    mat_destroy(x);
    mat_destroy(y);
}

void plot_add_parfunc(plot *p, mulvar_func *f, const mat *x_0, const size_t j,\
                      void *f_param, const double xmin, const double xmax,    \
                      const size_t npt, const strbuf title, const strbuf color)
{
    mat *x,*y,*X;
    size_t i;
    double x_i,y_i;
    
    x = mat_create(npt,1);
    y = mat_create(npt,1);
    X = mat_create_from_mat(x_0);
    
    for (i=0;i<npt;i++)
    {
        x_i = xmin + (xmax-xmin)*DRATIO(i,npt-1);
        mat_set(X,j,0,x_i);
        y_i = f(X,f_param);
        mat_set(x,i,0,x_i);
        mat_set(y,i,0,y_i);
    }
    plot_add_line(p,x,y,title,color);
    
    mat_destroy(x);
    mat_destroy(y);
    mat_destroy(X);
}

void plot_add_model(plot *p, model_func *f, const mat *x_0, const size_t j, \
                    const mat *par, void *f_param, const double xmin,       \
                    const double xmax, const size_t npt, const strbuf title,\
                    const strbuf color)
{
    mat *x,*y,*X;
    size_t i;
    double x_i,y_i;
    
    x = mat_create(npt,1);
    y = mat_create(npt,1);
    X = mat_create_from_mat(x_0);
    
    for (i=0;i<npt;i++)
    {
        x_i = xmin + (xmax-xmin)*DRATIO(i,npt-1);
        mat_set(X,j,0,x_i);
        y_i = f(X,par,f_param);
        mat_set(x,i,0,x_i);
        mat_set(y,i,0,y_i);
    }
    plot_add_line(p,x,y,title,color);
    
    mat_destroy(x);
    mat_destroy(y);
    mat_destroy(X);
}

void plot_add_histogram(plot *p, const mat *hist, const double xmin,\
                        const double xmax, const double w_tot,      \
                        const bool do_normalize, const strbuf title,\
                        const strbuf color)
{
    mat *nhist,*x;
    double dnbin;
    
    nhist = mat_create_from_mat(hist);
    x     = mat_create_from_dim(hist);
    
    dnbin = (double)nrow(hist);
    
    mat_set_step(x,xmin,(xmax-xmin)/dnbin);
    if (do_normalize)
    {
        mat_eqmuls(nhist,dnbin/(w_tot*(xmax-xmin)));
    }
    plot_add_points(p,x,nhist,title,color,"steps");
    
    mat_destroy(nhist);
    mat_destroy(x);
}

void plot_add_hline(plot *p, const double y, const strbuf color)
{
    strbuf plotcmd,colorcmd;

    if (strlen(color) == 0)
    {
        strbufcpy(colorcmd,"");
    }
    else
    {
        sprintf(colorcmd,"lc %s",color);
    }
    sprintf(plotcmd,"%.10e lt -1 %s notitle",y,colorcmd);
    plot_add_plot(p,plotcmd,"");
}

void plot_add_hlineerr(plot *p, const double y, const double err,\
                       const strbuf color)
{
    strbuf plotcmd,colorcmd;
    
    if (strlen(color) == 0)
    {
        strbufcpy(colorcmd,"");
    }
    else
    {
        sprintf(colorcmd,"lc %s",color);
    }
    sprintf(plotcmd,"%.10e %s w filledcurve y1=%.10e fs transparent solid 0.25 noborder notitle",\
            y+err,colorcmd,y-err);
    plot_add_plot(p,plotcmd,"");
    plot_add_hline(p,y,color);
}

void plot_add_vline(plot *p, const double x, const strbuf color)
{
    strbuf plotcmd,colorcmd;
    
    if (!(p->scale & YMANUAL))
    {
        LATAN_WARNING("vertical line will not work without manual y-scale setting",\
                      LATAN_EINVAL);
    }
    
    if (strlen(color) == 0)
    {
        strbufcpy(colorcmd,"");
    }
    else
    {
        sprintf(colorcmd,"lc %s",color);
    }
    sprintf(plotcmd,"sprintf(\"< printf ' %e %%e\\n%e %%e\\n'\",ymin,ymax) lt -1 %s w lines notitle",\
            x,x,colorcmd);
    plot_add_plot(p,plotcmd,"");
}

void plot_add_vlineaerr(plot *p, const double x, const double xerr[2],\
                        const strbuf color)
{
    strbuf plotcmd,colorcmd;
    
    if (strlen(color) == 0)
    {
        strbufcpy(colorcmd,"");
    }
    else
    {
        sprintf(colorcmd,"lc %s",color);
    }
    sprintf(plotcmd,"sprintf(\"< printf ' %e %%e\\n%e %%e\\n%e %%e\\n%e %%e\\n%e %%e\\n'\",ymin,ymax,ymax,ymin,ymin) %s w filledcurve fs transparent solid 0.25 noborder notitle",x-xerr[0],x-xerr[0],x+xerr[1],x+xerr[1],x-xerr[0],colorcmd);
    plot_add_plot(p,plotcmd,"");
    plot_add_vline(p,x,color);
}

void plot_add_vlineerr(plot *p, const double x, const double xerr,\
                       const strbuf color)
{
    double xerr_a[2];
    
    xerr_a[0] = xerr;
    xerr_a[1] = xerr;
    
    plot_add_vlineaerr(p,x,xerr_a,color);
}

void plot_add_fit(plot *p, const fit_data *d, const size_t ky, const mat *x_ex,\
                  const size_t kx, const mat *par, const double xmin,          \
                  const double xmax, const size_t npt, const bool do_sub,      \
                  const unsigned int obj,  const strbuf dat_title,             \
                  const strbuf fit_title, const strbuf dat_color,              \
                  const strbuf fit_color)
{
    mat *x,*x_err,*y,*y_err,*cor_data;
    bool have_x_err;
    size_t nfitpt,ndata;
    size_t i,j;
    
    nfitpt        = fit_data_fit_point_num(d);
    ndata      = fit_data_get_ndata(d);
    j          = 0;
    have_x_err = fit_data_have_x_covar(d,kx);
    
    x          = mat_create(nfitpt,1);
    y          = mat_create(nfitpt,1);
    x_err      = mat_create(nfitpt,1);
    y_err      = mat_create(nfitpt,1);
    cor_data   = mat_create(ndata,1);

    if (par != NULL)
    {
        if (obj & PF_FIT)
        {
            plot_add_model(p,d->model->func[ky],x_ex,kx,par,d->model_param,\
                           xmin,xmax,npt,fit_title,fit_color);
        }
        if ((obj & PF_DATA)&&(do_sub))
        {
            fit_partresidual(cor_data,d,ky,x_ex,kx,par);
        }
        else
        {
            fit_data_get_y_k(cor_data,d,ky);
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
            plot_add_dat(p,x,y,x_err,y_err,dat_title,dat_color);
        }
        else
        {
            plot_add_dat(p,x,y,NULL,y_err,dat_title,dat_color);
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
void plot_parse(FILE* outstr, const plot *p)
{
    strbuf begin, end, fname;
    size_t i;
    
    gnuplot_cmd(outstr,"set term %s",p->term);
    if (strlen(p->output) != 0)
    {
        gnuplot_cmd(outstr,"set output '%s'",p->output);
    }
    gnuplot_cmd(outstr,"xmin=%e",p->xmin);
    gnuplot_cmd(outstr,"xmax=%e",p->xmax);
    gnuplot_cmd(outstr,"ymin=%e",p->ymin);
    gnuplot_cmd(outstr,"ymax=%e",p->ymax);
    if (p->scale & XMANUAL)
    {
        gnuplot_cmd(outstr,"set xrange [xmin:xmax]");
    }
    if (p->scale & YMANUAL)
    {
        gnuplot_cmd(outstr,"set yrange [ymin:ymax]");
    }
    if (p->log & XLOG)
    {
        gnuplot_cmd(outstr,"set log x");
    }
    if (p->log & YLOG)
    {
        gnuplot_cmd(outstr,"set log y");
    }
    gnuplot_cmd(outstr,"set title '%s'",p->title);
    gnuplot_cmd(outstr,"set xlabel '%s'",p->xlabel);
    gnuplot_cmd(outstr,"set ylabel '%s'",p->ylabel);
    for (i=0;i<p->nhead;i++)
    {
        gnuplot_cmd(outstr,p->headbuf[i]);
    }
    for (i=0;i<p->nplot;i++)
    {
        if (i == 0)
        {
            strbufcpy(begin,"plot ");
        }
        else 
        {
            strbufcpy(begin,"     ");
        }
        if (i == p->nplot-1)
        {
            strbufcpy(end,"");
        }
        else
        {
            strbufcpy(end,",\\");
        }
        if (strlen(p->datfname[i]))
        {
            sprintf(fname,"'%s' ",p->datfname[i]);
        }
        else
        {
            strbufcpy(fname,"");
        }
        gnuplot_cmd(outstr,"%s%s%s%s",begin,fname,p->plotbuf[i],end);
    }
}

latan_errno plot_disp(const plot *p)
{
    FILE* ctrl;
    
    GNUPLOT_OPEN(ctrl);
    
    plot_parse(ctrl,p);
    
    GNUPLOT_CLOSE(ctrl);
                        
    return LATAN_SUCCESS;
}

latan_errno plot_save(const strbuf dirname, plot *p)
{
    strbuf path,buf,*datfname_bak,term_bak,output_bak;
    FILE *outf;
    size_t i,j,lc;
    mode_t mode_755;
    strbuf name,ver,err_msg;
    
    MALLOC_ERRVAL(datfname_bak,strbuf *,p->nplot,LATAN_ENOMEM);
    
    mode_755 = S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH;
    
    /* backup I/O parameters */
    strbufcpy(term_bak,p->term);
    strbufcpy(output_bak,p->output);
    for (i=0;i<p->nplot;i++)
    {
        strbufcpy(datfname_bak[i],p->datfname[i]);
    }
    /* generate directory */
    if (access(dirname,R_OK|W_OK|X_OK))
    {
        if (mkdir(dirname,mode_755))
        {
            sprintf(err_msg,"impossible to create directory %s",dirname);
            LATAN_ERROR(err_msg,LATAN_ESYSTEM);
        }
    }
    /* save pdf */
    sprintf(path,"%s/plot.pdf",dirname);
    plot_set_term(p,"pdf");
    plot_set_output(p,path);
    plot_disp(p);
    plot_set_term(p,term_bak);
    plot_set_output(p,output_bak);
    /* save script and datafiles */
    j = 0;
    for (i=0;i<p->nplot;i++)
    {
        if (strlen(p->datfname[i]))
        {
            sprintf(p->datfname[i],"points_%lu.dat",(long unsigned)(j));
            j++;
            sprintf(path,"%s/%s",dirname,p->datfname[i]);
            FOPEN(outf,path,"w");
            BEGIN_FOR_LINE(buf,datfname_bak[i],lc)
            {
                fprintf(outf,"%s\n",buf);
            }
            END_FOR_LINE
            fclose(outf);
        }
    }
    sprintf(path,"%s/disp.plt",dirname);
    FOPEN(outf,path,"w");
    fprintf(outf,"#!/usr/bin/env %s %s\n\n",GNUPLOT_CMD,GNUPLOT_CMD_ARGS);
    latan_get_name(name);
    latan_get_version(ver);
    fprintf(outf,"# script generated by %s v%s\n\n",name,ver);
    plot_parse(outf,p);
    fprintf(outf,"\n");
    fclose(outf);
    if(chmod(path,mode_755))
    {
        sprintf(err_msg,"impossible to set file %s in mode 755",path);
        LATAN_WARNING(err_msg,LATAN_ESYSTEM);
    }
    for (i=0;i<p->nplot;i++)
    {
        strbufcpy(p->datfname[i],datfname_bak[i]);
    }
    
    return LATAN_SUCCESS;
}
