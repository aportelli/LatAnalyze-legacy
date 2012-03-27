/* latan_plot.h, part of LatAnalyze library
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

#ifndef LATAN_PLOT_H_
#define LATAN_PLOT_H_

#include <latan/latan_globals.h>
#include <latan/latan_fit.h>
#include <latan/latan_mat.h>

/* system dependant gnuplot terminal setting */
#ifdef __APPLE__
#define DEFTERM "aqua"
#else
#define DEFTERM "x11"
#endif

__BEGIN_DECLS

/* flags */
enum
{
    PF_NOTHING = 0,
    PF_DATA    = 1 << 0,
    PF_FIT     = 1 << 1
};

/* function types */
typedef double univar_func(const double x, void *param);
typedef double mulvar_func(const mat *x, void *param);

/* the plot structure */
typedef struct 
{
    size_t nplot;       /* number of plot commands the buffer*/
    strbuf *plotbuf;    /* buffer of plot commands */
    size_t nhead;       /* number of head commands */
    strbuf *headbuf;    /* buffer of head commands */
    strbuf *tmpfname;   /* names of the temporary files opened for this plot */
    strbuf title;       /* title of the plot */
    strbuf term;        /* output terminal of the plot */
    strbuf output;      /* output file of the plot */
    unsigned int scale; /* scale mode of the plot */
    unsigned int log;   /* logarithmic mode of the axis */
    double xmin;        /* lower bound of x axis */
    double xmax;        /* upper bound of x axis */
    strbuf xlabel;      /* label of the x axis */
    double ymin;        /* lower bound of y axis */
    double ymax;        /* upper bound of y axis */
    strbuf ylabel;      /* label of the y axis */
} plot;

/* allocation */
plot *plot_create(void);
void plot_destroy(plot *p);

/* plot configuration */
void plot_set_scale_auto(plot *p);
void plot_set_scale_manual(plot *p, const double xmin, const double xmax,\
                           const double ymin, const double ymax);
void plot_set_scale_xmanual(plot *p, const double xmin, const double xmax);
void plot_set_scale_ymanual(plot *p, const double ymin, const double ymax);
void plot_set_scale_lin(plot *p);
void plot_set_scale_xlog(plot *p);
void plot_set_scale_ylog(plot *p);
void plot_set_scale_xylog(plot *p);
void plot_set_title(plot *p, const strbuf title);
void plot_set_xlabel(plot *p, const strbuf xlabel);
void plot_set_ylabel(plot *p, const strbuf xlabel);
void plot_set_term(plot *p, const strbuf term);
void plot_set_output(plot *p, const strbuf output);

/* plot functions */
void plot_add_plot(plot *p, const strbuf cmd);
void plot_add_head(plot *p, const strbuf cmd);
void plot_add_datpoint(plot *p, const double x, const double y,\
                       const double xerr, const double yerr,   \
                       const strbuf title, const strbuf color);
void plot_add_dat(plot *p, const mat *x, const mat *dat, const mat *xerr,\
                  const mat *yerr, const strbuf title, const strbuf color);
void plot_add_points(plot *p, const mat *x, const mat *y,const strbuf title,\
                     const strbuf color, const strbuf style);
#define plot_add_line(p,x,y,tit,col) plot_add_points(p,x,y,tit,col,"lines")
void plot_add_func(plot *p, univar_func *f, void *f_param, const double xmin,\
                   const double xmax, const size_t npt, const strbuf title,  \
                   const strbuf color);
void plot_add_parfunc(plot *p, mulvar_func *f, const mat *x_0, const size_t i,\
                      void *f_param, const double xmin, const double xmax,    \
                      const size_t npt, const strbuf title, const strbuf color);
void plot_add_model(plot *p, model_func *f, const mat *x_0, const size_t j, \
                    const mat *par, void *f_param, const double xmin,       \
                    const double xmax, const size_t npt, const strbuf title,\
                    const strbuf color);
void plot_add_histogram(plot *p, const mat *hist, const double xmin,\
                        const double xmax, const double w_tot,      \
                        const bool do_normalize, const strbuf title,\
                        const strbuf color);
#define plot_add_dat_xerr(p,x,dat,xerr,title,color)\
plot_add_dat(p,x,dat,xerr,NULL,title,color)
#define plot_add_dat_yerr(p,x,dat,yerr,title,color)\
plot_add_dat(p,x,dat,NULL,yerr,title,color)
#define plot_add_dat_xyerr plot_add_dat
void plot_add_hline(plot *p, const double y, const strbuf color);
void plot_add_hlineerr(plot *p, const double y, const double err,\
                       const strbuf color);
void plot_add_vline(plot *p, const double x, const strbuf color);
void plot_add_vlineaerr(plot *p, const double x, const double xerr[2],\
                        const strbuf color);
void plot_add_vlineerr(plot *p, const double x, const double xerr,\
                       const strbuf color);
void plot_add_fit(plot *p, const fit_data *d, const size_t ky, const mat *x_ex,\
                  const size_t kx, const mat *par, const double xmin,          \
                  const double xmax, const size_t npt, const bool do_sub,      \
                  const unsigned int obj,  const strbuf dat_title,             \
                  const strbuf fit_title, const strbuf dat_color,              \
                  const strbuf fit_color);

/* plot parsing */
void plot_parse(FILE* outstr, const plot *p);
#define plot_print(p) plot_parse(stdout,p);
latan_errno plot_disp(const plot *p);

__END_DECLS

#endif
