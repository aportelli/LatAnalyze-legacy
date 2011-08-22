/* latan_plot.h, part of LatAnalyze library
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

#ifndef LATAN_PLOT_H_
#define LATAN_PLOT_H_

#include <latan/latan_globals.h>
#include <latan/latan_fit.h>

/* system dependant gnuplot terminal setting */
#ifdef __APPLE__
#define DEFTERM "aqua"
#else
#define DEFTERM "x11"
#endif

/* flags */
#define AUTO 0
#define MANUAL 1
#define XMANUAL 2
#define YMANUAL 3

#define NOLOG 0
#define XLOG 1
#define YLOG 2
#define XYLOG 3

enum
{
    PF_NOTHING = 0,
    PF_DATA    = 1 << 0,
    PF_FIT     = 1 << 1
};

__BEGIN_DECLS

/* the plot structure */
typedef struct 
{
    size_t nplot;       /* number of plot commands the buffer*/
    strbuf *plotbuf;    /* buffer of plot commands */
    strbuf *tmpfname;   /* names of the temporary files opened for this plot */
    strbuf title;       /* title of the plot */
    strbuf term;        /* output terminal of the plot */
    strbuf output;      /* output file of the plot */
    int scale;          /* scale mode of the plot */
    int log;            /* logarithmic mode of the axis */
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
void plot_add_dat(plot *p, const mat *x, const mat *dat, const strbuf title,\
                  const strbuf color);
void plot_add_dat_yerr(plot *p, const mat *x, const mat *dat, const mat *yerr, \
                       const strbuf title, const strbuf color);
void plot_add_dat_xyerr(plot *p, const mat *x, const mat *dat, const mat *xerr,\
                        const mat *yerr, const strbuf title,                   \
                        const strbuf color);
void plot_add_hline(plot *p, const double y, const strbuf style,  \
                    const strbuf color);
void plot_add_hlineerr(plot *p, const double y, const double err, \
                       const strbuf style, const strbuf color1,   \
                       const strbuf color2);
void plot_add_fit(plot *p, const fit_data *d, const size_t k, const mat *x_ex,\
                  const mat *par, const bool do_sub, const unsigned int obj,  \
                  const strbuf title, const strbuf style, const strbuf pcolor,\
                  const strbuf lcolor);

/* plot parsing */
latan_errno plot_parse(FILE* outstr, const plot *p);
#define plot_print(p) plot_parse(stdout,p);
latan_errno plot_disp(const plot *p);

__END_DECLS

#endif
