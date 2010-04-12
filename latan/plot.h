/*!
 @file plot.h
 @brief <A HREF=http://www.gnuplot.info/>gnuplot</A> interface based on Nicolas Devillard <A HREF=http://ndevilla.free.fr/gnuplot/>gnuplot C API</A>.
 @author <A HREF=mailto:antonin.portelli@gmail.com>Antonin Portelli</A> (<A HREF=http://www.cpt.univ-mrs.fr/>CPT</A>)
*/

/*!
 @example ex_plot.c Example use of gnuplot interface.
*/

#ifndef LATAN_PLOT_H_
#define LATAN_PLOT_H_

#include <latan/globals.h>
#include <latan/mat.h>

/* system dependant gnuplot terminal setting */
/*!
 @def DEFTERM
 @brief Default gnuplot terminal. Automatically defined on \c "aqua" when compiled on a Mac OS X system.
*/
#define DEFTERM "x11"
#ifdef __APPLE__
#undef DEFTERM
#define DEFTERM "aqua"
#endif

/* flags */
/*!
 @def AUTO
 @brief Flag for plot::scale.
*/
#define AUTO 0
/*!
 @def MANUAL
 @brief Flag for plot::scale.
 */
#define MANUAL 1
#define XMANUAL 2
#define YMANUAL 3

#define NOLOG 0
#define XLOG 1
#define YLOG 2
#define XYLOG 3

__BEGIN_DECLS

/* the plot structure */
/*!
 @struct plot
 @brief Structure representing a gnuplot %plot.
 
 @remark It is not advised to modify directly the fields of this structure, use functions and macros provided in plot.h 
 to manipulate plot variables.
 @warning Don't forget to initialize (with #PLOT_INIT) plot variables before using them. Don't forget also to free them 
 (with #PLOT_FREE) after use to avoid memory leaks.
*/
typedef struct 
{
	size_t nplot;			/*!< number of %plot commands the buffer*/
	stringbuf* plotbuf;		/*!< buffer of %plot commands (see ::plot_disp function for more information) */
	size_t ntmpf;			/*!< number of temporary files opened for this %plot */
	stringbuf* tmpfname;	/*!< names of the temporary files opened for this %plot */
	stringbuf title;		/*!< title of the %plot */
	stringbuf output;		/*!< output terminal of the %plot */
	int scale;				/*!< scale mode of the %plot (see #PLOT_SET_MANSCALE and #PLOT_SET_AUTOSCALE)*/
	int log;				/*!< logarithmic mode of the axis */
	double xmin;			/*!< lower bound of x axis (see #PLOT_SET_MANSCALE) */
	double xmax;			/*!< upper bound of x axis (see #PLOT_SET_MANSCALE) */
	double ymin;			/*!< lower bound of y axis (see #PLOT_SET_MANSCALE) */
	double ymax;			/*!< upper bound of y axis (see #PLOT_SET_MANSCALE) */
}* plot;

/* allocation */
plot plot_create(void);
void plot_destroy(plot p);

/* plot configuration */
void plot_set_scale_auto(plot p);
void plot_set_scale_manual(plot p, const double xmin, const double xmax,\
						   const double ymin, const double ymax);
void plot_set_scale_xmanual(plot p, const double xmin, const double xmax);
void plot_set_scale_ymanual(plot p, const double ymin, const double ymax);
void plot_set_scale_lin(plot p);
void plot_set_scale_xlog(plot p);
void plot_set_scale_ylog(plot p);
void plot_set_scale_xylog(plot p);
void plot_set_title(plot p, const stringbuf title);

/* plot functions */
void plot_add_plot(plot p, const stringbuf cmd);
void plot_add_dat(plot p, const mat dat, const double start, const double step,\
				  const stringbuf title);
void plot_add_daterr(plot p, const mat dat, const mat err, const double start,\
					 const double step, const stringbuf title);
void plot_add_hline(plot p, const double y, const stringbuf style,\
					const stringbuf color);
void plot_add_hlineerr(plot p, const double y, const double err,		\
					   const stringbuf style, const stringbuf color1,	\
					   const stringbuf color2);

/* plot parsing */
latan_errno plot_parse(FILE* outstr, const plot p);
#define plot_print(p) plot_parse(stdout,p);
/*!
 @fn void plot_disp(const plot p)
 @brief Tell gnuplot to print the plot information contained in \b p on terminal \b p.output. 
 
 Plot command executed by gnuplot will be formated in the folowing way :<BR><BR>
 <TT>%plot&nbsp;</TT>\b p.plotbuf[0]<TT>,\ </TT><BR>
 <TT>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</TT>\b p.plotbuf[1]<TT>,\ </TT><BR>
 <TT>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;.</TT><BR>
 <TT>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;.</TT><BR>
 <TT>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;.</TT><BR>
 <TT>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</TT>\b p.plotbuf[nplot-1]<BR><BR>
 
 @param p plot to print
*/
latan_errno plot_disp(const plot p);

__END_DECLS

#endif