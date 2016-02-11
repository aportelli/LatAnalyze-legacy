/* ex_plot.c, part of LatAnalyze library
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

#include <stdlib.h>
#include <stdio.h>
#include <latan/latan_mat.h>
#include <latan/latan_plot.h>
#include <latan/latan_tabfunc.h>

const double sq_table[] =
{
    1.0,\
    4.0,\
    9.0,\
    16.0
};
tabfunc sq={1.0,4.0,sq_table,4};

int main(void)
{
    plot *p;
    double xerr[2];
    
    p = plot_create();

    xerr[0] = 0.2;
    xerr[1] = 0.5;
    
    plot_set_scale_ymanual(p,0.0,16.0);
    plot_add_vlineaerr(p,3,xerr,"rgb 'green'");
    plot_add_plot(p,"x**2 lc rgb 'red'","");
    plot_add_func(p,&tabfunc_lineval,&sq,sq.xmin,sq.xmax,1000,"x**2 table",\
                  "rgb 'blue'");
    plot_save("ex_plot_save",p);
    plot_disp(p);
    plot_print(p);
    
    
    plot_destroy(p);
    
    return EXIT_SUCCESS;
}
