/* latan_tabfunc.c, part of LatAnalyze library
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
#include <latan/latan_tabfunc.h>
#include <latan/latan_includes.h>

double tabfunc_lineval(const double x, void* vf)
{
    size_t i;
    double dx,x_i,y_i,y_ip1;
    tabfunc *f;
    
    f = (tabfunc *)vf;
    
    if ((x>f->xmax)||(x<f->xmin))
    {
        LATAN_ERROR_VAL("variable outside of table range",LATAN_EDOM,GSL_NAN);
    }
    
    dx    = (f->xmax-f->xmin)/((double)(f->nval-1));
    i     = (size_t)floor((x-f->xmin)/dx);
    x_i   = f->xmin + ((double)(i))*dx;
    y_i   = f->table[i];
    y_ip1 = f->table[i+1];
    
    return y_i + (x-x_i)*(y_ip1-y_i)/dx;
}
