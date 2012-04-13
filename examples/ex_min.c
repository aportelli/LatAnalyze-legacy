/* ex_min.c, part of LatAnalyze library
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

#include <latan/latan_minimizer.h>

double polynom(const mat *var, void *param);

double polynom(const mat *var, void *param)
{
    double x;
    
    param = NULL;
    x = mat_get(var,0,0);
    
    return (x-1.0)*(x-2.0)*(x-3.0)*(x-4.1);
}

int main(void)
{
    mat *res;
    double fval;
    
    res = mat_create(1,1);
    
    mat_set(res,0,0,2.5);
    
    minimize(res,&fval,&polynom,NULL);
    
    printf("%f\n",mat_get(res,0,0));
    
    return 0;
}
