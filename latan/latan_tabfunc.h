/* latan_tabfunc.h, part of LatAnalyze library
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
#ifndef LATAN_TABFUNC_H_
#define LATAN_TABFUNC_H_

#include <latan/latan_globals.h>

__BEGIN_DECLS

typedef struct tabfunc_s
{
    double xmin;
    double xmax;
    const double *table;
    size_t nval;
} tabfunc;

double tabfunc_lineval(const double x, void *vf);

__END_DECLS
               
#endif
