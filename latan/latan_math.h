/* latan_math.h, part of LatAnalyze library
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

#ifndef LATAN_MATH_H_
#define LATAN_MATH_H_

#include <latan/latan_globals.h>
#include <latan/latan_mat.h>

#define SQ(x) ((x)*(x))
#define DRATIO(a,b) (((double)(a))/((double)(b)))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define C_PI 3.1415926535897932384626433832795028841970

__BEGIN_DECLS

/* binomial coefficients */
unsigned int binomial(const unsigned int n, const unsigned int p);

/* discrete central derivative */
latan_errno finite_diff(mat *ddat, const mat *dat);

/* row major order indexing */
size_t coord_to_rowmaj(const size_t *x, const size_t *dim, const size_t ndim);
void rowmaj_to_coord(size_t *x, const size_t *dim, const size_t ndim,\
                     const size_t ind);
__END_DECLS

#endif
