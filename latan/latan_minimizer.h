/* latan_minimizer.h, part of LatAnalyze library
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

#ifndef LATAN_MINIMIZER_H_
#define LATAN_MINIMIZER_H_

#include <latan/latan_globals.h>
#include <latan/latan_mat.h>

__BEGIN_DECLS

/* minimize lib flags */
typedef enum
{
    GSL    = 0,\
    MINUIT = 1
} minlib_no;

/* minization algorithms */
#define NMINALG 6
typedef enum
{
    GSL_GRAD_FR    = 0,\
    GSL_GRAD_PR    = 1,\
    GSL_VEC_BFGS   = 2,\
    GSL_SIMPLEX_NM = 3,\
    MIN_MIGRAD     = 4,\
    MIN_SIMPLEX    = 5
} minalg_no;

minalg_no minalg_no_get(const strbuf m_id);
latan_errno minalg_id_get(strbuf m_id, const minalg_no n);

/* minimizer options */
minlib_no minimizer_get_lib(void);
minalg_no minimizer_get_alg(void);
latan_errno minimizer_set_alg(minalg_no alg);
latan_errno minimizer_get_alg_name(strbuf name);
unsigned int minimizer_get_max_iteration(void);
void minimizer_set_max_iteration(unsigned int max_iteration);

/* prototype of function to minimize */
typedef double min_func(const mat *x, void *param);

/* the minimizer */
latan_errno minimize(mat *x, const mat *x_limit, double *f_min, min_func *f,\
                     void *param);

__END_DECLS

#endif
