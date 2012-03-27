/* latan_mass.h, part of LatAnalyze library
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

#ifndef LATAN_MASS_H_
#define LATAN_MASS_H_

#include <latan/latan_globals.h>
#include <latan/latan_fit.h>
#include <latan/latan_mat.h>
#include <latan/latan_statistics.h>

__BEGIN_DECLS

typedef struct
{
    size_t start;
    size_t end;
    double mean;
    double sig;
} plat;

/* effective mass functions */
latan_errno effmass(mat *res, const mat *mprop, const int parity);
latan_errno rs_sample_effmass(rs_sample *s_res, const rs_sample *s_mprop,\
                              const int parity);
latan_errno effmass_PCAC(mat *res, const mat *mprop_AP, const mat *mprop_PP);
#define rs_sample_effmass_PCAC(s_res,s_mprop_AP,s_mprop_PP)\
rs_sample_binop(s_res,s_mprop_AP,s_mprop_PP,&effmass_PCAC)

/* interface to read experimental masses in latan_nunits.h */
latan_errno get_mass(double mass[2], const strbuf name);

/* mass fit parameter tuning functions */
plat *search_plat(size_t *nplat, mat *data, mat *sigdata,\
                  const size_t ntmax, const double nsig, const double tol);
latan_errno fit_data_mass_fit_tune(fit_data *d, mat *fit_init, mat *prop,\
                                   mat *em, mat *sigem, const int parity);

__END_DECLS

#endif
