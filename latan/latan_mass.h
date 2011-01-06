/* latan_mass.h, part of LatAnalyze library
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

#ifndef LATAN_MASS_H_
#define LATAN_MASS_H_

#include <latan/latan_globals.h>
#include <latan/latan_fit.h>

__BEGIN_DECLS

typedef struct
{
    size_t start;
    size_t end;
    double mean;
    double sig;
} plat;

latan_errno effmass(mat *res, mat *mprop, const int parity);
latan_errno effmass_PCAC(mat *res, mat *mprop_AP, mat *mprop_PP);
plat *search_plat(size_t *nplat, mat *data, mat *sigdata,\
                  const size_t ntmax, const double nsig, const double tol);
latan_errno fit_data_mass_fit_tune(fit_data *d, mat *fit_init, mat *prop,\
                                   mat *em, mat *sigem,            \
                                   const int parity);

__END_DECLS

#endif