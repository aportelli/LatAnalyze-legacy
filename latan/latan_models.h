/* latan_models.h, part of LatAnalyze library
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

#ifndef LATAN_MODELS_H_
#define LATAN_MODELS_H_

#include <latan/latan_globals.h>
#include <latan/latan_fit.h>

/** constant **/
extern fit_model fm_const;

/** exponential decay **/
extern fit_model fm_expdec;
extern fit_model fm_expdec_ex;
extern fit_model fm_expdec_splitsum;
extern fit_model fm_expdec_ex_splitsum;

/** hyperbolic cosine **/
extern fit_model fm_cosh;
extern fit_model fm_cosh_ex;
extern fit_model fm_cosh_splitsum;
extern fit_model fm_cosh_ex_splitsum;

#endif
