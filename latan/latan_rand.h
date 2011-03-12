/* latan_rand.h, part of LatAnalyze library
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

#ifndef LATAN_RAND_H_
#define LATAN_RAND_H_

#include <latan/latan_globals.h>

__BEGIN_DECLS

/* type for random generator state */
#define RLXG_STATE_SIZE 105
typedef int rg_state[RLXG_STATE_SIZE];

void randgen_init(const int seed);
void randgen_init_from_time(void);
void randgen_set_state(rg_state state);
void randgen_get_state(rg_state state);
double rand_u(double a, double b);
unsigned int rand_ud(const unsigned int n);
double rand_n(const double mean, const double sigma);

__END_DECLS

#endif
