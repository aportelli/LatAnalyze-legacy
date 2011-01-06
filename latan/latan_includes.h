/* latan_includes.h, part of LatAnalyze library
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

#ifndef LATAN_INCLUDES_H_
#define LATAN_INCLUDES_H_

#include "../config.h"

/* replacement functions */
#ifndef HAVE_ACOSH
#ifdef acosh
#undef acosh
#endif
#define acosh gsl_acosh
#endif

/* system includes */
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <unistd.h>
#include <stdarg.h>
#include <limits.h>
#include <float.h>
#include <gsl/gsl_math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#endif