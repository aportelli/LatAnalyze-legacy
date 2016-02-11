/* latan_includes.h, part of LatAnalyze library
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

#ifndef LATAN_INCLUDES_H_
#define LATAN_INCLUDES_H_

#include <config.h>

/* replacement functions */
#ifndef HAVE_ACOSH
#ifdef acosh
#undef acosh
#endif
#define acosh gsl_acosh
#else
#ifndef __cplusplus
extern double acosh(double x); /* acosh is not ANSI compliant */
#endif
#endif
#ifndef HAVE_STRTOK_R
#ifdef _OPENMP
#error "your system does not have strtok_r function which is needed for thread safety, try to compile without OpenMP support"
#else
#ifdef strtok_r
#undef strtok_r
#endif
#define strtok_r(str,sep,dumb) strtok(str,sep)
#endif
#endif

/* system includes */
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <float.h>
#include <unistd.h>
#include <sys/stat.h>
#include <gsl/gsl_math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* alias for status update */
#define USTAT(inst) LATAN_UPDATE_STATUS(status,inst)

/* memory allocation */
#define MALLOC(pt,typ,size)\
{\
    pt = (typ)(malloc((size_t)(size)*sizeof(*pt)));\
    if (pt == NULL)\
    {\
        LATAN_ERROR("memory allocation failed",LATAN_ENOMEM);\
    }\
}
#define MALLOC_ERRVAL(pt,typ,size,value)\
{\
    pt = (typ)(malloc((size_t)(size)*sizeof(*pt)));\
    if (pt == NULL)\
    {\
        LATAN_ERROR_VAL("memory allocation failed",LATAN_ENOMEM,value);\
    }\
}
#define MALLOC_NOERRET(pt,typ,size)\
{\
    pt = (typ)(malloc((size_t)(size)*sizeof(*pt)));\
    if (pt == NULL)\
    {\
        LATAN_ERROR_NORET("memory allocation failed",LATAN_ENOMEM);\
    }\
}
#define REALLOC(pt,pt_old,typ,size)\
{\
    pt = (typ)(realloc(pt_old,(size_t)(size)*sizeof(*pt)));\
    if (pt == NULL)\
    {\
        LATAN_ERROR("memory reallocation failed",LATAN_ENOMEM);\
    }\
}
#define REALLOC_ERRVAL(pt,pt_old,typ,size,value)\
{\
    pt = (typ)(realloc(pt_old,(size_t)(size)*sizeof(*pt)));\
    if (pt == NULL)\
    {\
        LATAN_ERROR_VAL("memory reallocation failed",LATAN_ENOMEM,value);\
    }\
}
#define REALLOC_NOERRET(pt,pt_old,typ,size)\
{\
    pt = (typ)(realloc(pt_old,(size_t)(size)*sizeof(*pt)));\
    if (pt == NULL)\
    {\
        LATAN_ERROR_NORET("memory reallocation failed",LATAN_ENOMEM);\
    }\
}
#define FREE(pt)\
{\
    if (pt != NULL)\
    {\
        free(pt);\
        pt = NULL;\
    }\
}
/* file opening */
#define FOPEN(f,fname,mode)\
{\
    f = fopen(fname,mode);\
    if (f == NULL)\
    {\
        strbuf _errmsg;\
        sprintf(_errmsg,"error opening file %s",fname);\
        LATAN_ERROR(_errmsg,LATAN_EFAULT);\
    }\
}
#define FOPEN_ERRVAL(f,fname,mode,value)\
{\
    f = fopen(fname,mode);\
    if (f == NULL)\
    {\
        strbuf _errmsg;\
        sprintf(_errmsg,"error opening file %s",fname);\
        LATAN_ERROR_VAL(_errmsg,LATAN_EFAULT,value);\
    }\
}
#define FOPEN_NOERRET(f,fname,mode)\
{\
    f = fopen(fname,mode);\
    if (f == NULL)\
    {\
        strbuf _errmsg;\
        sprintf(_errmsg,"error opening file %s",fname);\
        LATAN_ERROR_NORET(_errmsg,LATAN_EFAULT);\
    }\
}

#endif
