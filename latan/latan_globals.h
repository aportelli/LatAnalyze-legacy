/* latan_globals.h, part of LatAnalyze library
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

#ifndef LATAN_GLOBALS_H_
#define LATAN_GLOBALS_H_

/* macros for including from C++ code */
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS
#define __END_DECLS
#endif

__BEGIN_DECLS

/* attribute to switch off unused warnings with gcc */
#ifdef __GNUC__
#define __dumb __attribute__((unused))
#else
#define __dumb
#endif

/* boolean type */
#ifndef __cplusplus
#define false (1==0)
#define true  (1==1)
typedef int bool;
#endif

/* string buffers */
#ifndef STRING_LENGTH
#define STRING_LENGTH 1024
#endif

typedef char strbuf[STRING_LENGTH];

__END_DECLS

/* verbosity flags */
#define QUIET  0
#define VERB   1
#define DEBUG1 2
#define DEBUG2 3

/* size_t type */
#include <sys/types.h>

/* error handling */
#include <latan/latan_error.h>

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

__BEGIN_DECLS

/* LatAnalyze environment access */
void latan_get_name(strbuf name);
void latan_get_version(strbuf version);
void latan_get_msg_prefix(strbuf msg_prefix);
void latan_set_msg_prefix(const strbuf msg_prefix);
int  latan_get_verb(void);
#ifdef LATAN_ERRNO_DEF
latan_errno latan_set_verb(const int verb);
#endif
bool latan_get_warn(void);
void latan_set_warn(const bool warn);
bool latan_get_use_car_ret(void);
void latan_set_use_car_ret(const bool use_car_ret);


/* LatAnalyze message function */
void latan_printf(const int verb, const strbuf fmt, ...);

/* string operations */
char * strbufcat(strbuf a, const strbuf b);
int    strbufcmp(const strbuf a, const strbuf b);
char * strbufcpy(strbuf a, const strbuf b);


/* NaN */
double latan_nan(void);
bool   latan_isnan(const double x);

__END_DECLS

#endif
