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

/* endianness flags */
typedef enum
{
    LE = 0,\
    BE = 1
} endian_no;

/* size_t type */
#include <sys/types.h>

/* error handling */
#include <latan/latan_error.h>

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

/* nan/inf */
double latan_nan(void);
bool   latan_isnan(const double x);
double latan_inf(void);
bool   latan_isinf(const double x);

/* endianness management */
endian_no latan_get_endianness(void);
int       latan_swap_byte_i(int x);
double    latan_swap_byte_d(double x);
int       latan_conv_endianness_i(int x, endian_no source);
double    latan_conv_endianness_d(double x, endian_no source);

__END_DECLS

#endif
