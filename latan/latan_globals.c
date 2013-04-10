/* latan_globals.c, part of LatAnalyze library
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

#include <latan/latan_globals.h>
#include <latan/latan_includes.h>
#include <gsl/gsl_nan.h>

typedef struct
{
    strbuf name;
    strbuf version;
    strbuf msg_prefix;
    int verb;
    bool warn;
    bool use_car_ret;
} latan_env;

static latan_env env = 
{
    PACKAGE_NAME,       \
    PACKAGE_VERSION,    \
    "",                 \
    QUIET,              \
    true,               \
    true
};

/*                         LatAnalyze environment access                      */
/******************************************************************************/
void latan_get_name(strbuf name)
{
    strbufcpy(name,env.name);
}

void latan_get_version(strbuf version)
{
    strbufcpy(version,env.version);
}

void latan_get_msg_prefix(strbuf msg_prefix)
{
    strbufcpy(msg_prefix,env.msg_prefix);
}

void latan_set_msg_prefix(const strbuf msg_prefix)
{
    strbufcpy(env.msg_prefix,msg_prefix);
}

int latan_get_verb(void)
{
    return env.verb;
}

latan_errno latan_set_verb(const int verb)
{
    if ((verb < QUIET)||(verb > DEBUG2))
    {
        LATAN_ERROR("verbosity level invalid",LATAN_EINVAL);
    }
    
    env.verb = verb;
    
    return LATAN_SUCCESS;
}

bool latan_get_warn(void)
{
    return env.warn;
}

void latan_set_warn(const bool warn)
{
    env.warn = warn;
}

bool latan_get_use_car_ret(void)
{
    return env.use_car_ret;
}

void latan_set_use_car_ret(const bool use_car_ret)
{
    env.use_car_ret = use_car_ret;
}

/*                        LatAnalyze message function                         */
/******************************************************************************/
void latan_printf(const int verb, const strbuf fmt, ...)
{
    va_list args;
    strbuf head,tail,name,debug,version,prefix;
    
    if ((latan_get_verb() >= verb)&&(verb >= VERB))
    {
        latan_get_name(name);
        latan_get_version(version);
        latan_get_msg_prefix(prefix);
        if (verb >= DEBUG1)
        {
            strbufcpy(debug," - DEBUG");
        }
        else
        {
            strbufcpy(debug,"");
        }
        
        sprintf(head,"%s[%s v%s%s]",prefix,name,version,debug);
        va_start(args,fmt);
        vsprintf(tail,fmt,args);
        va_end(args);
        printf("%s %s",head,tail);
    }
}

/*                           string operations                                */
/******************************************************************************/
char * strbufcat(strbuf a, const strbuf b)
{
    if (strlen(a)+strlen(b)+1 > STRING_LENGTH)
    {
        LATAN_WARNING("string buffer overflow",LATAN_FAILURE);
    }
    
    return strncat(a,b,STRING_LENGTH-strlen(a)-1);
}

int strbufcmp(const strbuf a, const strbuf b)
{
    return strncmp(a,b,STRING_LENGTH-1);
}

char * strbufcpy(strbuf a, const strbuf b)
{
    return strncpy(a,b,STRING_LENGTH);
}

/*                                NaN                                         */
/******************************************************************************/
double latan_nan(void)
{
    return GSL_NAN;
}

bool   latan_isnan(const double x)
{
    return gsl_isnan(x);
}

double latan_inf(void)
{
    return GSL_POSINF;
}

bool   latan_isinf(const double x)
{
    return gsl_isinf(x);
}
/*                          endianness management                             */
/******************************************************************************/
endian_no latan_get_endianness(void)
{
    const int i = 1;
    
    if ((*(const char*)&i) == 0)
    {
        return BE;
    }
    else
    {
        return LE;
    }
}

int latan_swap_byte_i(int x)
{
    unsigned char b1, b2, b3, b4;
    
    b1 = ( x >> 0  ) & 0x000000FF;
    b2 = ( x >> 8  ) & 0x000000FF;
    b3 = ( x >> 16 ) & 0x000000FF;
    b4 = ( x >> 24 ) & 0x000000FF;
    
    return ((int)b1 << 24) + ((int)b2 << 16) + ((int)b3 << 8) + b4;
}

double latan_swap_byte_d(double x)
{
    double y;
    unsigned char *dst = (unsigned char *)&y;
    unsigned char *src = (unsigned char *)&x;
    
    dst[0] = src[7];
    dst[1] = src[6];
    dst[2] = src[5];
    dst[3] = src[4];
    dst[4] = src[3];
    dst[5] = src[2];
    dst[6] = src[1];
    dst[7] = src[0];
    
    return y;
}

int  latan_conv_endianness_i(int x, endian_no source)
{
    if (latan_get_endianness() == source)
    {
        return x;
    }
    else
    {
        return latan_swap_byte_i(x);
    }
}

double latan_conv_endianness_d(double x, endian_no source)
{
    if (latan_get_endianness() == source)
    {
        return x;
    }
    else
    {
        return latan_swap_byte_d(x);
    }
}
