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

typedef struct
{
    strbuf name;
    strbuf version;
    int verb;
} latan_env;

static latan_env env = 
{
    PACKAGE_NAME,       \
    PACKAGE_VERSION,    \
    QUIET               \
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

int latan_get_verb(void)
{
    return env.verb;
}

latan_errno latan_set_verb(const int verb)
{
    if ((verb < 0)||(verb > 2))
    {
        LATAN_ERROR("verbosity level invalid",LATAN_EINVAL);
    }
    
    env.verb = verb;
    
    return LATAN_SUCCESS;
}

/*                        LatAnalyze message function                         */
/******************************************************************************/
void latan_printf(const int verb, const strbuf fmt, ...)
{
    va_list args;
    strbuf head,tail,name,debug,version;
    
    if ((latan_get_verb() >= verb)&&(verb >= VERB))
    {
        latan_get_name(name);
        latan_get_version(version);
        if (verb >= DEBUG1)
        {
            strbufcpy(debug," - DEBUG");
        }
        else
        {
            strbufcpy(debug,"");
        }
        
        sprintf(head,"[%s v%s%s]",name,version,debug);
        va_start(args,fmt);
        vsprintf(tail,fmt,args);
        va_end(args);
        printf("%s %s",head,tail);
    }
}

/*                           string operations                                */
/******************************************************************************/
char * strbufcpy(strbuf a, const strbuf b)
{
    char *pt;

    pt = strncpy(a,b,STRING_LENGTH);
    
    return pt;
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

