/* latan_globals.c, part of LatAnalyze library
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

#include <latan/latan_globals.h>
#include <latan/latan_includes.h>

/* Prototypes for CUBLAS on/off, we don't want to include cublas.h here
 * because it is not ANSI compliant
 */
#ifdef HAVE_LIBCUBLAS
int cublasInit(void);
int cublasShutdown(void);
#endif

typedef struct
{
    strbuf name;
    strbuf version;
    int verb;
    int mat_op;
} latan_env;

static latan_env env = 
{
    PACKAGE_NAME,       \
    PACKAGE_VERSION,    \
    QUIET,              \
    CPU_MAT_OP,         \
};

static bool latan_is_cublas_run = false;

void latan_get_name(strbuf name)
{
    STRBUFCPY(name,env.name);
}

void latan_get_version(strbuf version)
{
    STRBUFCPY(version,env.version);
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

int latan_get_mat_op(void)
{
    return env.mat_op;
}

latan_errno latan_set_mat_op(const int mat_op)
{
    if ((mat_op < 0)||(mat_op > 1))
    {
        LATAN_ERROR("matrix operator flag invalid",LATAN_EINVAL);
    }
#ifndef HAVE_LIBCUBLAS
    if (mat_op == GPU_MAT_OP)
    {
        LATAN_ERROR("GPU operation support was not compiled",LATAN_EINVAL);
    }
    latan_is_cublas_run = false;
#else
    switch (mat_op)
    {
        case CPU_MAT_OP:
            if (latan_is_cublas_run)
            {
                cublasShutdown();
                latan_is_cublas_run = false;
            }
            break;
        case GPU_MAT_OP:
            if (!latan_is_cublas_run)
            {
                cublasInit();
                latan_is_cublas_run = true;
            }
            break;
    }
#endif
    env.mat_op = mat_op;
    
    return LATAN_SUCCESS;
}

void latan_printf(const int verb, const strbuf fmt, ...)
{
    va_list args;
    strbuf head,tail,name,debug,version;
    
    if ((latan_get_verb() >= verb)&&(verb >= VERB))
    {
        latan_get_name(name);
        latan_get_version(version);
        if (verb == DEBUG)
        {
            STRBUFCPY(debug," - DEBUG");
        }
        else
        {
            STRBUFCPY(debug,"");
        }
        
        sprintf(head,"[%s v%s%s]",name,version,debug);
        va_start(args,fmt);
        vsprintf(tail,fmt,args);
        va_end(args);
        printf("%s %s",head,tail);
    }
}