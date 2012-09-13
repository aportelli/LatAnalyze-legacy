/* latan_error.c, part of LatAnalyze library
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

#include <latan/latan_includes.h>
#include <latan/latan_error.h>

#define WARNING_DUR 5

static void no_error_handler(const char *reason, const char *file, int line,\
                             int no);

static latan_error_handler_t *latan_error_handler = NULL;
static unsigned int warn_count = 0;
static long unsigned int warn_time = 0;

void latan_error(const char *reason, const char *file, int line,\
                 int no)
{
    strbuf name,version,prefix;
    
    if (latan_error_handler) 
    {
        (*latan_error_handler)(reason,file,line,no);
        return;
    }
    
    latan_get_name(name);
    latan_get_version(version);
    latan_get_msg_prefix(prefix);
    fprintf(stderr,"%s%s v%s error %d: %s (%s:%d)\n",prefix,name,version,\
            no,reason,file,line);
    fflush(stderr);
    abort();
}

void latan_warning(const char *reason, const char *file, int line,\
                   int no)
{
    strbuf name,version,prefix;

    if (latan_get_warn())
    {
        if (warn_time == 0)
        {
            warn_count = 0;
            warn_time = (long unsigned int)time(NULL);
        }
        else if ((long unsigned int)time(NULL) - warn_time > WARNING_DUR)
        {
            warn_count = 0;
            warn_time = (long unsigned int)time(NULL);
        }
        warn_count++;
        latan_get_name(name);
        latan_get_version(version);
        latan_get_msg_prefix(prefix);
        if (warn_count <= LATAN_MAX_WARNING)
        {
            fprintf(stderr,"%s%s v%s warning %d: %s (%s:%d)\n",prefix,name,\
                    version,no,reason,file,line);
            fflush(stderr);
        }
        else if (warn_count == LATAN_MAX_WARNING+1)
        {
            fprintf(stderr,"%s%s v%s warning: to much warnings in the last %d seconds\n",\
                    prefix,name,version,WARNING_DUR);
        }
    }
}

latan_error_handler_t*\
latan_set_error_handler(latan_error_handler_t *new_handler)
{
    latan_error_handler_t *previous_handler = latan_error_handler;
    latan_error_handler = new_handler;
    return previous_handler;
}


latan_error_handler_t *latan_set_error_handler_off (void)
{
    latan_error_handler_t *previous_handler = latan_error_handler;
    latan_error_handler = no_error_handler;
    return previous_handler;
}

void no_error_handler(const char *reason __dumb, const char *file __dumb,
                      int line __dumb, int no __dumb)
{
}
