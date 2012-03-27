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

static void no_error_handler(const char *reason, const char *file, int line,\
                             int no);

latan_error_handler_t *latan_error_handler = NULL;

void latan_error(const char *reason, const char *file, int line,\
                 int no)
{
    strbuf name,version;
    
    if (latan_error_handler) 
    {
        (*latan_error_handler)(reason,file,line,no);
        return;
    }
    
    latan_get_name(name);
    latan_get_version(version);
    fprintf(stderr,"%s v%s error %d: %s (%s:%d)\n",name,version,\
            no,reason,file,line);
    fflush(stderr);
    abort();
}

void latan_warning(const char *reason, const char *file, int line,\
                   int no)
{
    strbuf name,version;

    latan_get_name(name);
    latan_get_version(version);
    fprintf(stderr,"%s v%s warning %d: %s (%s:%d)\n",name,version,\
            no,reason,file,line);
    fflush(stderr);
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
