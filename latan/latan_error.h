/* latan_error.h, part of LatAnalyze library
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

#ifndef LATAN_ERROR_H_
#define LATAN_ERROR_H_

#include <latan/latan_globals.h>

__BEGIN_DECLS

typedef enum
{
    LATAN_FAILURE   = -1,   /* generic failure statement            */
    LATAN_SUCCESS   = 0,    /* all is going well !                  */
    LATAN_EDOM      = 1,    /* input domain error, e.g sqrt(-1)     */
    LATAN_EFAULT    = 3,    /* invalid pointer                      */
    LATAN_EINVAL    = 4,    /* invalid argument supplied by user    */
    LATAN_ENOMEM    = 8,    /* malloc error                         */
    LATAN_EBADLEN   = 19,   /* dimension error                      */
    LATAN_ENOTSQR   = 20,   /* matrix is not square error           */
    LATAN_ELATSYN   = 33,   /* syntax error                         */
    LATAN_ESYSTEM   = 34    /* system error                         */
} latan_errno;

#ifndef LATAN_ERRNO_DEF
#define LATAN_ERRNO_DEF
#endif

typedef void latan_error_handler_t(const char *, const char *, int, int);

void latan_error(const char *reason, const char *file, int line, int no);
void latan_warning(const char *reason, const char *file, int line, int no);
latan_error_handler_t*\
latan_set_error_handler(latan_error_handler_t *new_handler);
latan_error_handler_t *latan_set_error_handler_off(void);

#define LATAN_ERROR(reason,no)\
{\
    strbuf _freason;\
    sprintf(_freason,"%s (function %s)",reason,__FUNCTION__);\
    latan_error(_freason,__FILE__,__LINE__,no);\
    return no;\
}

#define LATAN_ERROR_VAL(reason,no,value)\
{\
    strbuf _freason;\
    sprintf(_freason,"%s (function %s)",reason,__FUNCTION__);\
    latan_error(_freason,__FILE__,__LINE__,no);\
    return value;\
}

#define LATAN_ERROR_VOID(reason,no)\
{\
    strbuf _freason;\
    sprintf(_freason,"%s (function %s)",reason,__FUNCTION__);\
    latan_error(_freason,__FILE__,__LINE__,no);\
    return;\
}

#define LATAN_ERROR_NULL(reason,no)\
LATAN_ERROR_VAL(reason,no,NULL)

#define LATAN_ERROR_NORET(reason,no)\
{\
    strbuf _freason;\
    sprintf(_freason,"%s (function %s)",reason,__FUNCTION__);\
    latan_error(_freason,__FILE__,__LINE__,no);\
}\

#define LATAN_WARNING(reason,no)\
{\
    strbuf _freason;\
    sprintf(_freason,"%s (function %s)",reason,__FUNCTION__);\
    latan_warning(_freason,__FILE__,__LINE__,no);\
}

#define LATAN_ERROR_SELECT_2(a,b)       \
((a) != LATAN_SUCCESS ? (a) : ((b) != LATAN_SUCCESS ? (b) : LATAN_SUCCESS))
#define LATAN_ERROR_SELECT_3(a,b,c)     \
((a) != LATAN_SUCCESS ? (a) : LATAN_ERROR_SELECT_2(b,c))
#define LATAN_ERROR_SELECT_4(a,b,c,d)   \
((a) != LATAN_SUCCESS ? (a) : LATAN_ERROR_SELECT_3(b,c,d))
#define LATAN_ERROR_SELECT_5(a,b,c,d,e) \
((a) != LATAN_SUCCESS ? (a) : LATAN_ERROR_SELECT_4(b,c,d,e))
#define LATAN_UPDATE_STATUS(status,instruction)\
{\
    int _cstat;\
    _cstat = instruction;\
    status = LATAN_ERROR_SELECT_2(status,_cstat);\
}

__END_DECLS

#endif
