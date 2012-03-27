/* latan_io.h, part of LatAnalyze library
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

#ifndef LATAN_IO_H_
#define LATAN_IO_H_

#include <latan/latan_globals.h>
#include <latan/latan_hadron.h>
#include <latan/latan_mat.h>
#include <latan/latan_rand.h>
#include <latan/latan_statistics.h>

/* loop on the lines of a file */
#define BEGIN_FOR_LINE_F(str,f,lc)\
{\
    size_t _l;\
    (lc) = 0;\
    rewind(f);\
    while (!feof(f))\
    {\
        if (fgets(str,STRING_LENGTH,f))\
        {\
            _l = strlen(str);\
            if (str[_l-1] == '\n')\
            {\
                (lc)++;\
                (str)[_l-1] = '\0';\
            }\
            if ((_l > 0)&&(strspn(str," \n\t") < strlen(str)))\
            {

#define END_FOR_LINE_F\
            }\
        }\
    }\
}

#define BEGIN_FOR_LINE(str,f_name,lc)\
{\
    FILE* _f;\
    _f = fopen(f_name,"r");\
    BEGIN_FOR_LINE_F(str,_f,lc)\
    {

#define END_FOR_LINE\
    }\
    END_FOR_LINE_F;\
    fclose(_f);\
}

/* loop on the lines of a file with tokenization of each line (thread safe) */
/* field is assumed to be a non allocated strbuf*                           */
#define BEGIN_FOR_LINE_TOK_F(field,f,tok,nf,lc)\
{\
    strbuf _line;\
    BEGIN_FOR_LINE_F(_line,f,lc)\
    {\
        char *_line_ptr,*_save_ptr;\
        _line_ptr = strtok_r(_line,tok,&_save_ptr);\
        nf        = 0;\
        while(_line_ptr != NULL)\
        {\
            REALLOC_NOERRET(field,field,strbuf *,(size_t)(nf+1));\
            strbufcpy(field[nf],_line_ptr);\
            _line_ptr = strtok_r(NULL,tok,&_save_ptr);\
            nf++;\
        }\
        _save_ptr = NULL;

#define END_FOR_LINE_TOK_F(field)\
    }\
    END_FOR_LINE_F;\
    FREE(field);\
}

#define BEGIN_FOR_LINE_TOK(field,f_name,tok,nf,lc)\
{\
    FILE* _f;\
    _f = fopen(f_name,"r");\
    BEGIN_FOR_LINE_TOK_F(field,_f,tok,nf,lc)\
    {

#define END_FOR_LINE_TOK(field)\
    }\
    END_FOR_LINE_TOK_F(field);\
    fclose(_f);\
}

__BEGIN_DECLS

/* I/O format */
typedef enum
{
    IO_XML   = 0,\
    IO_ASCII = 1 \
} io_fmt_no;

latan_errno io_set_fmt(const io_fmt_no fmt);
io_fmt_no io_get_fmt(void);

/* I/O init/finish */
extern void (*io_init)(void);
extern void (*io_finish)(void);

/* general I/O */
int get_nfile(const strbuf manifestfname);
latan_errno get_firstfname(strbuf fname, const strbuf manifestfname);

/* mat I/O */
void mat_dump(FILE* stream, const mat *m, const strbuf fmt);
#define mat_print(m,fmt) mat_dump(stdout,m,fmt)
latan_errno mat_save_plotdat(mat *x, mat *m, const strbuf fname);
latan_errno mat_save_plotdat_yerr(mat *x, mat *dat, mat *yerr,\
                                  const strbuf fname);
latan_errno mat_save_plotdat_xyerr(mat *x, mat *dat, mat *xerr,\
                                   mat *yerr, const strbuf fname);
extern latan_errno (*mat_save)(const strbuf fname, const char mode,\
                               const mat *m, const strbuf name);
extern latan_errno (*mat_load_dim)(size_t dim[2], const strbuf fname,\
                                   const strbuf name);
extern latan_errno (*mat_load)(mat *m, const strbuf fname, const strbuf name);
latan_errno mat_load_subm(mat *m, const strbuf fname, const strbuf name,    \
                          const size_t k1, const size_t l1, const size_t k2,\
                          const size_t l2);

/* propagator I/O */
extern latan_errno (*prop_load_nt)(size_t *nt, const channel_no channel,\
                                   const quark_no q1, const quark_no q2,\
                                   const ss_no source, const ss_no sink,\
                                   strbuf fname);
extern latan_errno (*prop_load)(mat *prop, const channel_no channel, \
                                const quark_no q1, const quark_no q2,\
                                const ss_no source, const ss_no sink,\
                                strbuf fname);
extern latan_errno (*prop_save)(strbuf fname, const char mode, mat *prop, \
                                const strbuf channel,                     \
                                const quark_no q1, const quark_no q2,     \
                                const ss_no source, const ss_no sink,     \
                                const strbuf name);
latan_errno hadron_prop_load_bin(mat **prop, const hadron *h,              \
                                 const ss_no source, const ss_no sink,     \
                                 const strbuf manfname,const size_t binsize);
latan_errno hadron_prop_load_nt(size_t *nt, const hadron *h,               \
                                const ss_no source, const ss_no sink,      \
                                const strbuf manfname);

/* random generator state I/O */
extern latan_errno (*randgen_save_state)(const strbuf f_name, const char mode,\
                                         const rg_state state,                \
                                         const strbuf name);
extern latan_errno (*randgen_load_state)(rg_state state, const strbuf f_name,\
                                         const strbuf name);

/* resampled sample I/O */
extern latan_errno (*rs_sample_save)(const strbuf fname, const char mode,\
                                     const rs_sample *s);
latan_errno rs_sample_save_subsamp(const strbuf fname, const char mode,   \
                                   const rs_sample *s,                    \
                                   const size_t k1, const size_t k2);
extern latan_errno (*rs_sample_load_nrow)(size_t *nr, const strbuf fname,\
                                          const strbuf name);
extern latan_errno (*rs_sample_load_nsample)(size_t *nsample,   \
                                             const strbuf fname,\
                                             const strbuf name);
extern latan_errno (*rs_sample_load)(rs_sample *s, const strbuf fname,\
                                     const strbuf name);
latan_errno rs_sample_load_subsamp(rs_sample *s, const strbuf fname,  \
                                   const strbuf name, const size_t k1,\
                                   const size_t k2);

__END_DECLS

#endif
