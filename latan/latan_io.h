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
#include <latan/latan_mat.h>
#include <latan/latan_rand.h>
#include <latan/latan_statistics.h>

/* LatAnalayze separator in file paths */
#ifndef LATAN_PATH_SEP
#define LATAN_PATH_SEP ':'
#endif

/* include for text loop macros */
#include <stdio.h>
#include <string.h>

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
    FOPEN_NOERRET(_f,f_name,"r");\
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
    strbuf _line,_line_buf;\
    BEGIN_FOR_LINE_F(_line_buf,f,lc)\
    {\
        char *_line_ptr,*_save_ptr;\
        strbufcpy(_line,_line_buf);\
        _line_ptr = strtok_r(_line_buf,tok,&_save_ptr);\
        nf        = 0;\
        while(_line_ptr != NULL)\
        {\
            REALLOC_NOERRET(field,field,strbuf *,(size_t)(nf+1));\
            strbufcpy(field[nf],_line_ptr);\
            _line_ptr = strtok_r(NULL,tok,&_save_ptr);\
            nf++;\
        }\
        _save_ptr = NULL;\
        if (nf > 0)\
        {\

#define END_FOR_LINE_TOK_F(field)\
        }\
    }\
    END_FOR_LINE_F;\
    FREE(field);\
}

#define BEGIN_FOR_LINE_TOK(field,f_name,tok,nf,lc)\
{\
    FILE* _f;\
    FOPEN_NOERRET(_f,f_name,"r");\
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
void io_init(void);
void io_finish(void);

/* general I/O */
int  get_nfile(const strbuf manifestfname);
void get_firstfname(strbuf fname, const strbuf manifestfname);
void get_elname(strbuf fname, strbuf elname, const strbuf latan_path);

/* mat I/O */
void mat_dump(FILE* stream, const mat *m, const strbuf fmt);
#define mat_print(m,fmt) mat_dump(stdout,m,fmt)
latan_errno mat_save(const strbuf latan_path, const char mode, const mat *m);
latan_errno mat_save_subm(const strbuf latan_path, const char mode,      \
                          const mat *m, const size_t k1, const size_t l1,\
                          const size_t k2, const size_t l2);
latan_errno mat_load(mat *m, size_t *dim, const strbuf latan_path);
latan_errno mat_load_subm(mat *m, const strbuf latan_path, const size_t k1, \
                          const size_t l1, const size_t k2, const size_t l2);
latan_errno mat_ar_loadbin(mat **m, size_t *dim, const strbuf man_fname,\
                           const strbuf m_name, const size_t binsize);

/* random generator state I/O */
latan_errno randgen_save_state(const strbuf latan_path, const char mode,\
                               const rg_state state);
latan_errno randgen_load_state(rg_state state, const strbuf latan_path);

/* resampled sample I/O */
latan_errno rs_sample_save(const strbuf latan_path, const char mode,\
                           const rs_sample *s);
latan_errno rs_sample_save_subsamp(const strbuf latan_path, const char mode,\
                                   const rs_sample *s, const size_t k1,     \
                                   const size_t l1, const size_t k2,        \
                                   const size_t l2);
latan_errno rs_sample_load(rs_sample *s, size_t *nsample, size_t *dim,\
                           const strbuf latan_path);
latan_errno rs_sample_load_subsamp(rs_sample *s, const strbuf latan_path,\
                                   const size_t k1, const size_t l1,     \
                                   const size_t k2, const size_t l2);

__END_DECLS

#endif
