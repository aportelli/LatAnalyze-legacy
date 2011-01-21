/* latan_io.h, part of LatAnalyze library
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

#ifndef LATAN_IO_H_
#define LATAN_IO_H_

#include <latan/latan_globals.h>
#include <latan/latan_hadron.h>
#include <latan/latan_io_xml.h>
#include <latan/latan_rand.h>
#include <latan/latan_statistics.h>

/* loop on lines of a file */
#define BEGIN_FOR_LINE(str,f_name)\
{\
    FILE* _f;\
    strbuf _buf;\
    FOPEN(_f,f_name,"r");\
    while (!feof(_f))\
    {\
        if ((fgets(_buf,STRING_LENGTH,_f))&&(sscanf(_buf,"%s\n",str)>0))

#define END_FOR_LINE\
    }\
    fclose(_f);\
}

__BEGIN_DECLS

/* I/O init/finish */
extern void (*io_init)(void);
extern void (*io_finish)(void);

/* general I/O */
int get_nfile(const strbuf manifestfname);
latan_errno get_firstfname(strbuf fname, const strbuf manifestfname);

/* mat *I/O */
void mat_dump(FILE* stream, mat *m, const strbuf fmt);
#define mat_print(m,fmt) mat_dump(stdout,m,fmt)
latan_errno mat_save_plotdat(mat *x, mat *m, const strbuf fname);
latan_errno mat_save_plotdat_yerr(mat *x, mat *dat, mat *yerr,\
                                  const strbuf fname);
latan_errno mat_save_plotdat_xyerr(mat *x, mat *dat, mat *xerr,\
                                   mat *yerr, const strbuf fname);

/* propagator I/O */
extern latan_errno (*prop_load_nt)(size_t *nt, const channel_no channel,\
                                   const quark_no q1, const quark_no q2,\
                                   const ss_no source, const ss_no sink,\
                                   strbuf fname);
extern latan_errno (*prop_load)(mat *prop, const channel_no channel, \
                                const quark_no q1, const quark_no q2,\
                                const ss_no source, const ss_no sink,\
                                strbuf fname);
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
extern latan_errno (*rs_sample_load_nrow)(size_t *nr, const strbuf fname,\
                                          const strbuf name);
extern latan_errno (*rs_sample_load_nsample)(size_t *nsample,   \
                                             const strbuf fname,\
                                             const strbuf name);
extern latan_errno (*rs_sample_load)(rs_sample *s, const strbuf fname,\
                                     const strbuf name);

__END_DECLS

#endif