/* latan_io_ascii.h, part of LatAnalyze library
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

#ifndef LATAN_IO_ASCII_H_
#define	LATAN_IO_ASCII_H_

#include <latan/latan_globals.h>
#include <latan/latan_hadron.h>
#include <latan/latan_statistics.h>

__BEGIN_DECLS

/* I/O init/finish */
void io_init_ascii(void);
void io_finish_ascii(void);

/* matrix I/O */
latan_errno mat_save_ascii(const strbuf fname, const char mode, const mat *m,\
                           const strbuf name);
latan_errno mat_load_dim_ascii(size_t dim[2], const strbuf fname,\
                               const strbuf name);
latan_errno mat_load_ascii(mat *m, const strbuf fname, const strbuf name);

/* propagator I/O */
latan_errno prop_load_nt_ascii(size_t *nt, const channel_no channel,\
                               const quark_no q1, const quark_no q2,\
                               const ss_no source, const ss_no sink,\
                               strbuf fname);

latan_errno prop_load_ascii(mat *prop, const channel_no channel, \
                            const quark_no q1, const quark_no q2,\
                            const ss_no source, const ss_no sink,\
                            strbuf fname);
latan_errno prop_save_ascii(strbuf fname, const char mode, mat *prop,\
                            const strbuf channel,                    \
                            const quark_no q1, const quark_no q2,    \
                            const ss_no source, const ss_no sink,    \
                            const strbuf name);

/* random generator state I/O */
latan_errno randgen_save_state_ascii(const strbuf fname, const char mode,   \
                                     const rg_state state, const strbuf name);

latan_errno randgen_load_state_ascii(rg_state state, const strbuf fname,\
                                     const strbuf name);

/* resampled sample I/O */
latan_errno rs_sample_save_ascii(const strbuf fname, const char mode,\
                                 const rs_sample *s);

latan_errno rs_sample_load_nrow_ascii(size_t *nr, const strbuf fname,\
                                      const strbuf name);
latan_errno rs_sample_load_nsample_ascii(size_t *nsample, const strbuf fname,\
                                         const strbuf name);
latan_errno rs_sample_load_ascii(rs_sample *s, const strbuf fname,\
                                 const strbuf name);

__END_DECLS

#endif
