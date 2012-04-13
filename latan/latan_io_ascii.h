/* latan_io_ascii.h, part of LatAnalyze library
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

#ifndef LATAN_IO_ASCII_H_
#define	LATAN_IO_ASCII_H_

#include <latan/latan_globals.h>
#include <latan/latan_mat.h>
#include <latan/latan_statistics.h>

__BEGIN_DECLS

/* I/O init/finish */
void io_init_ascii(void);
void io_finish_ascii(void);

/* matrix I/O */
latan_errno mat_save_ascii(const strbuf fname, const char mode, const mat *m,\
                           const strbuf name);
latan_errno mat_load_ascii(mat *m, size_t *dim, const strbuf fname,\
                           const strbuf name);

/* random generator state I/O */
latan_errno randgen_save_state_ascii(const strbuf fname, const char mode,   \
                                     const rg_state state, const strbuf name);
latan_errno randgen_load_state_ascii(rg_state state, const strbuf fname,\
                                     const strbuf name);

/* resampled sample I/O */
latan_errno rs_sample_save_ascii(const strbuf fname, const char mode,\
                                 const rs_sample *s, const strbuf name);
latan_errno rs_sample_load_ascii(rs_sample *s, size_t *nsample, size_t *dim,\
                                 const strbuf fname, const strbuf name);

__END_DECLS

#endif
