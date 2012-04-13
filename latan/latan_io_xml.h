/* latan_io_xml.h, part of LatAnalyze library
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

#ifndef LATAN_IO_XML_H_
#define	LATAN_IO_XML_H_

#include <latan/latan_globals.h>
#include <latan/latan_mat.h>
#include <latan/latan_statistics.h>

/* I/O init/finish */
void io_init_xml(void);
void io_finish_xml(void);

/* matrix I/O */
latan_errno mat_save_xml(const strbuf fname, const char mode, const mat *m,\
                         const strbuf name);
latan_errno mat_load_dim_xml(size_t dim[2], const strbuf fname,\
                             const strbuf name);
latan_errno mat_load_xml(mat *m, size_t *dim, const strbuf fname,\
                         const strbuf name);

/* random generator state I/O */
latan_errno randgen_save_state_xml(const strbuf f_name, const char mode,
                                   const rg_state state, const strbuf name);
latan_errno randgen_load_state_xml(rg_state state, const strbuf f_name,\
                                   const strbuf name);

/* resampled sample I/O */
latan_errno rs_sample_save_xml(const strbuf fname, const char mode,\
                               const rs_sample *s, const strbuf name);
latan_errno rs_sample_load_xml(rs_sample *s, size_t *nsample, size_t *dim,\
                               const strbuf fname, const strbuf name);

#endif
