#ifndef LATAN_IO_H_
#define LATAN_IO_H_

#include <latan/latan_globals.h>
#include <latan/latan_hadron.h>
#include <latan/latan_rand.h>
#include <latan/latan_statistics.h>

__BEGIN_DECLS

/* general I/O */
int get_nfile(const stringbuf manifestfname);
latan_errno get_firstfname(stringbuf fname, const stringbuf manifestfname);

/* mat I/O */
void mat_dump(FILE* stream, mat m);
#define mat_print(m) mat_dump(stdout,m)
int mat_load_nrow(const stringbuf mark, const stringbuf matid,\
				  const stringbuf inputfname);
latan_errno mat_load(mat m, const stringbuf mark, const stringbuf matid,\
					 const stringbuf inputfname);
latan_errno mat_load_ar(mat* m, const stringbuf mark,\
						const stringbuf matid, const stringbuf manifestfname);
latan_errno mat_save_plotdat(const mat m, const double xstart,\
							 const double xstep, const stringbuf fname);
latan_errno mat_save_plotdaterr(const mat dat, const mat sig,\
								const double xstart, const double xstep,\
								const stringbuf fname);

/* propagator I/O */
int hadron_getnt(const hadron h, const int source, const int sink,\
				 const stringbuf manfname);
latan_errno hadron_prop(mat* prop, const hadron h, const int source,\
						const int sink, const stringbuf manfname);

/* random generator state I/O */
latan_errno randgen_save_state(const stringbuf prefname,\
							   const randgen_state state);
latan_errno randgen_load_state(randgen_state state, const stringbuf prefname);

/* reampled sample I/O */
latan_errno rs_sample_save(const rs_sample s, const stringbuf f_name);
latan_errno rs_sample_load(rs_sample s, const stringbuf f_name);

__END_DECLS

#endif