#ifndef LATAN_IO_H_
#define LATAN_IO_H_

#include <latan/globals.h>
#include <latan/mat.h>
#include <latan/rand.h>

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

/* random generator state I/O */
latan_errno randgen_save_state(const stringbuf prefname,\
							   const randgen_state state);
latan_errno randgen_load_state(randgen_state state, const stringbuf prefname);

__END_DECLS

#endif