#ifndef LATAN_IO_H_
#define LATAN_IO_H_

#include <latan/latan_globals.h>
#include <latan/latan_hadron.h>
#include <latan/latan_rand.h>
#include <latan/latan_statistics.h>

/* loop on lines of a file */
#define BEGIN_FOR_LINE(str,f_name)\
{\
	FILE* _f;\
	stringbuf _buf;\
	FOPEN(_f,f_name,"r");\
	while (!feof(_f))\
	{\
		if ((fgets(_buf,STRING_LENGTH,_f))&&(sscanf(_buf,"%s\n",str)>0))\
		{
#define END_FOR_LINE\
		}\
	}\
	fclose(_f);\
}

__BEGIN_DECLS

/* I/O options */
void io_get_prop_mark(stringbuf prop_mark);
void io_set_prop_mark(const stringbuf prop_mark);
void io_get_prop_idfmt(stringbuf prop_idfmt);
void io_set_prop_idfmt(const stringbuf prop_idfmt);

/* general I/O */
int get_nfile(const stringbuf manifestfname);
latan_errno get_firstfname(stringbuf fname, const stringbuf manifestfname);

/* mat *I/O */
void mat_dump(FILE* stream, const mat *m);
#define mat_print(m) mat_dump(stdout,m)
int mat_load_nrow(const stringbuf mark, const stringbuf matid,\
				  const stringbuf inputfname);
latan_errno mat_load(mat *m, const stringbuf mark, const stringbuf matid,\
					 const stringbuf inputfname);
latan_errno mat_load_ar(mat **m, const stringbuf mark,\
						const stringbuf matid, const stringbuf manifestfname);
latan_errno mat_save_plotdat(const mat *x, const mat *m, const stringbuf fname);
latan_errno mat_save_plotdat_yerr(const mat *x, const mat *dat, const mat *yerr,\
								  const stringbuf fname);
latan_errno mat_save_plotdat_xyerr(const mat *x, const mat *dat, const mat *xerr,\
								   const mat *yerr, const stringbuf fname);

/* propagator I/O */
int hadron_getnt(const hadron *h, const ss_no source, const ss_no sink,\
				 const stringbuf manfname);
latan_errno hadron_propbin(mat **prop, const hadron *h, const ss_no source,	\
						   const ss_no sink, const stringbuf manfname,	\
						   const size_t binsize);

/* random generator state I/O */
latan_errno randgen_save_state(const stringbuf prefname,\
							   const randgen_state state);
latan_errno randgen_load_state(randgen_state state, const stringbuf prefname);

/* reampled sample I/O */
latan_errno rs_sample_save(const rs_sample *s, const stringbuf f_name);
int rs_sample_load_nrow(const stringbuf f_name);
int rs_sample_load_nsample(const stringbuf f_name);
int rs_sample_load_method(const stringbuf f_name);
latan_errno rs_sample_load(rs_sample *s, const stringbuf f_name);

/* plot output */

__END_DECLS

#endif