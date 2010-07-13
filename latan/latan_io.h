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
	strbuf _buf;\
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
void io_get_prop_mark(strbuf prop_mark);
void io_set_prop_mark(const strbuf prop_mark);
void io_get_prop_idfmt(strbuf prop_idfmt);
void io_set_prop_idfmt(const strbuf prop_idfmt);

/* general I/O */
int get_nfile(const strbuf manifestfname);
latan_errno get_firstfname(strbuf fname, const strbuf manifestfname);

/* mat *I/O */
void mat_dump(FILE* stream, const mat *m);
#define mat_print(m) mat_dump(stdout,m)
int mat_load_nrow(const strbuf mark, const strbuf matid,\
				  const strbuf inputfname);
latan_errno mat_load(mat *m, const strbuf mark, const strbuf matid,\
					 const strbuf inputfname);
latan_errno mat_load_ar(mat **m, const strbuf mark,\
						const strbuf matid, const strbuf manifestfname);
latan_errno mat_save_plotdat(const mat *x, const mat *m, const strbuf fname);
latan_errno mat_save_plotdat_yerr(const mat *x, const mat *dat, const mat *yerr,\
								  const strbuf fname);
latan_errno mat_save_plotdat_xyerr(const mat *x, const mat *dat, const mat *xerr,\
								   const mat *yerr, const strbuf fname);

/* propagator I/O */
int hadron_getnt(const hadron *h, const ss_no source, const ss_no sink,\
				 const strbuf manfname);
latan_errno hadron_propbin(mat **prop, const hadron *h, const ss_no source,	\
						   const ss_no sink, const strbuf manfname,	\
						   const size_t binsize);

/* random generator state I/O */
latan_errno randgen_save_state(const strbuf prefname,\
							   const randgen_state state);
latan_errno randgen_load_state(randgen_state state, const strbuf prefname);

/* reampled sample I/O */
latan_errno rs_sample_save(const rs_sample *s, const strbuf f_name);
int rs_sample_load_nrow(const strbuf f_name);
int rs_sample_load_nsample(const strbuf f_name);
int rs_sample_load_method(const strbuf f_name);
latan_errno rs_sample_load(rs_sample *s, const strbuf f_name);

/* plot output */

__END_DECLS

#endif