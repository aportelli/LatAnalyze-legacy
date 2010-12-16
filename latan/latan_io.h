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
        if ((fgets(_buf,STRING_LENGTH,_f))&&(sscanf(_buf,"%s\n",str)>0))

#define END_FOR_LINE\
    }\
    fclose(_f);\
}

__BEGIN_DECLS

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
latan_errno prop_load(mat *prop, const channel_no channel, \
                      const quark_no q1, const quark_no q2,\
                      const ss_no source, const ss_no sink,\
                      strbuf fname);
latan_errno prop_load_nt(size_t *nt, const channel_no channel,\
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