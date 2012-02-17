/* latan_statistics.h, part of LatAnalyze library
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

#ifndef LATAN_STATISTICS_H_
#define LATAN_STATISTICS_H_

#include <latan/latan_globals.h>
#include <latan/latan_rand.h>

#define BOOT 0
#define JACK(depth) ((depth >= 1) ? depth : 1)

__BEGIN_DECLS

typedef latan_errno rs_func(mat *res, mat **dat, const size_t ndat,\
                            const size_t sampno, void *param);

/* elementary estimators */
double mat_elsum(const mat *m);
double mat_elmean(const mat *m);

latan_errno mat_mean(mat *mean, mat **m, const size_t size);
latan_errno mat_cov(mat *cov, mat **m, mat **n, const size_t size);
latan_errno mat_cov_m(mat *cov, mat **m, mat **n, const size_t size,\
                      mat *m_mean, mat *n_mean);
#define mat_var(var,m,size) mat_cov(var,m,m,size)
#define mat_var_m(var,m,size,mean) mat_cov(var,m,m,size,mean,mean)
latan_errno mat_covp(mat *cov, mat **m, mat **n, const size_t size);
latan_errno mat_covp_m(mat *cov, mat **m, mat **n, const size_t size,\
                       mat *m_mean, mat *n_mean);
#define mat_varp(var,m,size) mat_covp(var,m,m,size)
#define mat_varp_m(var,m,size,mean) mat_covp_m(var,m,m,size,mean,mean)

/* percentiles */
double ar_percentile(const double *data, const double *w, const size_t *sind,\
                     const size_t ndata, const double p);
double mat_elpercentile_with_sind(const mat *m, const mat *w,       \
                                  const size_t *sind, const double p);
double mat_elpercentile(const mat *m, const mat *w, const double p);

/* confidence interval */
double conf_int(double ci[2], const mat *m, const mat *w, const double nsig);

/* data binning */
latan_errno mat_ar_bin(mat **bindat, mat **dat, const size_t ndat,\
                       const size_t binsize);

/* histogram */
latan_errno histogram(mat *hist, mat *data, const double xmin,\
                      const double xmax, const size_t nint);

/* resampled sample type */
typedef struct
{
    strbuf name;
    mat *cent_val;
    mat **sample;
    size_t nsample;
    rg_state gen_state;
} rs_sample;

/** jackknife sample number calculation **/
size_t jackknife_nsample(const size_t ndat, const size_t jk_depth);

/** allocation **/
rs_sample *rs_sample_create(const size_t init_nrow, const size_t nsample);
void rs_sample_destroy(rs_sample *s);

/** access **/
size_t rs_sample_get_nrow(const rs_sample *s);
size_t rs_sample_get_nsample(const rs_sample *s);
void rs_sample_get_name(strbuf name, const rs_sample *s);
mat *rs_sample_pt_cent_val(const rs_sample *s);
mat *rs_sample_pt_sample(const rs_sample *s, const size_t i);
void rs_sample_set_name(rs_sample *s, const strbuf name);
latan_errno rs_sample_get_subsamp(rs_sample *s_a, const rs_sample *s_b,\
                                  const size_t k1, const size_t k2);
latan_errno rs_sample_set_subsamp(rs_sample *s_a, const rs_sample *s_b,\
                                  const size_t k1, const size_t k2);

/** estimators **/
latan_errno rs_sample_cov(mat *cov, const rs_sample *s, const rs_sample *t);
#define rs_sample_var(cov,s) rs_sample_cov(cov,s,s);
latan_errno rs_sample_covp(mat *cov, const rs_sample *s, const rs_sample *t);
#define rs_sample_varp(cov,s) rs_sample_covp(cov,s,s);

/* resampling function */
latan_errno resample(rs_sample *s, mat **dat, const size_t ndat, rs_func *f,\
                     unsigned int resamp_method, void *param);

/* useful rs_func */
latan_errno rs_mean(mat *res, mat **dat, const size_t ndat,\
                    const size_t sampno, void *nothing);

/* operations */
typedef latan_errno mat_unop(mat *a, const mat *b);
typedef latan_errno mat_unops(mat *a, const double s);
typedef latan_errno mat_binop(mat *a, const mat *b, const mat *c);
typedef latan_errno mat_binops(mat *a, const mat *b, const double s);

latan_errno rs_sample_unop(rs_sample *s_a, const rs_sample *s_b, mat_unop *f);
latan_errno rs_sample_unops(rs_sample *s_a, const double s, mat_unops *f);
latan_errno rs_sample_binop(rs_sample *s_a, const rs_sample *s_b, \
                            const rs_sample *s_c, mat_binop *f);
latan_errno rs_sample_binops(rs_sample *s_a, const rs_sample *s_b,\
                             const double s, mat_binops *f);
#define rs_sample_abs(s_a,s_b)      rs_sample_unop(s_a,s_b,&mat_abs)
#define rs_sample_eqabs(s_a)        rs_sample_abs(s_a,s_a)
#define rs_sample_add(s_a,s_b,s_c)  rs_sample_binop(s_a,s_b,s_c,&mat_add)
#define rs_sample_eqadd(s_a,s_b)    rs_sample_add(s_a,s_a,s_b)
#define rs_sample_adds(s_a,s_b,s)   rs_sample_binops(s_a,s_b,s,&mat_adds)
#define rs_sample_eqadds(s_a,s)     rs_sample_adds(s_a,s_a,s)
#define rs_sample_cst(s_a,s)        rs_sample_unops(s_a,s,&mat_cst)
#define rs_sample_divp(s_a,s_b,s_c) rs_sample_binop(s_a,s_b,s_c,&mat_divp)
#define rs_sample_eqdivp(s_a,s_b)   rs_sample_divp(s_a,s_a,s_b)
#define rs_sample_mulp(s_a,s_b,s_c) rs_sample_binop(s_a,s_b,s_c,&mat_mulp)
#define rs_sample_eqmulp(s_a,s_b)   rs_sample_mulp(s_a,s_a,s_b)
#define rs_sample_invp(s_a,s_b)     rs_sample_unop(s_a,s_b,&mat_invp)
#define rs_sample_eqinvp(s_a)       rs_sample_invp(s_a,s_a)
#define rs_sample_muls(s_a,s_b,s)   rs_sample_binops(s_a,s_b,s,&mat_muls)
#define rs_sample_eqmuls(s_a,s)     rs_sample_muls(s_a,s_a,s)
#define rs_sample_sqrt(s_a,s_b)     rs_sample_unop(s_a,s_b,&mat_sqrt)
#define rs_sample_eqsqrt(s_a)       rs_sample_sqrt(s_a,s_a)
#define rs_sample_expp(s_a,s_b)     rs_sample_unop(s_a,s_b,&mat_expp)
#define rs_sample_eqexpp(s_a)       rs_sample_expp(s_a,s_a)
#define rs_sample_sub(s_a,s_b,s_c)  rs_sample_binop(s_a,s_b,s_c,&mat_add)
#define rs_sample_eqsub(s_a,s_b)    rs_sample_sub(s_a,s_a,s_b)

__END_DECLS

#endif
