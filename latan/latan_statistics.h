/* latan_statistics.h, part of LatAnalyze library
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
double mat_elsum(mat *m);
double mat_elmean(mat *m);
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

/** estimators **/
latan_errno rs_sample_cov(mat *cov, const rs_sample *s, const rs_sample *t);
#define rs_sample_var(cov,s) rs_sample_cov(cov,s,s);
latan_errno rs_sample_covp(mat *cov, const rs_sample *s, const rs_sample *t);
#define rs_sample_varp(cov,s) rs_sample_covp(cov,s,s);

/* resampling function */
latan_errno resample(rs_sample *s, mat **dat, const size_t ndat,               \
                     const size_t nobs, rs_func *f, unsigned int resamp_method,\
                     void *param);

/* useful rs_func */
latan_errno rs_mean(mat *res, mat **dat, const size_t ndat,\
                    const size_t sampno, void *nothing);
latan_errno rs_finite_diff(mat *res, mat **dat, const size_t ndat,\
                           const size_t sampno, void *nothing);
latan_errno rs_effmass(mat *res, mat **dat, const size_t ndat,\
                       const size_t sampno, void *parity);
latan_errno rs_effmass_PCAC(mat *res, mat **dat, const size_t ndat,\
                            const size_t sampno, void *nothing);

__END_DECLS

#endif