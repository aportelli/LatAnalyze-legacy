/* latan_hadron.h, part of LatAnalyze library
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

#ifndef LATAN_HADRON_H_
#define LATAN_HADRON_H_

#include <latan/latan_globals.h>

#define NOMIX 0
#define SUM 1
#define MEAN 2

#ifndef MAXPROP
#define MAXPROP 2
#endif
#ifndef MAXQUARKPST
#define MAXQUARKPST 2
#endif
#ifndef MAXQUARKST
#define MAXQUARKST 2
#endif

__BEGIN_DECLS

/* channels */
#define NCHANNEL 9
typedef enum
{
    ch_SS       = 0,    /* scalar           -> scalar           */
    ch_VV       = 1,    /* vector           -> vector           */
    ch_PP       = 2,    /* pseudoscalar     -> pseudoscalar     */
    ch_PA       = 3,    /* pseudoscalar     -> axial            */
    ch_AP       = 4,    /* axial            -> pseudoscalar     */
    ch_AA       = 5,    /* axial            -> axial            */
    ch_N        = 6,    /* nucleon-like     -> nucleon-like     */
    ch_Lambda   = 7,    /* Lambda-like      -> Lambda-like      */
    ch_Delta    = 8     /* Delta-like       -> Delta-like       */
} channel_no;

channel_no channel_no_get(const strbuf label);
void channel_id_set(const channel_no i, const strbuf new_id);
void channel_id_get(strbuf str, const channel_no i);

/* quarks */
#define NQUARK 4
typedef enum
{
    qu_l    = 0,    /* light    */
    qu_u    = 1,    /* up       */
    qu_d    = 2,    /* down     */
    qu_s    = 3     /* strange  */
} quark_no;

void quark_id_set(const quark_no i, const strbuf new_id);
void quark_id_get(strbuf str, const quark_no i);

/* sources/sinks */
#define NSS 3
typedef enum
{
    ss_P    = 0,    /* point    */
    ss_W    = 1,    /* wall     */
    ss_G    = 2     /* gaussian */
} ss_no;

void ss_id_set(const ss_no i, const strbuf new_id);
void ss_id_get(strbuf str, const ss_no i);

/* hadron structure */
typedef struct
{
    strbuf name;
    int channel[MAXPROP];
    int quarkst[MAXQUARKST][MAXQUARKPST];
    int chmix;
    int stmix;
    int parity;
} hadron;

/** allocation **/
hadron *hadron_create(void);
void hadron_destroy(hadron *h);

/** access **/
void hadron_set_2q_nomix(hadron *h, const strbuf name, const int parity,    \
                         const channel_no channel, const quark_no q1,       \
                         const quark_no q2);
void hadron_set_2q_2stmean(hadron *h, const strbuf name, const int parity,\
                           const channel_no channel, const quark_no q11,    \
                           const quark_no q12, const quark_no q21,          \
                           const quark_no q22);
void hadron_get_name(strbuf str, const hadron *h);

/* spectrum type */
typedef struct
{
    hadron **particle;
    size_t nparticle;
} spectrum;

/** allocation **/
spectrum *spectrum_create(const size_t nparticle);
spectrum *spectrum_create_qcd(void);
spectrum *spectrum_create_qcdqed(void);
void spectrum_destroy(spectrum *s);

/** access **/
hadron *spectrum_get(const spectrum *s, const strbuf part_name);

__END_DECLS

#endif
