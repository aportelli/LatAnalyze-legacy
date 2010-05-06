#ifndef LATAN_HADRON_H_
#define LATAN_HADRON_H_

#include <latan/latan_globals.h>

#define NOMIX 0
#define SUM 1
#define MEAN 2

#define EVEN 0
#define ODD 1

#ifndef MAXPROP
#define MAXPROP 2
#endif
#ifndef MAXQUARKST
#define MAXQUARKST 2
#endif
#ifndef MAXISO
#define MAXISO 4
#endif

__BEGIN_DECLS

/* channels */
#define NCHANNEL 9
typedef enum
{
	ch_SS		= 0,
	ch_VV		= 1,
	ch_PP		= 2,
	ch_PA		= 3,
	ch_AP		= 4,
	ch_AA		= 5,
	ch_N		= 6,
	ch_Lambda	= 7,
	ch_Delta	= 8
} channel_no;

void channel_id_set(const channel_no i, const stringbuf new_id);
void channel_id_get(stringbuf str, const channel_no i);

/* quarks */
#define NQUARK 4
typedef enum
{
	qu_l	= 0,
	qu_u	= 1,
	qu_d	= 2,
	qu_s	= 3
} quark_no;

quark_no quark_no_get(const char c);
void quark_id_set(const quark_no i, const stringbuf new_id);
void quark_id_get(stringbuf str, const quark_no i);
void diquark_id_get(stringbuf str, const quark_no i1, const quark_no i2);
void triquark_id_get(stringbuf str, const quark_no i1, const quark_no i2,\
					 const quark_no i3);

/* hadron structure */
typedef struct
{
	stringbuf name;
	stringbuf channel[MAXPROP];
	stringbuf quarkst[MAXQUARKST];
	int chmix;
	int stmix;
	int parity;
}* hadron;

/** allocation **/
hadron hadron_create(void);
void hadron_destroy(hadron h);

/** access **/
void hadron_set_2q_nomix(hadron h, const stringbuf name, const int parity,	\
						 const channel_no channel, const quark_no q1,		\
						 const quark_no q2);
void hadron_set_2q_2stmean(hadron h, const stringbuf name, const int parity,\
						   const channel_no channel, const quark_no q11,	\
						   const quark_no q12, const quark_no q21,			\
						   const quark_no q22);

/* isohadron structure */
typedef struct
{
	stringbuf name;
	size_t isodim;
	const hadron* had[MAXISO];
} isohad;

/* spectrum type */
typedef struct
{
	hadron* particle;
	size_t nparticle;
}* spectrum;

/** allocation **/
spectrum spectrum_create(const size_t nparticle);
spectrum spectrum_create_qcd(void);
spectrum spectrum_create_qcdqed(void);
void spectrum_destroy(spectrum s);

/** access **/
hadron spectrum_get(const spectrum s, const stringbuf part_name);

__END_DECLS

#endif