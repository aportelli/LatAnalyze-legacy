#ifndef LATAN_HADRON_H_
#define LATAN_HADRON_H_

#include <latan/latan_globals.h>

#define NOMIX 0
#define SUM 1
#define MEAN 2

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
	ch_SS		= 0,	/* scalar			-> scalar			*/
	ch_VV		= 1,	/* vector			-> vector			*/
	ch_PP		= 2,	/* pseudoscalar		-> pseudoscalar		*/
	ch_PA		= 3,	/* pseudoscalar		-> axial			*/
	ch_AP		= 4,	/* axial			-> pseudoscalar		*/
	ch_AA		= 5,	/* axial			-> axial			*/
	ch_N		= 6,	/* nucleon-like		-> nucleon-like		*/
	ch_Lambda	= 7,	/* Lambda-like		-> Lambda-like		*/
	ch_Delta	= 8		/* Delta-like		-> Delta-like		*/
} channel_no;

channel_no channel_no_get(const stringbuf label);
void channel_id_set(const channel_no i, const stringbuf new_id);
void channel_id_get(stringbuf str, const channel_no i);

/* quarks */
#define NQUARK 4
typedef enum
{
	qu_l	= 0,	/* light	*/
	qu_u	= 1,	/* up		*/
	qu_d	= 2,	/* down		*/
	qu_s	= 3		/* strange	*/
} quark_no;

quark_no quark_no_get(const char c);
void quark_id_set(const quark_no i, const stringbuf new_id);
void quark_id_get(stringbuf str, const quark_no i);
void diquark_id_get(stringbuf str, const quark_no i1, const quark_no i2);
void triquark_id_get(stringbuf str, const quark_no i1, const quark_no i2,\
					 const quark_no i3);

/* sources/sinks */
#define NSS 3
typedef enum
{
	ss_P	= 0,	/* point	*/
	ss_W	= 1,	/* wall		*/
	ss_G	= 2		/* gaussian	*/
} ss_no;

ss_no ss_no_get(const char c);
void ss_id_set(const ss_no i, const stringbuf new_id);
void ss_id_get(stringbuf str, const ss_no i);

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
void hadron_get_name(stringbuf str, const hadron h);

/* spectrum type */
typedef struct
{
	hadron *particle;
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