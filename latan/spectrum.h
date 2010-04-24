#ifndef LATAN_SPECTRUM_H_
#define LATAN_SPECTRUM_H_

#include <latan/globals.h>
#include <latan/hadron.h>

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

void quark_id_set(const quark_no i, const stringbuf new_id);
void quark_id_get(stringbuf str, const quark_no i);
void diquark_id_get(stringbuf str, const quark_no i1, const quark_no i2);
void triquark_id_get(stringbuf str, const quark_no i1, const quark_no i2,\
					 const quark_no i3);

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