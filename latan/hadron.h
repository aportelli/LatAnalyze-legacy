#ifndef LATAN_HADRON_H_
#define LATAN_HADRON_H_

#include <latan/globals.h>
#include <latan/mat.h>

#define NOMIX 0
#define SUM 1
#define MEAN 2

#define EVEN 0
#define ODD 1

#ifndef MAXPROP
#define MAXPROP 4
#endif
#ifndef MAXQUARKST
#define MAXQUARKST 4
#endif
#ifndef MAXISO
#define MAXISO 4
#endif

__BEGIN_DECLS

typedef struct
{
	stringbuf name;
	stringbuf channel[MAXPROP];
	stringbuf quarkst[MAXQUARKST];
	int chmix;
	int stmix;
	int parity;
}* hadron;

typedef struct
{
	stringbuf name;
	size_t isodim;
	const hadron* had[MAXISO];
} isohad;

int hadron_getnt(const hadron h, const int source, const int sink,\
				 const stringbuf manfname);
latan_errno hadron_prop(mat* prop, const hadron h, const int source,\
						const int sink, const stringbuf manfname);

__END_DECLS

#endif