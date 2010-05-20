#include <stdlib.h>
#include <stdio.h>
#include "../config.h"
#include <latan/latan_hadron.h>

#define YES_NO(b) (((b) == 0) ? ("no") : ("yes"))
#define VERB_STR ((latan_get_verb() == QUIET) ? ("quiet") :\
					((latan_get_verb() == VERB) ? ("verbose") : ("debug")))
#define MIN_STR ((latan_get_minimize_lib() == GSL) ? ("GSL") : ("MINUIT"))
#define SEP printf("--------------------------------\n")

#ifdef HAVE_LIBCBLAS
#define LIBCBLAS 1
#else
#define LIBCBLAS 0
#endif
#ifdef HAVE_MINUIT2
#define MINUIT2 1
#else
#define MINUIT2 0
#endif
#ifdef HAVE_SSE
#define SSE 1
#else
#define SSE 0
#endif

int main(void)
{
	stringbuf prop_mark, prop_idfmt;
	
	latan_get_prop_mark(prop_mark);
	latan_get_prop_idfmt(prop_idfmt);
	
	printf("\n");
	printf("%s v%s\n",PACKAGE_NAME,PACKAGE_VERSION);
	SEP;
	printf("COMPILATION INFORMATION\n");
	SEP;
	printf("C   compiler              : %s v%s\n",C_COMP_VENDOR,GCC_VERSION);
	printf("C++ compiler              : %s v%s\n",CXX_COMP_VENDOR,GXX_VERSION);
	printf("CBLAS library outside GSL : %s\n",YES_NO(LIBCBLAS));
	printf("MINUIT library support    : %s\n",YES_NO(MINUIT2));
	printf("SSE optimizations         : %s\n",YES_NO(SSE));
	SEP;
	printf("DEFAULT OPTIONS\n");
	SEP;
	printf("verbosity level           : %s\n",VERB_STR);
	printf("minimization library      : %s\n",MIN_STR);
	printf("propagator mark           : %s\n",prop_mark);
	printf("propagator name format    : %s\n",prop_idfmt);
	printf("\n");

	return EXIT_SUCCESS;
}