#include <stdlib.h>
#include <stdio.h>
#include "../config.h"
#include <latan/latan_hadron.h>
#include <latan/latan_io.h>
#include <latan/latan_minimizer.h>

#define YES_NO(b) (((b) == 0) ? ("no") : ("yes"))
#define VERB_STR ((latan_get_verb() == QUIET) ? ("quiet") :\
					((latan_get_verb() == VERB) ? ("verbose") : ("debug")))
#define MIN_STR ((minimizer_get_lib() == GSL) ? ("GSL") : ("MINUIT"))
#define SEP printf("---------------------------------------------\n")

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
	
	io_get_prop_mark(prop_mark);
	io_get_prop_idfmt(prop_idfmt);
	
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
	printf("\n");
	SEP;
	printf("DEFAULT GLOBAL OPTIONS\n");
	SEP;
	printf("verbosity level           : %s\n",VERB_STR);
	printf("\n");
	SEP;
	printf("DEFAULT I/O OPTIONS\n");
	SEP;
	printf("propagator mark           : \"%s\"\n",prop_mark);
	printf("propagator name format *   : \"%s\"\n",prop_idfmt);
	printf("\n");
	SEP;
	printf("DEFAULT MINIMIZER OPTIONS\n");
	SEP;
	printf("library                   : %s\n",MIN_STR);
	printf("maximum iteration         : %u\n",minimizer_get_max_iteration());
	printf("\n");
	return EXIT_SUCCESS;
}