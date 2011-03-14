#include <stdlib.h>
#include <stdio.h>
#include "../config.h"
#include <latan/latan_hadron.h>
#include <latan/latan_io.h>
#include <latan/latan_minimizer.h>

#define YES_NO(b) (((b) == 0) ? ("no") : ("yes"))
#define VERB_STR ((latan_get_verb() == QUIET) ? ("quiet") :\
                    ((latan_get_verb() == VERB) ? ("verbose") : ("debug")))
#define FMT_STR ((io_get_fmt() == IO_ASCII) ? ("ASCII") : ("GSL"))
#define MIN_STR ((minimizer_get_lib() == GSL) ? ("GSL") : ("MINUIT"))
#define SEP printf("---------------------------------------------\n")

#if defined(CBLAS_NAME)||defined(HAVE_FRAMEWORK_ACCELERATE)
#define LIBCBLAS 1
#else
#define LIBCBLAS 0
#endif
#ifdef HAVE_MINUIT2
#define LIBMINUIT2 1
#else
#define LIBMINUIT2 0
#endif
#ifdef HAVE_LIBXML2
#define LIBXML2 1
#else
#define LIBXML2 0
#endif
#ifdef HAVE_SSE
#define SSE 1
#else
#define SSE 0
#endif
#ifdef _OPENMP
#define OMP 1
#else
#define OMP 0
#endif

int main(void)
{    
    printf("\n");
    printf("%s v%s\n",PACKAGE_NAME,PACKAGE_VERSION);
    SEP;
    printf("COMPILATION INFORMATION\n");
    SEP;
    printf("C   compiler              : %s v%s\n",C_COMP_VENDOR,GCC_VERSION);
    printf("C++ compiler              : %s v%s\n",CXX_COMP_VENDOR,GXX_VERSION);
#if (LIBCBLAS == 1)
    printf("CBLAS library outside GSL : yes (%s)\n",CBLAS_NAME);
#else
    printf("CBLAS library outside GSL : no\n");
#endif
    printf("MINUIT library support    : %s\n",YES_NO(LIBMINUIT2));
    printf("XML datafile support      : %s\n",YES_NO(LIBXML2));
    printf("SSE optimizations         : %s\n",YES_NO(SSE));
    printf("OpenMP support            : %s\n",YES_NO(OMP));
    printf("\n");
    SEP;
    printf("DEFAULT GLOBAL OPTIONS\n");
    SEP;
    printf("verbosity level           : %s\n",VERB_STR);
    printf("\n");
    SEP;
    printf("DEFAULT I/O OPTIONS\n");
    SEP;
    printf("datafile format           : %s\n",FMT_STR);
    printf("\n");
    SEP;
    printf("DEFAULT MINIMIZER OPTIONS\n");
    SEP;
    printf("library                   : %s\n",MIN_STR);
    printf("maximum iteration         : %u\n",minimizer_get_max_iteration());
    printf("\n");
    return EXIT_SUCCESS;
}
