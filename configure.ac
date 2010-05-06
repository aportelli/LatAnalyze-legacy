# Process this file with autoconf to produce a configure script.

# Initialization
AC_PREREQ([2.64])
AC_INIT([LatAnalyze],[2.0alpha],[antonin.portelli@gmail.com],[latan])
AC_CONFIG_AUX_DIR([.buildutils])
AC_CONFIG_SRCDIR([latan/globals.c])
AC_CONFIG_MACRO_DIR([.buildutils/m4])
AM_INIT_AUTOMAKE([-Wall -Werror])

AM_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CXX
AC_PROG_LIBTOOL

LT_INIT

# Configure parameters non-standard install prefix of libs
AC_ARG_WITH([gsl],
    [AS_HELP_STRING([--with-gsl=prefix],
	[try this for a non-standard install prefix of the GSL library])],
    [AM_CFLAGS="$AM_CFLAGS -I$with_gsl/include"]
	[AM_LDFLAGS="$AM_LDFLAGS -L$with_gsl/lib"],
	[]
)
AC_ARG_WITH([cblas],
    [AS_HELP_STRING([--with-cblas=prefix],
	[try this for a non-standard install prefix of the CBLAS library])],
    [AM_CFLAGS="$AM_CFLAGS -I$with_cblas/include"]
	[AM_LDFLAGS="$AM_LDFLAGS -L$with_cblas/lib"],
	[]
)
AC_ARG_WITH([Minuit2],
    [AS_HELP_STRING([--with-Minuit2=prefix],
	[try this for a non-standard install prefix of the Minuit2 library])],
    [AM_CFLAGS="$AM_CFLAGS -I$with_gsl/include"]
	[AM_LDFLAGS="$AM_LDFLAGS -L$with_gsl/lib"],
	[]
)

CFLAGS="$AM_CFLAGS $CFLAGS"
LDFLAGS="$AM_LDFLAGS $LDFLAGS"

# Checks for libraries.
AC_CHECK_LIB([m],[cos],[],[AC_MSG_ERROR([libm library not found])])
AC_CHECK_LIB([cblas],[cblas_dgemm],[],
	[AC_CHECK_LIB([gslcblas],[cblas_dgemm],[],
		[AC_MSG_ERROR([CBLAS library not found])])])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],[],[AC_MSG_ERROR([GSL library not found])])
AC_LANG([C++])
SAVED_LDFLAGS=$LDFLAGS
LDFLAGS="$LDFLAGS -lMinuit2"
AC_MSG_CHECKING([for ROOT::Minuit2::BasicMinimumError in -lMinuit2]);
AC_LINK_IFELSE(
	[AC_LANG_PROGRAM([#include <Minuit2/BasicMinimumError.h>],
	[ROOT::Minuit2::BasicMinimumError dummy(0)])],
	[LIBS="$LIBS -lMinuit2"]
	[AC_DEFINE([HAVE_MINUIT2],
				[1],
				[Define to 1 if you have the `Minuit2' library (-lMinuit2).])]
	[AC_MSG_RESULT([yes])],
	[AC_MSG_RESULT([no])])
LDFLAGS=$SAVED_LDFLAGS
AC_LANG([C])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([sqrt],[],[AC_MSG_ERROR([sqrt function not found])])
AC_CHECK_FUNCS([acosh])

AC_SUBST([LIBS])
AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_LDFLAGS])

AC_CONFIG_FILES([Makefile latan/Makefile doc/Makefile])
AC_OUTPUT
