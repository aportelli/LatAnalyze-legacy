# Process this file with autoconf to produce a configure script.

# Initialization
AC_PREREQ([2.64])
AC_INIT([LatAnalyze],[2.0rc1],[antonin.portelli@me.com],[latan])
AC_CONFIG_AUX_DIR([.buildutils])
AC_CONFIG_SRCDIR([latan/latan_globals.c])
AC_CONFIG_SRCDIR([utils/latan_info.c])
AC_CONFIG_SRCDIR([examples/ex_mat.c])
AC_CONFIG_SRCDIR([doc/latanalyze.texi])
AC_CONFIG_MACRO_DIR([.buildutils/m4])
AM_INIT_AUTOMAKE([-Wall -Werror])
AM_SILENT_RULES([yes])

AC_CONFIG_HEADERS([config.h])

# Checks for programs.
m4_ifdef([AM_PROG_AR],[AM_PROG_AR])
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CXX
AC_PROG_LIBTOOL

LT_INIT

AC_ARG_ENABLE([Accelerate],
	[AS_HELP_STRING([--enable-Accelerate],
		[link with CBLAS implentation in Accelerate.framework on Apple platforms])],
	[],
	[]
)
AC_ARG_ENABLE([SSE],
	[AS_HELP_STRING([--enable-SSE],
		[compiles SSE version of ranlux random generator])],
	[AC_DEFINE([HAVE_SSE],
				[1],
				[Define to 1 if your CPU support SSE instructions.])],
	[]
)

# Configure parameters non-standard install prefix and name of libs
AC_ARG_WITH([gsl],
    [AS_HELP_STRING([--with-gsl=prefix],
		[try this for a non-standard install prefix of the GSL library])],
    [AM_CFLAGS="$AM_CFLAGS -I$with_gsl/include"]
	[AM_CXXFLAGS="$AM_CXXFLAGS -I$with_gsl/include"]
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
AC_ARG_WITH([cblas-name],
    [AS_HELP_STRING([--with-cblas-name=name],
		[use this to specify the filename (without 'lib') of CBLAS library (default : 'cblas')])],
    [],[with_cblas_name="cblas"]
)
AC_ARG_WITH([Minuit2],
    [AS_HELP_STRING([--with-Minuit2=prefix],
		[try this for a non-standard install prefix of the Minuit2 library])],
    [AM_CFLAGS="$AM_CFLAGS -I$with_Minuit2/include"]
    [AM_CXXFLAGS="$AM_CXXFLAGS -I$with_Minuit2/include"]
	[AM_LDFLAGS="$AM_LDFLAGS -L$with_Minuit2/lib"],
	[]
)
AC_ARG_WITH([xml2],
    [AS_HELP_STRING([--with-xml2=prefix],
		[try this for a non-standard install prefix of the xml2 library])],
    [AM_CFLAGS="$AM_CFLAGS -I$with_xml2/include -I$with_xml2/include/libxml2"]
	[AM_LDFLAGS="$AM_LDFLAGS -L$with_xml2/lib"],
	[]
)

CFLAGS="$AM_CFLAGS $CFLAGS"
CXXFLAGS="$AM_CXXFLAGS $CXXFLAGS"
LDFLAGS="$AM_LDFLAGS $LDFLAGS"

# Get compilers informations
AX_COMPILER_VENDOR
AC_DEFINE_UNQUOTED([C_COMP_VENDOR],["$ax_cv_c_compiler_vendor"],
			[vendor of C compiler that will compile the code])
AM_CONDITIONAL([CC_GNU],[test $ax_cv_c_compiler_vendor = "gnu"])
AM_CONDITIONAL([CC_INTEL],[test $ax_cv_c_compiler_vendor = "intel"])
AC_LANG([C++])
AX_COMPILER_VENDOR
AC_DEFINE_UNQUOTED([CXX_COMP_VENDOR],["$ax_cv_cxx_compiler_vendor"],
			[vendor of C++ compiler that will compile the code])
AM_CONDITIONAL([CXX_GNU],[test $ax_cv_cxx_compiler_vendor = "gnu"])
AM_CONDITIONAL([CXX_INTEL],[test $ax_cv_cxx_compiler_vendor = "intel"])
AC_LANG([C])
AX_GCC_VERSION
AC_DEFINE_UNQUOTED([GCC_VERSION],["$GCC_VERSION"],
			[version of gcc that will compile the code])
AX_GXX_VERSION
AC_DEFINE_UNQUOTED([GXX_VERSION],["$GXX_VERSION"],
			[version of g++ that will compile the code])

# Check for OpenMP
AC_OPENMP
AC_LANG([C++])
AC_OPENMP
AC_LANG([C])
AM_CFLAGS="$AM_CFLAGS $OPENMP_CFLAGS"
CFLAGS="$AM_CFLAGS $CFLAGS"
AM_CXXFLAGS="$AM_CXXFLAGS $OPENMP_CXXFLAGS"
CXXFLAGS="$AM_CXXFLAGS $CXXFLAGS"

# Checks for libraries.
AC_CHECK_LIB([m],[cos],[],[AC_MSG_ERROR([libm library not found])])
AS_IF([test "x$enable_Accelerate" = "xyes"],
	[AC_CHECK_FRAMEWORK([Accelerate],[cblas_dgemm],
		[AC_DEFINE_UNQUOTED([CBLAS_NAME],["Accelerate framework"],[name of CBLAS library.])]
		[AC_DEFINE_UNQUOTED([HAVE_LIBCBLAS],[1],[Define to 1 if you have a CBLAS library outside the GSL.])]
		[LIBS="-framework Accelerate $LIBS"],
		[AC_MSG_ERROR([Accelerate framework not found])])],
	[AC_CHECK_LIB([$with_cblas_name],[cblas_dgemm],
		[AC_DEFINE_UNQUOTED([CBLAS_NAME],["lib$with_cblas_name"],[name of CBLAS library.])]
		[AC_DEFINE_UNQUOTED([HAVE_LIBCBLAS],[1],[Define to 1 if you have a CBLAS library outside the GSL.])]
		[LIBS="-l$with_cblas_name $LIBS"],
		[AC_CHECK_LIB([gslcblas],[cblas_dgemm],[],
				[AC_MSG_ERROR([CBLAS library not found])])])])
		
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],[],[AC_MSG_ERROR([GSL library not found])])
AC_CHECK_LIB([xml2],[xmlFree],[AM_CFLAGS="$AM_CFLAGS `xml2-config --cflags`"],[])
AC_LANG([C++])
AC_CHECK_LIB([stdc++],[main],[LIBS="-lstdc++ $LIBS"],[AC_MSG_ERROR([libstdc++ library not found])])
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
AC_CHECK_HEADERS([unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([sqrt],[],[AC_MSG_ERROR([sqrt function not found])])
AC_CHECK_FUNCS([acosh])
AC_CHECK_FUNCS([strtok_r])

AC_SUBST([LIBS])
AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_LDFLAGS])

AC_CONFIG_FILES([Makefile latan/Makefile utils/Makefile examples/Makefile \
                 doc/Makefile])
AC_OUTPUT
