if CC_GNU
    COM_CFLAGS = -Wall -Wmissing-prototypes -Wstrict-prototypes			\
        -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align               \
        -Wwrite-strings -Wnested-externs -W -fPIC -fno-common -ansi
else
if CC_INTEL
    COM_CFLAGS = -Wall -Wmissing-prototypes -Wstrict-prototypes         \
        -Wshadow -Wpointer-arith -Wcast-qual -Wwrite-strings -fPIC      \
        -fno-common -ansi
endif
endif

if CXX_GNU
    COM_CXXFLAGS = -Wall -Wshadow -Wpointer-arith -Wcast-qual           \
        -Wcast-align -Wwrite-strings -fPIC -fno-common
else
if CXX_INTEL
    COM_CXXFLAGS = -Wall -Wshadow -Wpointer-arith -Wcast-qual           \
        -Wwrite-strings -fPIC -fno-common
endif
endif

lib_LTLIBRARIES = liblatan.la
	
liblatan_la_SOURCES =       \
	latan_blas.h            \
	latan_blas.c            \
	latan_error.c           \
	latan_fit.c             \
	latan_globals.c         \
	latan_includes.h        \
	latan_io.c              \
	latan_io_ascii.h        \
	latan_io_ascii.c        \
	latan_io_xml.h          \
	latan_io_xml.c          \
	latan_mass.c            \
	latan_mat.c             \
	latan_math.c            \
	latan_minimizer.c       \
	latan_min_gsl.h         \
	latan_min_gsl.c         \
	latan_min_minuit2.h     \
	latan_min_minuit2.cpp   \
	latan_models.c          \
	latan_plot.c            \
	latan_rand.c            \
	latan_statistics.c      \
    latan_tabfunc.c         \
	latan_xml.h             \
	latan_xml.c             \
	../config.h
liblatan_ladir = $(includedir)/latan
liblatan_la_HEADERS =       \
	latan_error.h           \
	latan_fit.h             \
	latan_globals.h         \
	latan_io.h              \
	latan_mass.h            \
	latan_mat.h             \
	latan_math.h            \
	latan_minimizer.h       \
	latan_models.h          \
	latan_nunits.h          \
	latan_plot.h            \
	latan_statistics.h      \
    latan_tabfunc.h         \
	latan_rand.h			
liblatan_la_CFLAGS = $(COM_CFLAGS)
liblatan_la_CXXFLAGS = $(COM_CXXFLAGS)

ACLOCAL_AMFLAGS = -I .buildutils/m4