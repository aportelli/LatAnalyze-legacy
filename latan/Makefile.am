COM_CFLAGS = -Wall -Wmissing-prototypes -Wstrict-prototypes					\
	-Wshadow -Wpointer-arith -Wcast-qual -Wcast-align						\
	-Wwrite-strings -Wnested-externs -W -fPIC -fshort-enums -fno-common -O2	\
	-g -ansi -pedantic
	
COM_CXXFLAGS = -Wall -Wconversion -Wshadow -Wpointer-arith -Wcast-qual		\
	-Wcast-align -Wwrite-strings -W -fPIC -fshort-enums -fno-common -O2 -g	\
	-pedantic

noinst_LTLIBRARIES = liblatancublas.la
lib_LTLIBRARIES = liblatan.la

liblatancublas_la_SOURCES = \
	latan_mat_cublas.h		\
	latan_mat_cublas.c		\
	latan_globals.h			\
	latan_mat.h				\
	latan_error.h
	
liblatan_la_SOURCES =		\
	latan_error.c			\
	latan_fit.c				\
	latan_globals.c			\
	latan_hadron.c			\
	latan_includes.h		\
	latan_io.c				\
	latan_mass.c			\
	latan_mat.c				\
	latan_math.c			\
	latan_minimizer.c		\
	latan_min_gsl.h			\
	latan_min_gsl.c			\
	latan_min_minuit2.h		\
	latan_min_minuit2.cpp	\
	latan_models.c			\
	latan_plot.c			\
	latan_rand.c			\
	latan_statistics.c		\
	../config.h
liblatan_ladir = $(includedir)/latan
liblatan_la_HEADERS =		\
	latan_error.h			\
	latan_fit.h				\
	latan_globals.h			\
	latan_hadron.h			\
	latan_io.h				\
	latan_mass.h			\
	latan_mat.h				\
	latan_math.h			\
	latan_minimizer.h		\
	latan_models.h			\
	latan_nunits.h			\
	latan_plot.h			\
	latan_statistics.h		\
	latan_rand.h
liblatan_la_LIBADD = ./liblatancublas.la
liblatan_la_CFLAGS = $(COM_CFLAGS)
liblatan_la_CXXFLAGS = $(COM_CXXFLAGS)

ACLOCAL_AMFLAGS = -I .buildutils/m4