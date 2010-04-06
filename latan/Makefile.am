COM_CFLAGS = -Wall -Wmissing-prototypes -Wstrict-prototypes					\
	-Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align			\
	-Wwrite-strings -Wnested-externs -W -fPIC -fshort-enums -fno-common -O2	\
	-g -ansi -pedantic -Dinline=
	
COM_CXXFLAGS = -Wall -Wconversion -Wshadow -Wpointer-arith -Wcast-qual		\
	-Wcast-align -Wwrite-strings -W -fPIC -fshort-enums -fno-common -O2 -g	\
	-pedantic

lib_LTLIBRARIES = liblatan.la
liblatan_la_SOURCES =	\
	error.c				\
	globals.c			\
	io.c				\
	mat.c				\
	rand.c				\
	statistics.c		\
	../config.h
liblatan_ladir = $(includedir)/latan
liblatan_la_HEADERS =	\
	error.h				\
	globals.h			\
	io.h				\
	mat.h				\
	statistics.h		\
	rand.h
liblatan_la_CFLAGS = $(COM_CFLAGS)
liblatan_la_CXXFLAGS = $(COM_CXXFLAGS)