COM_CFLAGS = -Wall -Wmissing-prototypes -Wstrict-prototypes					\
	-Wshadow -Wpointer-arith -Wcast-qual -Wcast-align						\
	-Wwrite-strings -Wnested-externs -W -fshort-enums -fno-common -O2		\
	-g -ansi -pedantic

bin_PROGRAMS = latan_info

latan_info_SOURCES = latan_info.c
latan_info_CFLAGS = $(COM_CFLAGS)
latan_info_LDFLAGS = -llatan -L../latan/.libs

ACLOCAL_AMFLAGS = -I .buildutils/m4