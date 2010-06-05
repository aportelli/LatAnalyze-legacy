COM_CFLAGS = -Wall -Wmissing-prototypes -Wstrict-prototypes					\
	-Wshadow -Wpointer-arith -Wcast-qual -Wcast-align						\
	-Wwrite-strings -Wnested-externs -W -fshort-enums -fno-common -O2		\
	-g -ansi -pedantic

bin_PROGRAMS =				\
	latan_create_rg_state	\
	latan_info				\
	latan_rs_add			\
	latan_rs_mul			\
	latan_rs_sub

latan_create_rg_state_SOURCES = latan_create_rg_state.c
latan_create_rg_state_CFLAGS = $(COM_CFLAGS)
latan_create_rg_state_LDFLAGS = -llatan -L../latan/.libs

latan_info_SOURCES = latan_info.c
latan_info_CFLAGS = $(COM_CFLAGS)
latan_info_LDFLAGS = -llatan -L../latan/.libs

latan_rs_add_SOURCES = latan_rs_add.c
latan_rs_add_CFLAGS = $(COM_CFLAGS)
latan_rs_add_LDFLAGS = -llatan -L../latan/.libs

latan_rs_mul_SOURCES = latan_rs_mul.c
latan_rs_mul_CFLAGS = $(COM_CFLAGS)
latan_rs_mul_LDFLAGS = -llatan -L../latan/.libs

latan_rs_sub_SOURCES = latan_rs_sub.c
latan_rs_sub_CFLAGS = $(COM_CFLAGS)
latan_rs_sub_LDFLAGS = -llatan -L../latan/.libs

ACLOCAL_AMFLAGS = -I .buildutils/m4