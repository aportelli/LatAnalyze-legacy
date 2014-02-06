if CC_GNU
    COM_CFLAGS = -Wall -Wmissing-prototypes -Wstrict-prototypes			\
        -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align               \
        -Wwrite-strings -Wnested-externs -W -fPIC -fno-common -O2		\
        -g -ansi -pedantic
else
if CC_INTEL
    COM_CFLAGS = -Wall -Wmissing-prototypes -Wstrict-prototypes         \
        -Wshadow -Wpointer-arith -Wcast-qual -Wwrite-strings -fPIC      \
        -fno-common -O2 -g -ansi
endif
endif

bin_PROGRAMS =              \
	latan_create_rg_state   \
	latan_get_mass          \
	latan_info              \
	latan_mat_bench         \
	latan_rs_add            \
	latan_rs_adds           \
	latan_rs_div            \
	latan_rs_mul            \
	latan_rs_muls           \
    latan_rs_read           \
	latan_rs_sub            \
    latan_rs_sum            \
	latan_rs_sqrt

latan_create_rg_state_SOURCES = latan_create_rg_state.c
latan_create_rg_state_CFLAGS  = $(COM_CFLAGS)
latan_create_rg_state_LDFLAGS = -llatan -L../latan/.libs

latan_get_mass_SOURCES = latan_get_mass.c
latan_get_mass_CFLAGS  = $(COM_CFLAGS)
latan_get_mass_LDFLAGS = -llatan -L../latan/.libs

latan_info_SOURCES = latan_info.c
latan_info_CFLAGS  = $(COM_CFLAGS)
latan_info_LDFLAGS = -llatan -L../latan/.libs

latan_mat_bench_SOURCES = latan_mat_bench.c
latan_mat_bench_CFLAGS  = $(COM_CFLAGS)
latan_mat_bench_LDFLAGS = -llatan -L../latan/.libs

latan_rs_add_SOURCES = latan_rs_binop.c
latan_rs_add_CFLAGS  = $(COM_CFLAGS) -DBINOP=mat_add
latan_rs_add_LDFLAGS = -llatan -L../latan/.libs

latan_rs_adds_SOURCES = latan_rs_binops.c
latan_rs_adds_CFLAGS  = $(COM_CFLAGS) -DBINOPS=mat_adds
latan_rs_adds_LDFLAGS = -llatan -L../latan/.libs

latan_rs_div_SOURCES = latan_rs_binop.c
latan_rs_div_CFLAGS  = $(COM_CFLAGS) -DBINOP=mat_divp
latan_rs_div_LDFLAGS = -llatan -L../latan/.libs

latan_rs_mul_SOURCES = latan_rs_binop.c
latan_rs_mul_CFLAGS  = $(COM_CFLAGS) -DBINOP=mat_mulp
latan_rs_mul_LDFLAGS = -llatan -L../latan/.libs

latan_rs_muls_SOURCES = latan_rs_binops.c
latan_rs_muls_CFLAGS  = $(COM_CFLAGS) -DBINOPS=mat_muls
latan_rs_muls_LDFLAGS = -llatan -L../latan/.libs

latan_rs_read_SOURCES = latan_rs_read.c
latan_rs_read_CFLAGS  = $(COM_CFLAGS)
latan_rs_read_LDFLAGS = -llatan -L../latan/.libs

latan_rs_sub_SOURCES = latan_rs_binop.c
latan_rs_sub_CFLAGS  = $(COM_CFLAGS) -DBINOP=mat_sub
latan_rs_sub_LDFLAGS = -llatan -L../latan/.libs

latan_rs_sum_SOURCES = latan_rs_unop.c
latan_rs_sum_CFLAGS  = $(COM_CFLAGS) -DUNOP=mat_sum
latan_rs_sum_LDFLAGS = -llatan -L../latan/.libs

latan_rs_sqrt_SOURCES = latan_rs_unop.c
latan_rs_sqrt_CFLAGS  = $(COM_CFLAGS) -DUNOP=mat_sqrt
latan_rs_sqrt_LDFLAGS = -llatan -L../latan/.libs

ACLOCAL_AMFLAGS = -I .buildutils/m4