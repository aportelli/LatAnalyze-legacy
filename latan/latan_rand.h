#ifndef LATAN_RAND_H_
#define LATAN_RAND_H_

#include <latan/latan_globals.h>
#include <latan/latan_statistics.h>

__BEGIN_DECLS

void randgen_init(const int seed);
void randgen_init_from_time(void);
latan_errno randgen_set_state_from_rs_sample(rs_sample s);
void randgen_set_state(randgen_state state);
void randgen_get_state(randgen_state state);
double rand_u(double a, double b);
unsigned int rand_ud(const unsigned int n);
double rand_n(const double mean, const double sigma);

__END_DECLS

#endif