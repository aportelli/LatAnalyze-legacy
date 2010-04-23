#ifndef LATAN_RAND_H_
#define LATAN_RAND_H_

#include <latan/globals.h>

#define RLXG_STATE_SIZE 105

__BEGIN_DECLS

typedef int randgen_state[RLXG_STATE_SIZE];

void randgen_init(const int seed);
void randgen_timeinit(void);
void randgen_set_state(randgen_state state);
void randgen_get_state(randgen_state state);
double rand_u(double a, double b);
unsigned int rand_ud(const unsigned int n);
double rand_n(const double mean, const double sigma);

__END_DECLS

#endif