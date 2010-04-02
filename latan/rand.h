/*!
 @file rand.h
 @brief Random number generators based on the <A HREF=http://luscher.web.cern.ch/luscher/ranlux/>ranlux random number generator</A>
 by Martin Lüscher.

 @author <A HREF=mailto:antonin.portelli@gmail.com>Antonin Portelli</A> (<A HREF=http://www.cpt.univ-mrs.fr/>CPT</A>)
*/

#ifndef LATAN_RAND_H_
#define LATAN_RAND_H_

#include <latan/globals.h>

#define RLXG_STATE_SIZE 105

__BEGIN_DECLS

typedef int rand_gen_state[RLXG_STATE_SIZE];

/*!
 @fn void rand_init(const int seed)
 @brief Initialize the generator with a given seed.
 
 @param seed a integer between \f$1\f$ and \f$2^{31}\f$
*/
void rand_init(const int seed);
/*!
 @fn void rand_timeinit(void)
 @brief Initialize the generator with the number of seconds since the 1st January of 1970 12:00am as a seed.
*/
void rand_timeinit(void);
void rand_set_gen_state(rand_gen_state state);
void rand_get_gen_state(rand_gen_state state);
double rand_u(double a, double b);
/*!
 @fn int rand_ud(int n)
 @brief Random number generator with a discrete uniform distribution.
 
 @return a uniformly distributed random integer between \f$0\f$ and \b n.
*/
int rand_ud(const int n);

__END_DECLS

#endif