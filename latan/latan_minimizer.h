#ifndef LATAN_MINIMIZER_H_
#define LATAN_MINIMIZER_H_

#include <latan/latan_globals.h>

__BEGIN_DECLS

/* minimize lib flags */
typedef enum
{
	GSL    = 0,\
	MINUIT = 1
} minlib_no;

/* minization algorithms */
#define NMINALG 6
typedef enum
{
	GSL_GRAD_FR    = 0,\
	GSL_GRAD_PR    = 1,\
	GSL_VEC_BFGS   = 2,\
	GSL_SIMPLEX_NM = 3,\
	MIN_MIGRAD     = 4,\
	MIN_SIMPLEX    = 5
} minalg_no;

minalg_no minalg_no_get(const stringbuf m_id);
latan_errno minalg_id_get(stringbuf m_id, const minalg_no n);

/* minimizer options */
minlib_no minimizer_get_lib(void);
minalg_no minimizer_get_alg(void);
latan_errno minimizer_set_alg(minalg_no alg);
latan_errno minimizer_get_alg_name(stringbuf name);
unsigned int minimizer_get_max_iteration(void);
void minimizer_set_max_iteration(unsigned int max_iteration);

/* prototype of function to minimize */
typedef double min_func(const mat *x, void *param);

/* the minimizer */
latan_errno minimize(mat *x, double *f_min, min_func *f, void *param);

__END_DECLS

#endif