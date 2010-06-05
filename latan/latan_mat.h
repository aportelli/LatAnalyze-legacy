#ifndef LATAN_MAT_H_
#define LATAN_MAT_H_

#include <latan/latan_globals.h>
#include <gsl/gsl_matrix.h>

__BEGIN_DECLS

/* definition of the matrix type */
typedef gsl_matrix* mat;

/* loop */
#define FOR_VAL(m,i,j)\
for (i=0;i<nrow(m);i++)\
for (j=0;j<ncol(m);j++)\

/* functions */
/** allocation **/
mat mat_create(const size_t init_nrow, const size_t init_ncol);
#define mat_create_from_dim(n) mat_create(nrow(n),ncol(n))
#define mat_create_from_trdim(n) mat_create(ncol(n),nrow(n))
mat mat_create_from_mat(const mat n);
mat mat_create_from_ar(const double* ar, const size_t init_nrow,\
					   const size_t init_ncol);
mat* mat_ar_create(const size_t nmat, const size_t init_nrow,\
				   const size_t init_ncol);
#define mat_ar_create_from_dim(nmat,n) mat_ar_create(nmat,nrow(n),ncol(n))
void mat_destroy(mat m);
void mat_ar_destroy(mat* m, const size_t nmat);

/** access **/
size_t nrow(const mat m);
size_t ncol(const mat m);
double mat_get(const mat m, const size_t i, const size_t j);
void mat_set(mat m, const size_t i, const size_t j, const double val);
latan_errno mat_set_subm(mat m, const mat n, const size_t k1, const size_t l1, \
						 const size_t k2, const size_t l2);
latan_errno mat_set_diag(mat m, const mat diag);
latan_errno mat_set_step(mat m, const double x0, const double step);
latan_errno mat_set_from_ar(mat m, const double* ar);
#define mat_inc(m,i,j,val) mat_set(m,i,j,mat_get(m,i,j)+val)
#define mat_pp(m,i,j) mat_inc(m,i,j,1.0)

/** tests **/
bool mat_issamedim(mat m, mat n);
bool mat_issquare(mat m);

/** operations **/
void mat_zero(mat m);
void mat_cst(mat m, const double x);
void mat_rand_u(mat m, const double a, const double b);
void mat_id(mat m);
latan_errno mat_cp(mat m, const mat n);
latan_errno mat_cp_subm(mat m, const mat n, const size_t k1, const size_t l1, \
						const size_t k2, const size_t l2);
latan_errno mat_eqadd(mat m, const mat n);
latan_errno mat_add(mat m, const mat n, const mat o);
latan_errno mat_eqsub(mat m, const mat n);
latan_errno mat_sub(mat m, const mat n, const mat o);
#define mat_eqmul_l(m,n) mat_mul(m,n,m);
#define mat_eqmul_r(m,n) mat_mul(m,m,n);
#define mat_eqmul_l_t(m,n) mat_mul_tn(m,n,m);
#define mat_eqmul_r_t(m,n) mat_mul_nt(m,m,n);
latan_errno mat_mul_nn(mat m, const mat n, const mat o);
latan_errno mat_mul_nt(mat m, const mat n, const mat o);
latan_errno mat_mul_tn(mat m, const mat n, const mat o);
latan_errno mat_mul_tt(mat m, const mat n, const mat o);
#define mat_mul(m,n,o) mat_mul_nn(m,n,o);
latan_errno mat_eqtranspose(mat m);
latan_errno mat_transpose(mat m, const mat n);
latan_errno mat_eqmulp(mat m, const mat n);
latan_errno mat_mulp(mat m, const mat n, const mat o);
latan_errno mat_eqmuls(mat m, const double s);
latan_errno mat_muls(mat m, const mat n, const double s);
latan_errno mat_eqdivp(mat m, const mat n);
latan_errno mat_divp(mat m, const mat n, const mat o);
#define mat_eqabs(m) mat_abs(m,m)
latan_errno mat_abs(mat m, const mat n);
#define mat_eqsqrt(m) mat_sqrt(m,m)
latan_errno mat_sqrt(mat m, const mat n);

/** linear algebra **/
#define mat_eqinv(m) mat_inv(m,m);
latan_errno mat_inv(mat m, const mat n);
latan_errno mat_sym_diagonalize(mat eigval, mat eigvec, const mat m);

__END_DECLS

#endif