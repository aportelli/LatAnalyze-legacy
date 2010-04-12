/* Main Doxygen page */
/*!
 @mainpage LatAnalyze - C library for lattice QCD results analysis
 @author <A HREF=mailto:antonin.portelli@gmail.com>Antonin Portelli</A>
 
 @section seci_nst 1 - Installation
 To compile LatAnalyze you need only the <A HREF=http://gcc.gnu.org/>GNU C compiler</A> and C standard library. If you want to compile the documentation,
 you will need <A HREF=http://www.doxygen.nl>Doxygen</A> and a <A HREF=http://www.latex-project.org/>LaTeX</A> distribution. 
 <A HREF=http://www.gnuplot.info/>gnuplot</A> is not needed to compile the library, but it must be present at run time if you use functions from plot.h
 header. Once this environment is set :
 -# download the source of the library at http://dakant.hd.free.fr/latanalyze/latan-1.1.tar.gz :<BR>
 <TT>$ wget http://dakant.hd.free.fr/latanalyze/latan-1.1.tar.gz</TT>
 -# untar the <TT>latan-*.tar.gz</TT> archive :<BR>
 <TT>$ tar -xvzf latan-*.tar.gz</TT>
 -# go in the created directory :<BR>
 <TT>$ cd latan-*</TT>
 -# check that the head of <TT>Makefile</TT> match your system configuration
 -# compile the library :<BR>
 <TT>$ make lib</TT>
 -# finally install the library :<BR>
 <TT>$ make install</TT>
 -# (optional) compile the documentation, it will be generated in three formats, HTML \htmlonly (this document) \endhtmlonly in <TT>doc/html</TT>, PDF
 \latexonly (this document) \endlatexonly in <TT>doc/latex</TT> and Unix man \manonly (this document) \endmanonly in <TT>doc/man</TT> :<BR>
 <TT>$ make doc</TT>
 -# (optional) compile the examples given with the documentation, it's a good way to test that the library is working, the compiled example
 executables will be in <TT>doc/examples</TT> :<BR>
 <TT>$ make examples</TT>
 
 You can uninstall LatAnalyze library with the <TT>make uninstall</TT> command.
 
 @section sec_usage 2 - Usage
 Include the needed LatAnalyze headers in your C code with include macros like :
 \code #include <latan/xxx.h> \endcode
 where <TT>xxx.h</TT> is the header you want to use. Eventualy use gcc <TT>-I</TT> option to indicate the LatAnalyze header files location at compilation. 
 You also need to link your final executable with LatAnalyze with the gcc <TT>-llatan</TT> option, eventualy use gcc <TT>-L</TT> option to indicate 
 the LatAnalyze library location for linking.
 @section sec_history 3 - History
 @subsection v2_0a v2.0a (2010-??-??)
*/
/*!
 @file mat.h
 @brief Matrix of double precision float numbers structure implementation.
 
 See io.h for reading and saving matrices on files.
 @author <A HREF=mailto:antonin.portelli@gmail.com>Antonin Portelli</A> (<A HREF=http://www.cpt.univ-mrs.fr/>CPT</A>)
*/
/*!
 @example ex_mat.c Manipulating mat variables.
*/
#ifndef LATAN_MAT_H_
#define LATAN_MAT_H_

#include <latan/globals.h>
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
mat* mat_create_ar(const size_t nmat, const size_t init_nrow,\
				   const size_t init_ncol);
#define mat_create_ar_from_dim(nmat,n) mat_create_ar(nmat,nrow(n),ncol(n))
void mat_destroy(mat m);
void mat_destroy_ar(mat* m, const size_t nmat);

/** access **/
size_t nrow(const mat m);
size_t ncol(const mat m);
double mat_get(const mat m, const size_t i, const size_t j);
void mat_set(mat m, const size_t i, const size_t j, const double val);
latan_errno mat_set_subm(mat m, const mat n, const size_t k1, const size_t l1, \
						 const size_t k2, const size_t l2);
latan_errno mat_set_from_ar(mat m, const double* ar);
#define mat_inc(m,i,j,val) mat_set(m,i,j,mat_get(m,i,j)+val)
#define mat_pp(m,i,j) mat_inc(m,i,j,1.0)

/** tests **/
/*!
 @fn BOOL mati_ssamedim(mat m, mat n)
 @brief Compare dimensions of two matrices.
 
 @return #TRUE is \b m got the same dimensions than \b n, #FALSE else.
*/
bool mat_issamedim(mat m, mat n);
/*!
 @fn BOOL mati_ssquare(mat m)
 @brief Check if a matrix is square.
 
 @return #TRUE if \b m is square, #FALSE else.
*/
bool mat_issquare(mat m);

/** operations **/
void mat_zero(mat m);
void mat_cst(mat m, const double x);
void mat_rand_u(mat m, const double a, const double b);
void mat_id(mat m);
/*!
 @fn void mat_cp(mat m, const mat n)
 @brief Copy the content of the matrix \b n in \b m.
 
 @warning \b m and \b n must have the same dimensions, fatal error will be generated if it is not the case.
 */
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
/*!
 @fn void mat_mul(mat m, const mat n, const mat o)
 @brief Store in \b m the matrix product of \b n and \b o.
 @see #mat_eqmul
 
 @warning \b NCOL(n) and \b NROW(o) must be equals, \b m must have the same number of rows than \b n and the same
 number of columns than \b o. Fatal error will be generated if precedent conditions are not satisfied.
 */
latan_errno mat_mul_nn(mat m, const mat n, const mat o);
latan_errno mat_mul_nt(mat m, const mat n, const mat o);
latan_errno mat_mul_tn(mat m, const mat n, const mat o);
latan_errno mat_mul_tt(mat m, const mat n, const mat o);
#define mat_mul(m,n,o) mat_mul_nn(m,n,o);
latan_errno mat_eqtranspose(mat m);
latan_errno mat_transpose(mat m, const mat n);
#define mat_eqinv(m) mat_inv(m,m);
latan_errno mat_inv(mat m, const mat n);
/*!
 @fn void mat_eqmulp(mat m, const mat n)
 @brief Multiply \b m by \b n "coefficient by coefficient".
 @see #mat_mulp
 
 @warning \b m and \b n must have the same dimensions, fatal error will be generated if it is not the case.
*/
latan_errno mat_eqmulp(mat m, const mat n);
latan_errno mat_mulp(mat m, const mat n, const mat o);
latan_errno mat_eqmuls(mat m, const double s);
latan_errno mat_muls(mat m, const mat n, const double s);
latan_errno mat_eqdivp(mat m, const mat n);
latan_errno mat_divp(mat m, const mat n, const mat o);
#define mat_eqabs(m) mat_abs(m,m)
latan_errno mat_abs(mat m, const mat n);
/*!
 @def mat_eqsqrt(m)
 @brief Compute the "coefficient by coefficient" square root of \b m.
 @see ::mat_sqrt
 */
#define mat_eqsqrt(m) mat_sqrt(m,m)
/*!
 @fn void mat_sqrt(mat m, const mat n)
 @brief Store in \b m the "coefficient by coefficient" square root of \b n.
 @see #mat_eqsqrt
 
 @warning \b m and \b n must have the same dimensions, fatal error will be generated if it is not the case.
*/
latan_errno mat_sqrt(mat m, const mat n);


__END_DECLS

#endif