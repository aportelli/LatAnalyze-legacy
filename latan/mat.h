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
 d6rtdddfrdrf
 @section sec_history 3 - History
 @subsection v1_0 v1.0 (2009-11-09)
 First effective and tested version including initial features :
 \arg dynamically allocated matrix type with basic arithmetic operations (see mat.h)
 \arg I/O functions to read and save matrices in files (see io.h)
 \arg bootstrap error evaluation function using array of matrices as data (see bootstrap.h)
 @subsection v1_01 v1.01 (2009-11-10)
 \arg FIX : ::mat_meansig was calculating variance and not standard deviation
 \arg ADD : "coefficient by coefficient" square root operation on matrix added (see ::mat_sqrt and #mat_eqsqrt)
 @subsection v1_1 v1.1 (2009-11-14)
 \arg FIX : several memory leaks and allocation bugs
 \arg ADD : <A HREF=http://www.gnuplot.info/>gnuplot</A> interface based on Nicolas Devillard <A HREF=http://ndevilla.free.fr/gnuplot/>gnuplot C API</A>
 (see plot.h)
 \arg ADD : function to compute effective mass from propagator datas (see ::effmass)
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
void mat_create(mat* m, const size_t init_nrow, const size_t init_ncol);
void mat_create_from_mat(mat* m, const mat n);
void mat_create_ar(mat** m, const size_t nmat, const size_t init_nrow,\
				   const size_t init_ncol);
void mat_destroy(mat* m);
void mat_destroy_ar(mat** m, const size_t nmat);

/** access **/
size_t nrow(const mat m);
size_t ncol(const mat m);
double mat_get(const mat m, const size_t i, const size_t j);
void mat_set(mat m, const size_t i, const size_t j, const double val);
int mat_set_subm(mat m, const mat n, const size_t k1, const size_t l1, \
				 const size_t k2, const size_t l2);

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
int mat_cp(mat m, const mat n);
int mat_cp_subm(mat m, const mat n, const size_t k1, const size_t l1, \
				const size_t k2, const size_t l2);
int mat_eqadd(mat m, const mat n);
int mat_add(mat m, const mat n, const mat o);
int mat_eqsub(mat m, const mat n);
int mat_sub(mat m, const mat n, const mat o);
#define mat_eqmul_l(m,n) mat_mul(m,n,m);
#define mat_eqmul_r(m,n) mat_mul(m,m,n);
/*!
 @fn void mat_mul(mat m, const mat n, const mat o)
 @brief Store in \b m the matrix product of \b n and \b o.
 @see #mat_eqmul
 
 @warning \b NCOL(n) and \b NROW(o) must be equals, \b m must have the same number of rows than \b n and the same
 number of columns than \b o. Fatal error will be generated if precedent conditions are not satisfied.
 */
int mat_mul(mat m, const mat n, const mat o);
#define mat_eqinv(m) mat_inv(m,m);
int mat_inv(mat m, const mat n);
/*!
 @fn void mat_eqmulp(mat m, const mat n)
 @brief Multiply \b m by \b n "coefficient by coefficient".
 @see #mat_mulp
 
 @warning \b m and \b n must have the same dimensions, fatal error will be generated if it is not the case.
*/
int mat_eqmulp(mat m, const mat n);
int mat_mulp(mat m, const mat n, const mat o);
int mat_eqmuls(mat m, const double s);
int mat_muls(mat m, const mat n, const double s);
int mat_eqdivp(mat m, const mat n);
int mat_divp(mat m, const mat n, const mat o);
#define mat_eqabs(m) mat_abs(m,m)
int mat_abs(mat m, const mat n);
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
int mat_sqrt(mat m, const mat n);
/*!
 @fn void mat_mean(mat mean, const mat *m, const int size)
 @brief Compute the "coefficient by coefficient" mean of an array of matrices.
 @see #mat_print
 
 @param mean matrix to store the mean
 @param m input array of matrices
 @param size number of matrices in the array
 @warning \b mean and all the elements of \b m must have the same dimensions, fatal error will be generated if it is not the case.
 */
int mat_mean(mat mean, const mat *m, const size_t size);
/*!
 @fn void mat_meansig(mat mean, mat sig, const mat *m, const int size)
 @brief Compute the "coefficient by coefficient" mean and standard deviation of an array of matrices.
 
 @param mean matrix to store the mean
 @param sig matrix to store the standard deviation
 @param m input array of matrices
 @param size number of matrices in the array
 @warning \b mean, \b sig and all the elements of \b m must have the same dimensions, fatal error will be generated if it is not the case.
 */
int mat_meansig(mat mean, mat sig, const mat *m, const size_t size);

/** I/O **/
/*!
 @fn void mat_dump(FILE *stream, mat m)
 @brief Print the coefficients of a matrix on a character stream. Columns are separated by a space character and rows 
 are separated by a new line character. Value are printed in scientific notation with 11 digits.
 
 @param stream character stream where the matrix will be printed
 @param m matrix to print
*/
double mat_elsum(mat m);
double mat_elmean(mat m);
void mat_dump(FILE *stream, mat m);
/*!
 @def mat_print(m)
 @brief Print \b m on the standard output.
 @see ::mat_dump
 */
#define mat_print(m) mat_dump(stdout,m)

__END_DECLS

#endif