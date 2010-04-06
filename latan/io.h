/*!
 @file io.h
 @brief I/O functions.
 
 Matrices in data files are expected to be formated in the following way :<BR><BR>
 <TT>
 START_mark name<BR>
 x x ... x<BR>
 .&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;.<BR>
 .&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;.<BR>
 .&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;.<BR>
 x x ... x<BR>
 END_mark
 </TT><BR><BR>
 where \c mark and \c name are strings of your choice to identify the matrix and
 the \c x are the matrix coefficients. Columns are separated by space
 characters and lines by Unix new line characters.<BR>
 For example, propagator data could be writed like :<BR><BR>
 <TT>
 START_PROP PP_00_2_0<BR>
 4.884724e+00<BR>
 3.962649e-01<BR>
 5.028446e-02<BR>
 5.765858e-03<BR>
 6.111912e-04<BR>
 7.593253e-05<BR>
 8.888678e-06<BR>
 1.023063e-06<BR>
 2.396448e-07<BR>
 1.035865e-06<BR>
 8.654427e-06<BR>
 7.408905e-05<BR>
 6.450945e-04<BR>
 5.609604e-03<BR>
 4.925793e-02<BR>
 4.257967e-01<BR>
 END_PROP
 </TT>
 
 @author <A HREF=mailto:antonin.portelli@gmail.com>Antonin Portelli</A> 
*/

/*!
 @example ex_getmat.c Reading mat variables from files.
*/

#ifndef LATAN_IO_H_
#define LATAN_IO_H_

#include <latan/globals.h>
#include <latan/mat.h>
#include <latan/rand.h>

__BEGIN_DECLS

/* general I/O */
/*!
 @fn int get_nfile(const stringbuf manifestfname)
 @brief Get the number of files referenced in a manifest file.
 
 @param manifestfname name of the manifest file
 
 @return number of file in the manifest file.
 @remark in reality this function only count the number of non blank lines in 
 the file \c manifestfname.
 */
int get_nfile(const stringbuf manifestfname);
/*!
 @fn void get_firstfname(stringbuf fname, const stringbuf manifestfname);
 @brief Get the first file name referenced in a manifest file.
 
 @param fname string to store the result
 @param manifestfname name of the manifest file
 */
int get_firstfname(stringbuf fname, const stringbuf manifestfname);

/* mat I/O */
/*!
 @fn void mat_dump(FILE *stream, mat m)
 @brief Print the coefficients of a matrix on a character stream. Columns are separated by a space character and rows 
 are separated by a new line character. Value are printed in scientific notation with 11 digits.
 
 @param stream character stream where the matrix will be printed
 @param m matrix to print
 */
void mat_dump(FILE* stream, mat m);
/*!
 @def mat_print(m)
 @brief Print \b m on the standard output.
 @see ::mat_dump
 */
#define mat_print(m) mat_dump(stdout,m)
/*!
 @fn int mat_load_nrow(const stringbuf mark, const stringbuf matid, const stringbuf inputfname)
 @brief Get the number of rows of a matrix stored in a file.
 
 @param mark mark of the matrix
 @param matid name of the matrix
 @param inputfname name of the file where the matrix is stored
 
 @return number of rows of the designed matrix.
 */
int mat_load_nrow(const stringbuf mark, const stringbuf matid,\
				  const stringbuf inputfname);
/*!
 @fn void mat_get(mat m, const stringbuf mark, const stringbuf matid, 
 const stringbuf inputfname)
 @brief Get a matrix stored in a file.
 
 @param m matrix to store read one
 @param mark mark of the matrix
 @param matid name of the matrix
 @param inputfname name of the file where the matrix is stored
*/
int mat_load(mat m, const stringbuf mark, const stringbuf matid,\
			 const stringbuf inputfname);
/*!
 @fn void mat_get_ar(mat* m, const stringbuf mark, const stringbuf matid,
 const stringbuf manifestfname)
 @brief Get an array of matrices from files. All matrices got the same mark an 
 name and are stored one per file.
 
 @param m array of matrice to store read ones
 @param mark mark of the matrices
 @param matid name of the matrices
 @param manifestfname name of the manifest file referencing files where matrices
 are stored
 */
int mat_load_ar(mat* m, const stringbuf mark, const stringbuf matid,\
				const stringbuf manifestfname);
int mat_save_plotdat(const mat m, const double xstart, const double xstep,\
					 const stringbuf fname);
int mat_save_plotdaterr(const mat dat, const mat sig, const double xstart,\
						const double xstep, const stringbuf fname);

/* random generator state I/O */
int rand_save_gen_state(const stringbuf prefname, const rand_gen_state state);
int rand_load_gen_state(rand_gen_state state, const stringbuf prefname);

__END_DECLS

#endif