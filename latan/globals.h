/*!
 @file macros.h
 @brief Miscellaneous utility macros.
 @author <A HREF=mailto:antonin.portelli@gmail.com>Antonin Portelli</A>
*/
#ifndef LATAN_GLOBALS_H_
#define LATAN_GLOBALS_H_

/* macros for including from C++ code */
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS
#define __END_DECLS
#endif

/* portability of exit status */
#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif

/* verbosity flags */
#define QUIET false
#define VERB true
#define VVERB 2

/* error handling */
#include <latan/error.h>

/* memory allocation */
/*!
 @def MALLOC(pt,typ,size)
 @brief Allocate memory for an array, generates a fatal error in case of failure.
 
 @param pt pointer on the first element of the array
 @param typ pointer type of the elements of the array
 @param size number of elements of the array
 */
#define MALLOC(pt,typ,size)\
{\
	pt = (typ)(malloc((size_t)(size)*sizeof(*pt)));\
	if (pt == NULL)\
	{\
		LATAN_ERROR("memory allocation failed",LATAN_EFAULT);\
	}\
}
#define MALLOC_ERRVAL(pt,typ,size,value)\
{\
	pt = (typ)(malloc((size_t)(size)*sizeof(*pt)));\
	if (pt == NULL)\
	{\
		LATAN_ERROR_VAL("memory allocation failed",LATAN_EFAULT,value);\
	}\
}
#define MALLOC_NOERRET(pt,typ,size)\
{\
	pt = (typ)(malloc((size_t)(size)*sizeof(*pt)));\
	if (pt == NULL)\
	{\
		LATAN_ERROR_NORET("memory allocation failed",LATAN_EFAULT);\
	}\
}
/*!
 @def REALLOC(pt,pt_old,typ,size)
 @brief Reallocate memory for an array, backup of the old values is made. Generates a fatal error in case of failure.
 
 @param pt pointer on the first element of the array
 @param pt_old pointer on the first element of the array to backup
 @param typ pointer type of the elements of the array
 @param size number of elements of the array
 */
#define REALLOC(pt,pt_old,typ,size)\
{\
	pt = (typ)(realloc(pt_old,(size_t)(size)*sizeof(*pt)));\
	if (pt == NULL)\
	{\
		LATAN_ERROR("memory allocation failed",LATAN_EFAULT);\
	}\
}
#define REALLOC_ERRVAL(pt,pt_old,typ,size,value)\
{\
	pt = (typ)(realloc(pt_old,(size_t)(size)*sizeof(*pt)));\
	if (pt == NULL)\
	{\
		LATAN_ERROR_VAL("memory allocation failed",LATAN_EFAULT,value);\
	}\
}
#define REALLOC_NOERRET(pt,pt_old,typ,size)\
{\
	pt = (typ)(realloc(pt_old,(size_t)(size)*sizeof(*pt)));\
	if (pt == NULL)\
	{\
		LATAN_ERROR_NORET("memory allocation failed",LATAN_EFAULT);\
	}\
}
/*!
 @def FREE(pt)
 @brief Desallocate a previoulsy memory allocated array.
 
 @param pt pointer on the first element of the array
 */
#define FREE(pt)\
{\
	if (pt != NULL)\
	{\
		free(pt);\
		pt = NULL;\
	}\
}
/* file opening */
/*!
 @def FOPEN(f,fname,mode)
 @brief Macro opening a file, generates a fatal error in case of failure.
 
 @param f pointer on a \b FILE structure (standard C) to access the file
 @param fname string for the file name
 @param mode string for the opening mode of the file (\e cf. \b fopen standard C function)
 */
#define FOPEN(f,fname,mode)\
{\
	f = fopen(fname,mode);\
	if (f == NULL)\
	{\
		stringbuf _errmsg;\
		sprintf(_errmsg,"error opening file %s",fname);\
		LATAN_ERROR(_errmsg,LATAN_EFAULT);\
	}\
}
#define FOPEN_ERRVAL(f,fname,mode,value)\
{\
	f = fopen(fname,mode);\
	if (f == NULL)\
	{\
		stringbuf _errmsg;\
		sprintf(_errmsg,"error opening file %s",fname);\
		LATAN_ERROR_VAL(_errmsg,LATAN_EFAULT,value);\
	}\
}
#define FOPEN_NOERRET(f,fname,mode)\
{\
	f = fopen(fname,mode);\
	if (f == NULL)\
	{\
		stringbuf _errmsg;\
		sprintf(_errmsg,"error opening file %s",fname);\
		LATAN_ERROR_NORET(_errmsg,LATAN_EFAULT);\
	}\
}

__BEGIN_DECLS

/* boolean type */
typedef enum
{
	false = 0,
	true = 1
} bool;

/* string buffers */
/*!
 @def STRING_LENGTH
 @brief Number of characters of a string defined by the ::stringbuf type.
*/
#define STRING_LENGTH 128
/*!
 @typedef stringbuf
 @brief String type of fixed length #STRING_LENGTH
*/
typedef char stringbuf[STRING_LENGTH];

extern const stringbuf latan_name;
extern const stringbuf latan_version;

/* math and physics things */
/*!
 @def SQ(x)
 @brief Macro for square operation.
 */
#define SQ(x) ((x)*(x))
#define DRATIO(a,b) (((double)(a))/((double)(b)))
/*!
 @def C_PI
 @brief \f$\pi\f$ math constant.
 */
#define C_PI 3.1415926535897932384626433832795028841970
/** 1 fm = 0.005067731 MeV^(-1) **/
#define C_FM_IMEV 0.005067731

unsigned int latan_binomial(const unsigned int n, const unsigned int p);

__END_DECLS

#endif