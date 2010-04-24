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

/* verbosity flags */
#define QUIET false
#define VERB true
#define VVERB 2

/* error handling */
#include <latan/error.h>

/* size_t type */
#include <sys/types.h>

/* memory allocation */
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
#define FREE(pt)\
{\
	if (pt != NULL)\
	{\
		free(pt);\
		pt = NULL;\
	}\
}
/* file opening */
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
#ifndef STRING_LENGTH
#define STRING_LENGTH 128
#endif

typedef char stringbuf[STRING_LENGTH];

extern const stringbuf latan_name;
extern const stringbuf latan_version;

/* type for random generator state */
#define RLXG_STATE_SIZE 105
typedef int randgen_state[RLXG_STATE_SIZE];

/* math and physics things */
#define SQ(x) ((x)*(x))
#define DRATIO(a,b) (((double)(a))/((double)(b)))
#define C_PI 3.1415926535897932384626433832795028841970
/** 1 fm = 0.005067731 MeV^(-1) **/
#define C_FM_IMEV 0.005067731

unsigned int latan_binomial(const unsigned int n, const unsigned int p);

__END_DECLS

#endif