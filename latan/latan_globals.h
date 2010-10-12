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

__BEGIN_DECLS

/* boolean type */
#ifndef __cplusplus
typedef enum
{
	false = 0,
	true = 1
} bool;
#endif

/* string buffers */
#ifndef STRING_LENGTH
#define STRING_LENGTH 512
#endif

typedef char strbuf[STRING_LENGTH];

__END_DECLS

/* verbosity flags */
#define QUIET 0
#define VERB 1
#define DEBUG 2

/* matrix operation flags */
#define CPU_MAT_OP 0
#define GPU_MAT_OP 1

/* parity flags */
#define EVEN 0
#define ODD 1

/* size_t type */
#include <sys/types.h>

/* error handling */
#include <latan/latan_error.h>

/* matrices */
#ifdef LATAN_ERRNO_DEF
#include <latan/latan_mat.h>
#endif

/* memory allocation */
#define MALLOC(pt,typ,size)\
{\
	pt = (typ)(malloc((size_t)(size)*sizeof(*pt)));\
	if (pt == NULL)\
	{\
		LATAN_ERROR("memory allocation failed",LATAN_ENOMEM);\
	}\
}
#define MALLOC_ERRVAL(pt,typ,size,value)\
{\
	pt = (typ)(malloc((size_t)(size)*sizeof(*pt)));\
	if (pt == NULL)\
	{\
		LATAN_ERROR_VAL("memory allocation failed",LATAN_ENOMEM,value);\
	}\
}
#define MALLOC_NOERRET(pt,typ,size)\
{\
	pt = (typ)(malloc((size_t)(size)*sizeof(*pt)));\
	if (pt == NULL)\
	{\
		LATAN_ERROR_NORET("memory allocation failed",LATAN_ENOMEM);\
	}\
}
#define REALLOC(pt,pt_old,typ,size)\
{\
	pt = (typ)(realloc(pt_old,(size_t)(size)*sizeof(*pt)));\
	if (pt == NULL)\
	{\
		LATAN_ERROR("memory allocation failed",LATAN_ENOMEM);\
	}\
}
#define REALLOC_ERRVAL(pt,pt_old,typ,size,value)\
{\
	pt = (typ)(realloc(pt_old,(size_t)(size)*sizeof(*pt)));\
	if (pt == NULL)\
	{\
		LATAN_ERROR_VAL("memory allocation failed",LATAN_ENOMEM,value);\
	}\
}
#define REALLOC_NOERRET(pt,pt_old,typ,size)\
{\
	pt = (typ)(realloc(pt_old,(size_t)(size)*sizeof(*pt)));\
	if (pt == NULL)\
	{\
		LATAN_ERROR_NORET("memory allocation failed",LATAN_ENOMEM);\
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
		strbuf _errmsg;\
		sprintf(_errmsg,"error opening file %s",fname);\
		LATAN_ERROR(_errmsg,LATAN_EFAULT);\
	}\
}
#define FOPEN_ERRVAL(f,fname,mode,value)\
{\
	f = fopen(fname,mode);\
	if (f == NULL)\
	{\
		strbuf _errmsg;\
		sprintf(_errmsg,"error opening file %s",fname);\
		LATAN_ERROR_VAL(_errmsg,LATAN_EFAULT,value);\
	}\
}
#define FOPEN_NOERRET(f,fname,mode)\
{\
	f = fopen(fname,mode);\
	if (f == NULL)\
	{\
		strbuf _errmsg;\
		sprintf(_errmsg,"error opening file %s",fname);\
		LATAN_ERROR_NORET(_errmsg,LATAN_EFAULT);\
	}\
}

__BEGIN_DECLS

/* LatAnalyze environment access */
void latan_get_name(strbuf name);
void latan_get_version(strbuf version);
int latan_get_verb(void);
int latan_get_mat_op(void);
#ifdef LATAN_ERRNO_DEF
latan_errno latan_set_verb(int verb);
latan_errno latan_set_mat_op(int mat_op);
#endif

/* LatAnalyze message function */
void latan_printf(const int verb, const strbuf fmt, ...);

__END_DECLS

#endif