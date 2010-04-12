#ifndef LATAN_ERROR_H_
#define LATAN_ERROR_H_

#include <latan/globals.h>

__BEGIN_DECLS

enum
{
	LATAN_FAILURE	= -1,	/* generic failure statement			*/
	LATAN_SUCCESS	= 0,	/* all is going well !					*/
	LATAN_EDOM		= 1,	/* input domain error, e.g sqrt(-1)		*/
	LATAN_EFAULT	= 3,	/* invalid pointer						*/
	LATAN_EINVAL	= 4,    /* invalid argument supplied by user	*/
	LATAN_ENOMEM	= 8,	/* malloc error							*/
	LATAN_EBADLEN	= 19,	/* matrix dimension error				*/
	LATAN_ENOTSQR	= 20,	/* matrix is not square error			*/
	LATAN_ELATSYN	= 33,	/* syntax error reading input file		*/
	LATAN_ESYSTEM	= 34	/* system error							*/
};

typedef void latan_error_handler_t(const char*, const char*, int, int);

void latan_error(const char* reason, const char *file, int line, \
				 int latan_errno);

latan_error_handler_t*\
latan_set_error_handler(latan_error_handler_t* new_handler);

latan_error_handler_t* latan_set_error_handler_off(void);

#define LATAN_ERROR(reason, latan_errno)\
{\
	latan_error(reason,__FILE__,__LINE__,latan_errno);\
	return latan_errno;\
}

#define LATAN_ERROR_VAL(reason, latan_errno, value)\
{\
	latan_error(reason,__FILE__,__LINE__,latan_errno);\
	return value;\
}

#define LATAN_ERROR_VOID(reason, latan_errno)\
{\
	latan_error(reason,__FILE__,__LINE__,latan_errno);\
	return;\
}

#define LATAN_ERROR_NULL(reason, latan_errno)\
LATAN_ERROR_VAL(reason,latan_errno,NULL)

#define LATAN_ERROR_NORET(reason, latan_errno)\
latan_error(reason,__FILE__,__LINE__,errno)

#define LATAN_ERROR_SELECT_2(a,b)       \
((a) != LATAN_SUCCESS ? (a) : ((b) != LATAN_SUCCESS ? (b) : LATAN_SUCCESS))
#define LATAN_ERROR_SELECT_3(a,b,c)     \
((a) != LATAN_SUCCESS ? (a) : LATAN_ERROR_SELECT_2(b,c))
#define LATAN_ERROR_SELECT_4(a,b,c,d)   \
((a) != LATAN_SUCCESS ? (a) : LATAN_ERROR_SELECT_3(b,c,d))
#define LATAN_ERROR_SELECT_5(a,b,c,d,e) \
((a) != LATAN_SUCCESS ? (a) : LATAN_ERROR_SELECT_4(b,c,d,e))
#define LATAN_UPDATE_STATUS(status,instruction)\
{\
	int cstat;\
	cstat=instruction;\
	status = LATAN_ERROR_SELECT_2(status,cstat);\
}
__END_DECLS

#endif