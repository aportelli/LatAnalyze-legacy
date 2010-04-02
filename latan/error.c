#include <latan/includes.h>
#include <latan/error.h>

static void no_error_handler(const char* reason, const char* file, int line,\
							 int latan_errno);

latan_error_handler_t* latan_error_handler = NULL;

void latan_error(const char* reason, const char* file, int line,\
				 int latan_errno)
{
	if (latan_error_handler) 
    {
		(*latan_error_handler)(reason,file,line,latan_errno);
		return;
    }
	
	fprintf(stderr,"%s v%s error %d: %s (%s:%d)\n",latan_name,latan_version,\
			latan_errno,reason,file,line);
	fflush(stderr);
	abort();
}

latan_error_handler_t*\
latan_set_error_handler(latan_error_handler_t* new_handler)
{
	latan_error_handler_t* previous_handler = latan_error_handler;
	latan_error_handler = new_handler;
	return previous_handler;
}


latan_error_handler_t* latan_set_error_handler_off (void)
{
	latan_error_handler_t* previous_handler = latan_error_handler;
	latan_error_handler = no_error_handler;
	return previous_handler;
}

static void no_error_handler(const char* reason, const char* file, int line,\
							 int latan_errno)
{
	/* do nothing */
	reason = NULL;
	file = NULL;
	line = 0;
	latan_errno = 0;
}
