#include <latan/latan_includes.h>
#include <latan/latan_error.h>

static void no_error_handler(const char* reason, const char* file, int line,\
							 int no);

latan_error_handler_t* latan_error_handler = NULL;

void latan_error(const char* reason, const char* file, int line,\
				 int no)
{
	stringbuf name,version;
	
	if (latan_error_handler) 
    {
		(*latan_error_handler)(reason,file,line,no);
		return;
    }
	
	latan_get_name(name);
	latan_get_version(version);
	fprintf(stderr,"%s v%s error %d: %s (%s:%d)\n",name,version,\
			no,reason,file,line);
	fflush(stderr);
	abort();
}

void latan_warning(const char* reason, const char* file, int line,\
				   int no)
{
	stringbuf name,version;

	latan_get_name(name);
	latan_get_version(version);
	fprintf(stderr,"%s v%s warning %d: %s (%s:%d)\n",name,version,\
			no,reason,file,line);
	fflush(stderr);
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
							 int no)
{
	/* do nothing */
	reason = NULL;
	file = NULL;
	line = 0;
	no = LATAN_SUCCESS;
}
