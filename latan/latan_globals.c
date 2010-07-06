#include <latan/latan_globals.h>
#include <latan/latan_includes.h>

typedef struct
{
	stringbuf name;
	stringbuf version;
	int verb;
} latan_env;

static latan_env env = 
{
	PACKAGE_NAME,		\
	PACKAGE_VERSION,	\
	QUIET,				\
};

void latan_get_name(stringbuf name)
{
	strcpy(name,env.name);
}

void latan_get_version(stringbuf version)
{
	strcpy(version,env.version);
}

int latan_get_verb(void)
{
	return env.verb;
}

latan_errno latan_set_verb(int verb)
{
	if ((verb < 0)||(verb > 2))
	{
		LATAN_ERROR("verbosity level invalid",LATAN_EINVAL);
	}
	
	env.verb = verb;
	
	return LATAN_SUCCESS;
}

void latan_printf(const int verb, const stringbuf fmt, ...)
{
	va_list args;
	stringbuf head,tail,name,debug,version;
	
	if ((latan_get_verb() >= verb)&&(verb >= VERB))
	{
		latan_get_name(name);
		latan_get_version(version);
		if (verb == DEBUG)
		{
			strcpy(debug," - DEBUG");
		}
		else
		{
			strcpy(debug,"");
		}
		
		sprintf(head,"[%s v%s%s]",name,version,debug);
		va_start(args,fmt);
		vsprintf(tail,fmt,args);
		va_end(args);
		printf("%s %s",head,tail);
	}
}