#include <latan/latan_globals.h>
#include <latan/latan_includes.h>

typedef struct
{
	stringbuf name;
	stringbuf version;
	int verb;
	int minimize_lib;
	stringbuf prop_mark;
	stringbuf prop_idfmt;
} latan_env;

static latan_env env = 
{
	PACKAGE_NAME,		\
	PACKAGE_VERSION,	\
	QUIET,				\
	GSL,				\
	"PROP",				\
	"%s_%s_%d_%d"
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

int latan_get_minimize_lib(void)
{
	return env.minimize_lib;
}

latan_errno latan_set_minimize_lib(int minimize_lib)
{
	if ((minimize_lib < 0)||(minimize_lib > 1))
	{
		LATAN_ERROR("minimize library flag invalid",LATAN_EINVAL);
	}
	
	env.minimize_lib = minimize_lib;
	
	return LATAN_SUCCESS;
}

void latan_get_prop_mark(stringbuf prop_mark)
{
	strcpy(prop_mark,env.prop_mark);
}

void latan_set_prop_mark(const stringbuf prop_mark)
{
	strcpy(env.prop_mark,prop_mark);
}

void latan_get_prop_idfmt(stringbuf prop_idfmt)
{
	strcpy(prop_idfmt,env.prop_idfmt);
}

void latan_set_prop_idfmt(const stringbuf prop_idfmt)
{
	strcpy(env.prop_idfmt,prop_idfmt);
}