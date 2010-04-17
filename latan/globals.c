#include <latan/globals.h>
#include <latan/includes.h>

const stringbuf latan_name = PACKAGE_NAME;
const stringbuf latan_version = PACKAGE_VERSION;

unsigned int latan_binomial(const unsigned int n, const unsigned int p)
{
	unsigned int* b;
	unsigned int i,j;
	unsigned int res;
	
	MALLOC_NOERRET(b,unsigned int*,n+1);
	
	b[0]=1;
	for (i=1;i<=n;i++)
	{
		b[i] = 1;
		for (j=i-1;j>0;j--)
		{
			b[j] += b[j-1];
		}
	}
	res = b[p];
	
	FREE(b);
	
	return res;
}