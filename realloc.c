#include <config.h>
#undef realloc
#include <sys/types.h>

void *realloc(void *pt, size_t n);

/* Reallocate an N-byte block of memory from the heap. If N is zero, allocate a 1-byte block. */
void *rpl_realloc (void *pt, size_t n) 
{
	if (n == 0)
	{
		n = 1;
	}
	if (pt == NULL)
	{
		return malloc(n);
	}
	
	return realloc(pt,n);
}