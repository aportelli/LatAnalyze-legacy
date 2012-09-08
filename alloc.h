#ifndef ALLOC_H_
#define ALLOC_H_

extern void *malloc(size_t n);
extern void *realloc(void *pt, size_t n);
void *rpl_malloc(size_t n);
void *rpl_realloc (void *pt, size_t n);

#endif
