#ifndef _PMALLOC_H_
#define _PMALLOC_H_ 1

#define malloc(x)      pmalloc((x), __FILE__, __LINE__)
#define calloc(n, x)   pcalloc((n), (x), __FILE__, __LINE__)
#define realloc(p, n)  prealloc((p), (n), __FILE__, __LINE__)
#define free(p)        pfree((p), __FILE__, __LINE__)

void *pmalloc(size_t size, char *f, int nline);
void *pcalloc(size_t n, size_t size, char *f, int nline);
void *prealloc(void *p, size_t size, char *f, int nline);
void pfree(void *p, char *f, int nline);

void pmalloc_check(void);

#endif
