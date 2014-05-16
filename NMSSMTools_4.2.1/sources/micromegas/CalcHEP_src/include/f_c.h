

#ifndef F2C_INCLUDE
#define F2C_INCLUDE

#ifndef ABS
#define ABS(x) ((x) >= 0 ? (x) : -(x)) 
#endif

#ifndef MAX
#define MAX(a,b)        (((a) > (b)) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b)        (((a) < (b)) ? (a) : (b))
#endif

extern double pow_dl(double ap, long bp);

#endif
