
#ifndef __LHA__ 
#define __LHA__

extern double  alphaspdf_(double *);

extern void   getdatapath_(char* dirpath, int len);
extern void   initpdfsetbynamem_(int *nSet,char *name, int len);
extern void   numberpdfm_(int* nSet,int *);
extern void   evolvepdfm_(int* nSet,double *x,double *Q,double *f);
extern void   initpdfm_(int* nSet,int *);

#endif
