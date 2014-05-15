#ifndef __PS_TENSOR_
#define __PS_TENSOR_

typedef struct Etens_str
{
   struct Etens_str * next;
   tensor tcoef;
   char eps[4];
} Etens_str;

#define X_MARK -121
typedef struct Etens_str  * Etens;

extern void    delEtens(Etens s);
extern void    multEtensInt(Etens *t , long l);
extern void    multEtensPoly(Etens*, poly);
extern void    multEtensTens(Etens * spn, tensor tns);
extern void    addEtens(Etens * t1, Etens t2);

extern Etens   mult2Etens(Etens t1, Etens t2);
extern Etens   copyEtens(Etens s);
extern Etens   etens1(void);

#endif
