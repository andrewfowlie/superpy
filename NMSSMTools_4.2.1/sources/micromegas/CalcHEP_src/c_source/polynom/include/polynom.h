#ifndef __POLYNOM_
#define __POLYNOM_

#include"lnum.h"

#define VAR_NAME_SIZE 12

#define poly struct monom *
typedef struct monom
   {  
      poly       next;
      NUM_TYPE   num;
      unsigned long  power[2];
   } monom;
#undef poly
typedef struct monom *poly;


typedef struct tensor_str
   {  
      struct tensor_str * next;
      poly re;
      poly im;
      char tens[2*sizeof(long)];
   } tensor_str;

typedef struct tensor_str * tensor;

typedef struct  varinfo
{
   char  name[VAR_NAME_SIZE+6];  
   unsigned long maxdeg,  zerodeg;
   int  wordpos;
   int  num;
}varinfo;
                                    



extern poly    garbage;
extern poly   *contracts;
extern int    monomLength,maxLength;



extern void  delpoly(poly  * p);

extern poly  plusone(void);
extern poly  copypoly(poly  p);
extern void  sewpoly(poly  * p1,	poly  * p2);
extern void  multpolyint(poly     * p, long    i);
extern poly  multtwopoly(poly  q1, poly  q2);
extern poly  scalarmult(int  p1,	int  p2);
extern void  assignsclmult(int  p1,	int  p2,	poly      p);
extern void  deltensor(tensor  * t);
extern tensor  copytens(tensor  t);
extern void  sewtens(tensor * t1,tensor *t2);
extern void  multtensint(tensor    * t,	long i);
extern void  multtensComplexpoly(tensor* t,poly re,poly im);
extern void  multtenspoly(tensor  * t,  poly  p);
extern void (*memoryInfo) (int);
extern void tensRealPart(tensor *t);
extern void tensImPart(tensor *t);
extern void  symb_start(int nvar, varinfo * Vars, int nSpin, int nIndex, int nMom);
extern void  makeNewGarbage(void);

/*
#define delunit(p) {((poly)(p))->next=garbage; garbage=(poly)(p); (p)=NULL;}
*/

#define NewUnit(p) {if(!garbage) makeNewGarbage(); p=(void*)garbage; garbage=garbage->next;}

extern void  delunit(void* p);

#endif
