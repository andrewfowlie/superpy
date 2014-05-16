#ifndef __OPTIMISE_
#define __OPTIMISE_

#include"syst.h"
#include"lnum.h"

#define infoptr struct inforec *
#define varptr struct var_rec *
typedef struct var_rec
   {
      varptr      next;
      int         sgn;
      infoptr     coef;
      int       vars[STRSIZ];
   }  var_rec;
#undef varptr
#undef infoptr
typedef struct var_rec *varptr;


#define infoptr struct inforec *
typedef struct inforec
   {
      infoptr    next;
      unsigned   count;
      char       name[20];
      enum{numb,expr,rnumb} consttype;
      varptr     const_;
      NUM_TYPE   ival;
      double     rval;
   }  inforec;
#undef infoptr
typedef struct inforec *infoptr;


typedef void * ( * smplemit)(varptr ex);
typedef void * ( * vfact   )(int ch, int deg,  void * pmult, void * psum);
typedef void * ( * cfact   )(infoptr c, void * pmult, void * psum);
	
extern void  initinfo(void);
extern void  readpolynom(varptr  * expr_);
extern void * emitexpr(varptr    ex,smplemit  smplemitfun,
	vfact     vfactfun,  cfact     cfactfun);
extern int equalexpr(varptr  v1,	varptr  v2);

#define minvarrec (sizeof(struct var_rec) - sizeof(int)*(STRSIZ-1))

extern infoptr info;

extern int firstVar;

#endif
