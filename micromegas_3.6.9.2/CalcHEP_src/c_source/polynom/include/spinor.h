#ifndef __SPINOR_
#define __SPINOR_

#include "ps_tensor.h"

typedef struct SpinTensor_str
{
   struct SpinTensor_str * next;
   tensor tcoef;
   char g5;
   char l;
   char g[4];
} SpinTensor_str;

typedef struct SpinTensor_str  * SpinTensor;

extern int     spinLength;



extern void    delSpin(SpinTensor s);
extern void    multSpinInt(SpinTensor *t , long l);
extern void    multSpinPoly(SpinTensor*, poly);
extern void    multSpinTens(SpinTensor * spn, tensor tns);
extern void    addSpin(SpinTensor * t1, SpinTensor t2);

extern SpinTensor mult2Spin(SpinTensor t1, SpinTensor t2, int forspur);
extern SpinTensor copySpin(SpinTensor s);
extern SpinTensor spin1(void); 

extern Etens calcspur(SpinTensor );
#endif
