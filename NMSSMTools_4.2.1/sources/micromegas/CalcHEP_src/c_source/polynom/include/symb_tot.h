#ifndef  __SYMB_TOT_
#define __SYMB_TOT_

#include"polynom.h"
#include"spinor.h"
#include"tensor.h"
#include"ps_tensor.h"


#define errortp  0
#define numbertp 1
#define polytp   2
#define rationtp 3
#define vectortp 4
#define indextp  5
#define tenstp   6
#define spintp   7
#define etenstp  8


typedef struct symb_data
{ union
  { poly p;
    tensor t;
    SpinTensor s;
    Etens et;
  } expr;
  int type;
} symb_data;

                  
extern void symb_clean(symb_data S);
extern symb_data symb_copy(symb_data S);
extern symb_data symb_imult(symb_data S,int del,int factor);
extern symb_data symb_mult(symb_data S1, int del1, symb_data S2, int del2);
extern symb_data symb_sum(symb_data S1, int del1, symb_data S2, int del2);
extern symb_data symb_typeUp(symb_data S, int del,int type);
extern symb_data symb_spur(symb_data S, int del);
extern int symb_testindex(symb_data S);
extern int symb_iszero(symb_data S);

extern symb_data symb_real(symb_data S, int del);
extern symb_data symb_delEps(symb_data S, int del);
#endif
