
#ifndef __TEST_WRT_
#define __TEST_WRT_
#include"spinor.h"
#include"symb_reader.h"

extern void  writepoly(poly  p);
extern void  writetens(tensor  p);
extern void  writespinor(SpinTensor p);
extern void  writeEtens(Etens p);
extern void  symb_print(char*txt,symb_data m);

#endif
