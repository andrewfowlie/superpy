#ifndef __READER0_
#define __READER0_

#include "sets.h"

char*parseVertex(int v, int forReduce);
char*fermPropagTxt(int v, int l, int forReduce);
char*tPropagator6mass4(int v,int l,set*indexs);
char * spin3_2_propagator(int v,int l,int forReduce);
#endif
