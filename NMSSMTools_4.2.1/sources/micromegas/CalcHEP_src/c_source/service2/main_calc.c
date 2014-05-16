/*
Copyright (C) 2002, Alexander Pukhov
*/

#include <ctype.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "parser.h"
#include "calc.h"


int main(int argn, char** argc)
{ char s[500];
  double r;
  int i;
  for(i=1, s[0]=0; i < argn;i++) strcat(s,argc[i]);
  userFuncTxt=s;
  r = userFuncNumC();
  if(!rderrcode)  printf("%.10G\n",r);
  return rderrcode;
}          

