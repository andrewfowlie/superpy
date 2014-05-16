/*
Copyright (C) 2006, Alexander Pukhov
*/

#include <ctype.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "parser.h"
#include "calc.h"
#include "simpson.h"

int err_code=0;
int main(int argn, char** argc)
{ 
  double a,b,r;
 
  if(argn!=4) 
  { printf("this programs needs 3 argumnents:\n"
           "1. integrand\n"
           "2. down limit\n"
           "3. up limit\n");
    exit(100);       
  }           
  userFuncTxt=argc[2];
  a=userFuncNumC();
  if(rderrcode) 
  { printf( "Can not calculate down limit\n");
    return rderrcode;
  }  
  userFuncTxt=argc[3];
  b=userFuncNumC();
  if(rderrcode) 
  { printf( "Can not calculate up limit\n");
    return rderrcode;
  }  
  
  userFuncTxt=argc[1];
  
  r=gauss345(userFuncNumI,a,b,1.E-6,NULL);
  if(rderrcode) { printf("Can not calculate function\n");}
  else  printf("%.10G\n",r);
  return rderrcode;
}          

