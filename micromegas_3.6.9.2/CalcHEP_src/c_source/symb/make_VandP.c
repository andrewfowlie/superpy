#include<stdio.h>
#include <string.h>
#include <math.h>
#include "model.h"
#include "read_mdl.h"
#include "procvar.h"
#include "reader_c.h"
#include "parser.h"

#include "nType.h"

int menulevel; // for compilation 

int main(int argv, char**argc)
{
  int i,i10,nv,nLn,L;
  FILE*f,*fExt;  
  int nVar=0,nFunc=0,first;
  int mode; 
  char path[200];
  char * CalcHEP=NULL;
  
  if(argv<3) { printf("Arguments expected: 1)path to model files; 2) model number.\n"); return 1;}

  L=strlen(argc[0]);
  CalcHEP=malloc(L+10);
  strcpy(CalcHEP,argc[0]);
  CalcHEP[L-15]=0;
  if(sscanf(argc[2],"%d",&L)!=1) { printf("Second argument should be a number\n"); return 1;}
   if(argv>=4)sscanf(argc[3],"%d",&mode); else mode=0;
   makeVandP(1,argc[1],L,mode,CalcHEP);  
   return 0;
}

