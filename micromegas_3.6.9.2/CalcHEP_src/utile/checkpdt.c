#include <stdio.h>
#include "pdt.h"

int main(int np, char ** par)
{
  pdtStr S1;
  double x=0.2,q=1.E3,f;
  int ip=0;
  int err;
  long pNum;
  if(np!=5)
  {  printf("This function needs 4 parameters:\n"
                   "1 name of pdt-file\n"
                   "2 parton code\n"
                   "3 x parameter\n"
                   "4 q parameter\n");
     exit(0);
  }
  pdtList *list=NULL;
  sscanf(par[2],"%ld",&pNum);
  sscanf(par[3],"%lf",&x);
  sscanf(par[4],"%lf",&q); 

  makePdtList(par[1], pNum, &list);
  
  for(ip=0;list;list=list->next)
  { 
   ip=list->position;
   if(getPdtData(par[1], list->position,&S1)) continue;
   printf("%s   f=%E\n", list->name,interFunc(x,q, &S1));
   if(!list->next && S1.alpha)  printf("alpha=%.3E\n", interAlpha(q, &S1)); 
   freePdtData(&S1);
  }

  delPdtList(list);  
}
