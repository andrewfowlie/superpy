/*
Copyright (C) 1997, Alexander Pukhov
*/

#include <ctype.h>
#include <math.h>
#include <string.h>
#include "syst.h"
#include "parser.h"
#include "getmem.h"
#include "f_c.h"
#include "read_func.h"


static int (*nameToVal)(char *,double *);
static int isAble;

static void*  rd_num(char* s)
{
   double    * p;
   int rc;

   p = (double *) getmem_(sizeof(double)); 
   
   if(!nameToVal){ printf("Error in programming,nameToVal==NULL\n"); sortie(90);}
   
   rc=(*nameToVal)(s,p);
   if(!rc)      
   { 
      if(isdigit(*s))  rderrcode=typemismatch;else rderrcode=unknownidentifier; 
      return NULL; 
   } 
   else if(rc==-1) isAble=0; 
   return (void*)p;     
}


static void*  act_num(char* ch,int n, void**args)
{  double  p1, p2,p3;
   if(!isAble) { *(double*)args[0] = 0.; return args[0];}
   p1= *(double*)args[0];
   if(n>=2) p2=*(double*)args[1];
   if(n>=3) p3=*(double*)args[2];
   
   switch (ch[0])
   { 
     case '-': if(n==1) p1*=-1;else p1-=p2; break;
     case '+': p1+=p2;  break;
     case '*': p1*=p2;  break;
     case '/': if(p2==0.)rderrcode=naninoperation; else p1/=p2;  break;
     case '^': p1=pow(p1,p2);
               break;
     case '.': rderrcode=typemismatch; 
               break;
     default : switch(n) 
               {      
                   case 1: 
                      if(!strcmp(ch,"sqrt")) p1=sqrt(p1);
                 else if(!strcmp(ch,"sin"))  p1=sin(p1);
                 else if(!strcmp(ch,"cos"))  p1=cos(p1);
                 else if(!strcmp(ch,"tan"))  p1=tan(p1);
                 else if(!strcmp(ch,"asin")) p1=asin(p1);
                 else if(!strcmp(ch,"acos")) p1=acos(p1);
                 else if(!strcmp(ch,"atan")) p1=atan(p1);
                 else if(!strcmp(ch,"exp"))  p1=exp(p1);
                 else if(!strcmp(ch,"log"))  p1=log(p1);
                 else if(!strcmp(ch,"fabs")) p1=fabs(p1);
                 else isAble =0;
                   break;
                   case 2:
                      if(!strcmp(ch,"atan2")) p1=atan2(p1,p2); 
                 else isAble=0; 
                   break;
                   case 3:
                      if(!strcmp(ch,"if")) {if(p1>0) p1=p2; else p1=p3;}
                      else isAble=0; 
                   break;   
                 default :isAble=0;
               }
   } 
   if(rderrcode) return NULL;
   if(!isAble)
   {  int i;
      if(strcmp(ch,"min")==0)
      { for(i=1; i<n; i++) if(p1> *(double*)args[i]) p1=*(double*)args[i];
        isAble=1;
      }else if(strcmp(ch,"max")==0)
      { for(i=1; i<n; i++) if(p1< *(double*)args[i]) p1=*(double*)args[i];
        isAble=1;
      }else p1=0;
   }

   *(double*)args[0] = p1; return args[0];
}

int calcExpression(char *s,int(*nameToDouble)(char *,double *), double *p)
{
  marktp  heapbeg;
  double * r;
  isAble=1;

  nameToVal=nameToDouble;

  mark_(&heapbeg);
  r = (double *)readExpression(s,rd_num,act_num,NULL);
  if(rderrcode==0 && !isAble) rderrcode=unknownfunction; 
 
  if(!rderrcode) *p=*r;
  if(!rderrcode && !isfinite(*r)) rderrcode=cannotevaluate;
 
  release_(&heapbeg); 
  return rderrcode;
}          
