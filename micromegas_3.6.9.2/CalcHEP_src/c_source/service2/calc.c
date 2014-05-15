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

static void*  rd_num(char* s)
{
   double    * p;
   int rc;
   p = (double *) malloc(sizeof(double)); 
   if(strcmp(s,"PI")==0) *p=M_PI;
   else       
   { rc= sscanf(s,"%lf",p);
     if(rc!=1)   {rderrcode=unknownidentifier;}
   }
   return (void*)p;     
}


static void*  act_num(char* ch,int n, void**args)
{  double  p1, p2, p3;
   p1= *(double*)args[0];
   if(n>=2) { p2=*(double*)args[1]; free(args[1]);}
   if(n>=3) { p3=*(double*)args[2]; free(args[2]);}
   
   switch (ch[0])
   { 
     case '-': if(n==1) p1*=-1;else p1-=p2; break;
     case '+': p1+=p2;  break;
     case '*': p1*=p2;  break;
     case '/': if(p2==0.)rderrcode=naninoperation; else p1/=p2;  break;
     case '^':  p1=pow(p1,p2);
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
                 else if(!strcmp(ch,"floor"))p1=floor(p1);
                 else if(!strcmp(ch,"log10"))p1=log10(p1);
                 else rderrcode=unknownfunction;
                   break;
                   case 2:
                      if(!strcmp(ch,"atan2")) p1=atan2(p1,p2); 
                 else if(!strcmp(ch,"min"))   p1= p1>p2 ? p2:p1;
                 else if(!strcmp(ch,"max"))   p1= p1<p2 ? p2:p1; 
                 else rderrcode=unknownfunction; ; 
                   break;
                   case 3:
                   if(!strcmp(ch,"if")) {if(p1>0) p1=p2; else p1=p3;}
                   break;
                 default :rderrcode=unknownfunction;;
               }
   } 
   if(rderrcode) return NULL;
   *(double*)args[0] = p1; return args[0];
}

char * userFuncTxt;
static double X;

static void*  rd_num_int(char* s)
{
   double    * p;
   int rc;
   p = (double *) malloc(sizeof(double)); 
   if(strcmp(s,"PI")==0) *p=M_PI;
   else if(strcmp(s,"x")==0) *p=X;  
   else       
   { rc= sscanf(s,"%lf",p);
     if(rc!=1)   {rderrcode=unknownidentifier;}
   }
   return (void*)p;     
}

double userFuncNumI(double x)
{ double *r;
  double rr;
  X=x;
  r = (double *)readExpression(userFuncTxt,rd_num_int,act_num,NULL); 
  rr=*r;
  free(r);
  return rr;
}  

double userFuncNumC(void)
{ double *r;
  double rr;
  r = (double *)readExpression(userFuncTxt,rd_num,act_num,NULL); 
  if(!rderrcode)
  { rr=*r;
    free(r);
    return rr;
  }else return 0;  
}  



