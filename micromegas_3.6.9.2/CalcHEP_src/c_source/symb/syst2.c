/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include "syst.h"
#include "syst2.h"





void lShift(char *  s,int  l)
{ int i;
  if (l>0)
  {  int  m=strlen(s);
     for(i=0;i<= m-l; i++) s[i]=s[i+l];
  }
  if (l<0)
  {    i=strlen(s);
		 while ( i >=0 ) {s[i-l]=s[i]; i--;}
		 for (i=0;i<-l;i++) s[i]=' ';
  }
}



void  revers(void ** list)
{ void **q, **p, **r;

   if (*list == NULL) return;
   r = (void * *)(*list);
   q = (void * *)(*r);
   *r = NULL;
   while(q != NULL)
   {  p = (void * *)(*q);
      *q = (void *)r;
      r = q;
      q = p;
   }
   *list = (void *)r;
}


