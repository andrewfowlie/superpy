/*
 Copyright (C) 1997, Alexander Pukhov, e-mail pukhov@theory.npi.msu.su
*/
 
 
#include"model.h"
#include<stdlib.h>
#include<string.h>

 int       nmodelvar;
 varlist    modelvars=NULL;
 char *    EXTLIB=NULL;
 char *    EXTFunc=NULL;
 prtcl_base* prtclbase=NULL;
 prtcl_base* prtclbase1=NULL;
  
 int       nparticles;   /*  Number particles in model */
 algvertptr lgrgn=NULL;

void  locateinbase(char* name,int * number)
{int   i; 
 char name_[20];
 strcpy(name_,name); 
 i=strlen(name_); if(i) i--;
 if(name_[i]=='%') name_[i]=0; 

   for (i = 1; i <= nparticles; i++) 
      if (strcmp(name_,prtclbase[i-1].name) == 0) 
      {  *number = i; 
         return;
      }
   *number = 0; 
} 

int  pseudop(int np)
{ 
   return (np == 0 || np > nparticles || prtclbase[np-1].hlp == '*');
} 

int  fermionp(int p)
{ 
  return (prtclbase[p-1].spin % 2 == 1 && p <= prtclbase[p-1].anti);
} 

int  a_fermionp(int p)
{ 
  return (prtclbase[p-1].spin % 2 == 1 && p >= prtclbase[p-1].anti);
} 



int  bosonp(int p)
{ 
   return (prtclbase[p-1].spin % 2 == 0);
} 


int  vectorp(int p)
{ 
   return (prtclbase[p-1].spin == 2);
} 


int zeromass(int p)
{ 
   return (strcmp(prtclbase[p-1].massidnt,"0") == 0);
} 


int photonp(int p)
{ 
   return (vectorp(p) && zeromass(p) || pseudop(p));
} 

int ghostp(int p)
{ 
   return (prtclbase[p-1].hlp == 'c' || 
           prtclbase[p-1].hlp == 'C' || 
           prtclbase[p-1].hlp == 'f' );
} 

int gaugep(int j)
{ 
   return (j != 0 && prtclbase[j-1].hlp == 'G'); 
} 

int ghostmother(int  j)
{
   if (j == 0) return 0;
      switch (prtclbase[j-1].hlp)
      {
         case 'c': return j - 1;
         case 'C': return j - 2;
         case 'f': return j - 3;
         case 't': return j + 1;
         case 'T': return j + 2;
         default : return j;
      }
}
