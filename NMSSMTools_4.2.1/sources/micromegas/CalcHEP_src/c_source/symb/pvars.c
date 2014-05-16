/*
 Copyright (C) 1997,2006, Alexander Pukhov
*/
#include <limits.h>

#include "physics.h"
#include "polynom.h"
#include "sos.h"
#include "syst2.h"
#include "spinor.h"

#include "pvars.h"

#ifdef STRACE
#include "test_wrt.h"
#endif
 polyvars *vardef;


void increaseVars( polyvars * v)
{  int oldsize=ALIG(v->nvar);
   if(++v->nvar > oldsize) 
   v->vars=re_alloc(v->vars,ALIG(v->nvar)*sizeof(* v->vars));
}

void clearVars( polyvars* v)
{
  if(v->vars){ free(v->vars); v->vars=NULL;}
  v->nvar=0;
}  
  
void  unite_vardef(polyvars *vardef_s,polyvars *vardef)
{int  n, nn;
 char  s[STRSIZ];

   for (nn = 0; nn < vardef->nvar; nn++)
   {
      strcpy(s,vardef->vars[nn].name);

      n = 0;
      while (n < vardef_s->nvar &&
         strcmp(vardef->vars[nn].name,vardef_s->vars[n].name)) n++;
      if (n < vardef_s->nvar) vardef_s->vars[n].maxdeg =
            MAX(vardef_s->vars[n].maxdeg,vardef->vars[nn].maxdeg);
      else
      {  
         increaseVars(vardef_s);
         strcpy(vardef_s->vars[n].name,s);
         vardef_s->vars[n].maxdeg = vardef->vars[nn].maxdeg;
      }
   }
}


void  addvar(char* varname,int deg)
{int    n;

   n = 0;
   while (n < vardef->nvar && strcmp(varname,vardef->vars[n].name)) n++;
   if(n < vardef->nvar) vardef->vars[n].maxdeg += deg;
   else
   { 
      increaseVars(vardef);
      strcpy(vardef->vars[n].name,varname);
      vardef->vars[n].maxdeg = deg + 1;
   }
}


int scalarProductPos(int p1,int p2)
{
  if(p1>p2) {int pp=p1;p1=p2;p2=pp;}
  return nmodelvar + PPSHIFT + (p1-1) + ((p2-1)*(p2-2))/ 2;  
}


int  modelVarPos(char* s)
{int      bt;
   for(bt=0; bt<=nmodelvar;bt++)
   {  if(strcmp(modelvars[bt].varname,s) == 0)return bt; }
return 0;   
    save_sos(14);
}



void  sortvar( )
{  int  i, p1, p2;

   for (i = 0; i < vardef->nvar; i++)   /*  numeration  */
      if(strchr(vardef->vars[i].name,'.'))
      {
         p1 = vardef->vars[i].name[1] - '0';
         p2 = vardef->vars[i].name[4] - '0';
	 vardef->vars[i].num = scalarProductPos(p1,p2);
      } 
      else if(!strcmp(vardef->vars[i].name,"Helicity1"))vardef->vars[i].num=nmodelvar+2;
      else if(!strcmp(vardef->vars[i].name,"Helicity2"))vardef->vars[i].num=nmodelvar+3;
      else if(!strcmp(vardef->vars[i].name,"HelicityN1"))vardef->vars[i].num=nmodelvar+4;
      else if(!strcmp(vardef->vars[i].name,"HelicityN2"))vardef->vars[i].num=nmodelvar+5;
      else  vardef->vars[i].num = modelVarPos(vardef->vars[i].name);

   for(i=0;i< vardef->nvar; i++) 
           if(!strcmp(vardef->vars[i].name,"i")) { vardef->nvar--;break;}
   for(;i< vardef->nvar; i++) vardef->vars[i]=vardef->vars[i+1];
   
                         
   if (vardef->nvar > 1)   /*  Sorting  */
   { 
      i = 1;
      while (i < vardef->nvar)
      if (vardef->vars[i-1].num > vardef->vars[i ].num) i++;
      else
      { varinfo tmpv=vardef->vars[i-1];
        vardef->vars[i-1]=vardef->vars[i]; 
        vardef->vars[i]= tmpv;
        if (i == 1) ++(i); else --(i);
      }
   }
}
