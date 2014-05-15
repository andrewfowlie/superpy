/*
 Copyright (C) 2003, Alexander Pukhov 
*/

#include"syst.h"
#include"sets.h"

set set_constr(int i, ...)
{  va_list v;
   int k;
   set r;

   va_start(v,i);

   for(k=0;k<SETLEN;k++) r.field[k]=0;
   
   while(i != _E)
   {  if(i>255)
      {  fprintf(stderr,"Too large setofb element %d!! Hit enter:",i);
         getchar(); return r ;
      }      
      if(i == UpTo) i=va_arg(v,int); else k=i;
      for(;k<=i;k++){int z=k-1; r.field[z >> 3] |= 1 << ((z & 7));}
      i=va_arg(v,int);
   }
   va_end(v);
   return r;
}

int set_in(int a,set sp)
{  

   if (a<1 || a>8*SETLEN)
   {  fprintf(stderr,"Too large %d setofb element!! Hit enter:",a);
/*      getchar();
      exit(-1);
*/  
     return 0;      
   }
   a--;
   return ((1 << (a&7)) & sp.field[a >> 3]) != 0;
}

set set_or(set a,set b)
{ set res;
  int k;
   for( k=0; k<SETLEN; k++) res.field[k] = a.field[k] | b.field[k];
   return res;
}

set set_aun(set a,set b)
{ set res;
  int k;
  for( k=0; k<SETLEN; k++) res.field[k] = a.field[k] & (~ b.field[k]);
  return res;
}

set set_and(set a,set b)
{ set res;
  int k;
  for( k=0; k<SETLEN; k++) res.field[k] = a.field[k] & b.field[k];
  return res;
}

int set_eq0(set a)
{  int k;
   for(k=0; k<SETLEN; k++) if(a.field[k]) return 0;
   return 1;
}


int set_first(set* a,int i)
{ int k,l,p;
  unsigned char cc;

  if(i>8*SETLEN) return 0;
  if(i<1) i=0; else i--;
  
  k=(i)/8;
  l=(i)&7;
  p=1<<l;

  if(a->field[k]<p) 
  { for(k++; k<SETLEN &&(a->field[k]==0) ; k++) {;}
    if( k==SETLEN) return 0;
    p=1;
    l=0;
  }
  for(cc=a->field[k]; ;p=p<<1,l++) if(p&cc) return 8*k+l+1;
}


void set_add1(set *a, int k) 
{ if(k<1||k>8*SETLEN) return;     
  k--; 
  a->field[k >> 3] |= 1 << (k & 7);
}

void set_del1(set *a, int k) 
{
  if(k<1||k>8*SETLEN) return;
  k--; 
  a->field[k >> 3] &= ~(1 << (k & 7));
}    


void set_print(set s)
{  int i,k;  
   printf("{");
   for(k=0;k<SETLEN;k++) for(i=0;i<8;i++) if(s.field[k]&(1<<i)) printf(" %d ",1+8*k+i);
   printf("}\n"); 
   
/*   for(k=0;k<SETLEN;k++) printf(" %u ",s.field[k]);
   printf("\n");
*/
}

int set_eq(set a,set b){ return !memcmp(&a,&b,SETLEN);} 


