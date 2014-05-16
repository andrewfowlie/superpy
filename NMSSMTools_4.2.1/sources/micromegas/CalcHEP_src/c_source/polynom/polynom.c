/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include <limits.h>

#include "syst.h" 
#include"getmem.h"
#include "polynom.h"
#include "spinor.h"
#include "tensor.h"
#include <string.h>

 poly    garbage;
 poly   *contracts;
 int    monomLength,maxLength;

static int symb_size=0;

void (*memoryInfo) (int)=NULL;

void  delunit(void * p){((poly)p)->next = garbage; garbage = (poly)p;}

void  makeNewGarbage(void)
{  int brest;
   poly q;
   
   garbage=getmem_(symb_size);
   brest = blockrest(symb_size);

   for(q=garbage; brest; brest--, q=q->next) q->next=(poly) getmem_(symb_size);
   q->next=NULL;
   
   if(memoryInfo) (*memoryInfo)(usedmemory);
}

void  delpoly(poly* p)
{  poly   m, mm; 
   if (*p == NULL) return;
   m = *p; 
   
   mm = m->next; 
   while (mm != NULL) { m = mm; mm = mm->next; } 
   m->next = garbage; 
   garbage = *p; 
   *p = NULL; 
} 

poly  plusone(void)
{  poly  p;
   int i;
   NewUnit(p);
   p->next = NULL;
   p->num = NUM_ONE;
   for (i = 0; i < monomLength ; i++)  p->power[i] = 0;
   return p;
}

poly  copypoly(poly p)
{  poly  pp, qq,  copypoly1; 
   int  i; 

   if (p == NULL) return NULL; 
   else 
   { 
      NewUnit(pp);
      copypoly1 = pp; 
label_1: 
      for(i=0; i<monomLength; i++) pp->power[i] = p->power[i]; 
      pp->num = p->num; 
      p = p->next; 
      if (p) 
      { 
         NewUnit(qq);
         pp->next = qq; 
         pp = qq; 
         goto label_1;
      } 
      pp->next = NULL; 
   } 
   return copypoly1; 
} 


void  sewpoly(poly*p1, poly*p2)
{  poly  m, mm, m1, m2;
   monom ans_str;

   if (*p2 == NULL) return;
   if (*p1 == NULL) { *p1 = *p2; *p2 = NULL; return;} 

   m1 = *p1; 
   m2 = *p2; 

   m = &ans_str; 
   ans_str.next=NULL;
   *p2 = NULL; 

   for(;;)
   { int neq;

      unsigned long *pw1=m1->power, *pw2=m2->power, *pw_end=pw1+monomLength;
      for(; *pw1==*pw2 &&  pw1!=pw_end ;pw1++,pw2++ ){;}

      if(pw1==pw_end) neq=0; else if(*pw1>*pw2) neq= 1; else neq= -1;

     if(neq>0) 
     { 
        m->next = m1; 
        m=m1;
        m1 = m1->next; 
        if(!m1) {m->next = m2; break;} 
     } else if(neq<0)
     {
        m->next = m2;
        m=m2;
        m2 = m2->next;               
        if (!m2) {m->next = m1; break;}   
     }  else
     {  
        m1->num += m2->num; 
        mm = m2; 
        m2 = m2->next; 
        delunit(mm); 

        if (m1->num == NUM_ZERO)
        { 
           mm = m1; 
           m1 = m1->next; 
           delunit(mm); 
           if(!m1) { m->next = m2; break;} 
           m->next = m1; 
        } else 
        {
           m->next = m1;   
           m=m1;
           m1 = m1->next;  
           if(!m1) {m->next = m2; break;}                        
        }
        if(!m2) break;
      }  
   } 
   *p1 = ans_str.next;  
} 

void  multpolyint(poly* p,long i)
{
   if(i)
   { poly pp=*p;
     for(pp=*p;  pp;  pp = pp->next) pp->num *= i;
   } else delpoly(p); 
} 

static poly  multpolymono(poly plnm,poly mono)
{  poly  pp, qq, multpolymono1; 
   int   i; 

   if (plnm == NULL || mono == NULL) return NULL; 
   NewUnit(pp);
   multpolymono1 = pp;

label_1: 
   pp->num = plnm->num * mono->num; 
   for (i = 0; i < monomLength; i++) 
      pp->power[i] = plnm->power[i] + mono->power[i]; 
   plnm = plnm->next; 
   if (plnm != NULL) 
   {  
       NewUnit(qq);
      pp->next = qq; 
      pp = qq; 
      goto label_1;
   } 
   pp->next = NULL; 
   return multpolymono1; 
} 

/* -------------------------------------------------- */ 

poly  multtwopoly(poly q1,poly q2)
{  poly    mlttwpl, mltplmn, p1, p2; 
   mlttwpl = NULL; 

 
   if (q1 != NULL && q2 != NULL) 
   { 
      p1 = q1; 
      p2 = q2; 
      for(;;) 
      {
         p1 = p1->next; 
         p2 = p2->next; 
         if (p1 == NULL) { p1 = q2; p2 = q1; break;} 
         if (p2 == NULL) { p1 = q1; p2 = q2; break;} 
      }

      mlttwpl = multpolymono(p1,p2); 
      for(p2 = p2->next;   p2; p2 = p2->next) 
      { 
         mltplmn = multpolymono(p1,p2); 
         sewpoly(&mlttwpl,&mltplmn);
      } 
   } 
   return mlttwpl; 
} 



/* ---------- Common ----------- */ 

poly  scalarmult(int p1,int p2)
{unsigned         c, cc, n; 

   if (p1 < p2) 
   	{ c = -p1; cc = -p2; } 
   else 
   	{ c = -p2; cc = -p1; }
   n = cc + c * (c - 1) / 2; 
   return contracts[n-1]; 
}   /*  ScalarMult  */ 

void  assignsclmult(int p1,int p2,poly p)
{  int c, cc, n; 

   if (p1 < p2) 
      { c = -p1; cc = -p2; } 
   else 
      { c = -p2; cc = -p1; } 
   n = cc + c * (c - 1) / 2; 
   contracts[n-1] = p; 
}   /*  ScalarMult  */

void  deltensor(tensor * t)
{tensor   m, mm; 

   if (*t == NULL) return;
   mm = *t;

   do
   {
      m = mm;
      delpoly(&m->re);
      delpoly(&m->im);
      mm = m->next;
   }  while (mm != NULL);
   m->next = (tensor)garbage;
   garbage =(poly)(*t);
   *t = NULL;
}

tensor  copytens(tensor t)
{  tensor tt, qq, copytens1; 

   if(t == NULL) return NULL; 
   NewUnit(tt);
   copytens1 = tt; 

label_1: 
   memcpy(tt->tens,t->tens,maxIndex);
   tt->re = (poly)copypoly(t->re); 
   tt->im = (poly)copypoly(t->im);
   t = t->next; 
   if(t) 
   { 
      NewUnit(qq);
      tt->next = qq; 
      tt = qq; 
      goto label_1;
   } 
   tt->next = NULL; 
   return copytens1; 
} 


void  sewtens(tensor* t1,tensor* t2)
{ tensor         m, mm, m1, m2; 
  tensor_str ans_str;
 
   if (*t2 == NULL) return;
   if (*t1 == NULL) { *t1 = *t2; *t2 = NULL; return; } 
   m1 = *t1; 
   m2 = *t2; 

   ans_str.next=NULL;
   m=&ans_str;
   
   *t2 = NULL; 

   for(;;)
   { int neq;
   
     neq=memcmp(m1->tens,m2->tens,maxIndex);
     
     if(neq>0) 
     { 
        m->next = m1; 
        m=m1;
        m1 = m1->next; 
        if(!m1) {m->next = m2; break;} 
     } else if(neq<0)
     {
        m->next = m2;
        m=m2;
        m2 = m2->next;               
        if (!m2) {m->next = m1; break;}   
     }  else
     {  
        sewpoly(&m1->re,&m2->re);
        sewpoly(&m1->im,&m2->im);
                    
        mm = m2; 
        m2 = m2->next; 
        delunit(mm); 

        if (m1->re==NULL && m1->im==NULL)
        { 
           mm = m1; 
           m1 = m1->next; 
           delunit(mm); 
           if(!m1) { m->next = m2; break;} 
           m->next = m1; 
        } else 
        {
           m->next = m1;   
           m=m1;
           m1 = m1->next;  
           if(!m1) {m->next = m2; break;}                        
        }
        if(!m2) break;
      }  
   } 
   *t1 = ans_str.next;  
} 



void  multtensint(tensor * t,long i)
{ tensor      tt;

   if (i == 0)
      deltensor(t);
   else
   {
      tt = *t;
      while (tt != NULL)
      {
         multpolyint(&tt->re,i);
         multpolyint(&tt->im,i);
         tt = tt->next;
      }
   }
}

void  multtenspoly(tensor* t,poly p)
{
   if (p) 
   {  tensor tt;
      for(tt=*t;  tt; tt=tt->next) 
      {  poly pp=multtwopoly(tt->re,p); 
         delpoly(&tt->re); 
         tt->re = pp; 

         pp = multtwopoly(tt->im,p); 
         delpoly(&tt->im); 
         tt->im = pp; 
      } 
   } else deltensor(t); 
} 

void  multtensComplexpoly(tensor* t,poly re,poly im)
{  tensor tt;
   poly pRe,pIm,qq;

   if(re==NULL && im==NULL ) deltensor(t);  else 
   for(tt=*t; tt; tt = tt->next) 
   { 
     pRe = (poly)multtwopoly(tt->re,re);
     qq = (poly)multtwopoly(tt->im,im);
     multpolyint(&qq,-1);
     sewpoly(&pRe,&qq);

     pIm = multtwopoly(tt->re,im);
     qq =  multtwopoly(tt->im,re);
     
     sewpoly(&pIm,&qq);
     delpoly(&tt->re);  tt->re = pRe;
     delpoly(&tt->im);  tt->im = pIm;
   } 
} 

void tensRealPart(tensor * t)
{
  tensor first,*pred,tt;

  if(!t || !(*t)) return;
  tt=*t;
  first=tt;
  pred=&first;       
  while(tt)
  {
     delpoly( &(tt->im));        
     if(tt->re) { pred=&(tt->next); tt=tt->next;}
           else { *pred=tt->next;  delunit(tt); tt=*pred;}
  }
  *t= first;
}


void tensImPart(tensor * t)
{
  tensor first,*pred,tt;

  if(!t || !(*t)) return;
  tt=*t;
  first=tt;
  pred=&first;       
  while(tt)
  {
     delpoly( &(tt->re));        
     if(tt->im) { pred=&(tt->next); tt=tt->next;}
           else { *pred=tt->next;  delunit(tt); tt=*pred;}
  }
  *t= first;
}


void  symb_start(int nvar, varinfo * Vars, int nSpin, int nIndex,int nMom)
{
  unsigned long  z;
  int i,wp;

  monom M;
  tensor_str T;
  SpinTensor_str S;
  int width;
  char * e,*b;
 
  monomLength = 1;
  z = 1;
  for(i=0;i<nvar;i++)
  {
     if(z >= ULONG_MAX/Vars[i].maxdeg ) { monomLength++; z=Vars[i].maxdeg; }
     else   z*=Vars[i].maxdeg;

     Vars[i].wordpos = monomLength;
  }

   z = 1;
   wp = monomLength;
   for(i=nvar-1; i >= 0; i--) if (Vars[i].wordpos == wp)
      { Vars[i].zerodeg=z; z *= Vars[i].maxdeg; }
   else
      { Vars[i].zerodeg=1; z = Vars[i].maxdeg; wp--; }


   b=(char*) &M;
   e=(char*) &(M.power[monomLength]);
   symb_size=(e-b);
   
   maxIndex=nIndex;
   if(maxIndex == 0) maxIndex=1;
   tensLength=(maxIndex+ sizeof(long) -1)/sizeof(long);
             
   b=(char*)&T;
   e=(char*)&T.tens[maxIndex];   
   width=e-b;
   if(width>symb_size) symb_size=width;

   b=(char*)&S;
   e=(char*)&(S.g[nSpin]);
   width=e-b;     
   if(width>symb_size) symb_size=width;
         

   garbage = NULL;

   contracts=(poly*)getmem_((sizeof(poly)*nMom*(nMom+1))/2);
   
   for(i=0;i<nMom*(nMom+1)/2;i++) contracts[i]=NULL;
#ifdef STRACE 
	printf("vars position \n");
	for (i=0;i< nvar;i++)
	{  printf(" name= %s pos= %d maxdeg=  %lu zerodeg= %lu\n",
          Vars[i].name,Vars[i].wordpos,Vars[i].maxdeg,Vars[i].zerodeg);
	}
	printf("monomLength=  %d\n", monomLength);
#endif


} 
