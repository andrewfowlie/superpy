/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include "syst.h"
#include "polynom.h"
/*#include "../symb/include/test_wrt.h" */
#include "tensor.h"
#include "ps_tensor.h"
#include "sets.h"
#include <string.h>


void delEtens(Etens s)
{
   if(s)
   {  Etens    m,  mm;
      for(mm = s; mm;)
      {  
        m = mm; mm=m->next;
        deltensor(&(m->tcoef));
        delunit(m);
      }  
   }     
}        

Etens etens1(void)
{ Etens new;
  int i;
  
  NewUnit(new);
  new->next=NULL;
  new->tcoef=newtensor1();
  for(i=0;i<4;i++) new->eps[i]=X_MARK;
  return new;
}
Etens copyEtens(Etens s)
{
  Etens  pred, first=NULL;
  
  for(; s;s=s->next)
  {  Etens mm;
     NewUnit(mm);
     if(first) pred->next=mm;  else  first=mm;
     pred=mm;
     memcpy(&(mm->eps),&(s->eps),4);
     mm->tcoef=copytens(s->tcoef);
  }
  if(first) pred->next=NULL;
     
  return first;
}



void  addEtens(Etens* t1, Etens t2)
{  Etens    m, mm, m1, m2;
   Etens_str ans_str;
   
   if (t2 == NULL) return;
   if (*t1 == NULL) { *t1 = t2; return; }
   m1 = *t1;
   m2 = t2;

   ans_str.next=NULL;
   m = &ans_str;
    
   for(;;)
   { 
     char*c1=m1->eps;
     char*c2=m2->eps;
     char*c_end=c1+4;
     int neq=0;
 
     for(;*c1==*c2 && c1!=c_end;c1++,c2++){;}
     if(c1!=c_end){ if(*c1>*c2) neq=-1; else neq=1;}
     
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
        if(!m2) {m->next = m1; break;}   
     }else
     {
         sewtens(&(m1->tcoef),&(m2->tcoef));

         mm = m2;
         m2 = m2->next;
         delunit(mm);

         if(m1->tcoef==NULL)
         {
            mm = m1;
            m1 = m1->next;
            delunit(mm);
            if(m1) m->next = m1; else { m->next = m2; break;}
         } else
         {   m->next=m1;
             m=m1;
             m1 = m1->next;
             if (!m1) { m->next = m2; break;}
         }
         if (!m2) break;
     }
   }
   *t1 = ans_str.next;
}

void multEtensInt(Etens *t , long l)
{ Etens tt;
  if(l==0) { delEtens(*t); *t=NULL;} else 
  for(tt=*t;tt;tt=tt->next) multtensint(&(tt->tcoef),l);
}

void  multEtensPoly(Etens * t,poly p)
{  Etens   tt;
   if (p == NULL)  {delEtens(*t); *t=NULL;} 
   else  for (tt = *t;tt;tt = tt->next) multtenspoly(&(tt->tcoef),p);
}


void  multEtensTens(Etens * spn,tensor  tns)
{
  Etens s,sum,pred,sum0;
  char pos[120];
  int Ni,i,i0,l,m,n;
  char * c, *v;
  tensor t;  

/*
printf("multEtens\n");
writeEtens(*spn);
printf("\ntens\n");
writetens(tns);
*/  
  if(!tns) {delEtens(*spn); *spn=NULL; return;}
  if(!(*spn)) return;

  bzero(pos,maxIndex+1);
  
  v=tns->tens-1;
  for(i=1;i<=maxIndex;i++) if(v[i]) pos[i]++; 
  c=(*spn)->eps;
  for(i=0;i<4;i++) { n=c[i]; if(n>0 && n<=maxIndex) pos[n]++;}
  v=(*spn)->tcoef->tens-1;
  for(i=1;i<=maxIndex;i++) if(v[i]) pos[i]++;
  for(Ni=0,i=1;i<=maxIndex;i++) if(pos[i]>1) Ni++;  


  for(s=*spn,pred=NULL;s; )
  { t=multtwotens(s->tcoef,tns);
    deltensor(&s->tcoef);
    if(t) {s->tcoef=t;pred=s; s=s->next; } else 
    { Etens s_=s;
      s=s->next;
      if(pred) pred->next=s; else *spn=s;
      delunit(s_);
    } 
  } 

  if(!Ni) return;

  sum0=NULL;
  pred=NULL;
  sum=NULL;
  for(s=*spn;s; )
  { 
    Etens s_=s;
    s=s->next;
    s_->next=NULL;

    c=s_->eps;
    l=4;
           
    for(i=0;i<4;i++) {n=c[i]; if(n>0 && n<=maxIndex && pos[n]>1) break;}
       
    if(i==l) 
    { if(pred) pred->next=s_; else sum0=s_;
      pred=s_;
      continue;
    }
    i0=i;
    for(t=s_->tcoef; t;)
    { Etens s2;
      int kah=0;
      if(t->next) 
      { NewUnit(s2);
        s2->next=NULL;
        memcpy(&(s2->eps),&(s_->eps),4);
      } else s2=s_;
      
      s2->tcoef=t;
      v=t->tens-1;
      t=t->next;
      s2->tcoef->next=NULL;      
      c=s2->eps;

      for(i=i0;i<4;i++) 
      {  n=c[i]; 
         if(n>0)
         { m=v[n];
           if(m) 
           { v[n]=0; c[i]=m;
             if(m>0) v[m]=0;
             kah=1; 
           }
         }         
      } 
      if(kah)
      { int   sgn=1;
        for(i=0;i<3; )
        { 
          if(c[i]==c[i+1])
          {  deltensor(&s2->tcoef);
             delunit(s2);
             s2=NULL;
             break;
          }   
          if(c[i]<c[i+1])
          { char cc=c[i];
            c[i]=c[i+1];
            c[i+1]=cc;
            sgn*=-1;
            if(i)i--; else i++;
          } else i++;
        }
        if(s2&&sgn==-1) multtensint(&s2->tcoef,-1); 
      } 
      addEtens(&sum,s2);
    }
  }
  addEtens(&sum,sum0);    
/*
printf("\nresult\n");
writeEtens(sum);
printf("\nOK\n");  
*/
  *spn=sum;
} 

static tensor  multeps(char * e1,char *  e2)
{
  int i1,i2,j1,j2,k,L=4,sgn=1;
  char c1[4],c2[4];
  tensor sum=NULL;

  static struct epsdata
  {  int      contr[4];
     int  sgn;
  }  epsdata[24]=
  {
     {{0, 1, 2, 3}, 1}, {{1, 0, 2, 3},-1}, {{1, 2, 0, 3}, 1},
     {{2, 1, 0, 3},-1}, {{2, 0, 1, 3}, 1}, {{0, 2, 1, 3},-1},
     {{0, 1, 3, 2},-1}, {{1, 0, 3, 2}, 1}, {{1, 2, 3, 0},-1},
     {{2, 1, 3, 0}, 1}, {{2, 0, 3, 1},-1}, {{0, 2, 3, 1}, 1},
     {{0, 3, 2, 1},-1}, {{1, 3, 2, 0}, 1}, {{1, 3, 0, 2},-1},
     {{2, 3, 0, 1}, 1}, {{2, 3, 1, 0},-1}, {{0, 3, 1, 2}, 1},
     {{3, 1, 0, 2}, 1}, {{3, 0, 1, 2},-1}, {{3, 2, 1, 0}, 1},
     {{3, 1, 2, 0},-1}, {{3, 0, 2, 1}, 1}, {{3, 2, 0, 1},-1}
  };
  static int   epsnumb[5]  = {1, 1, 2, 6, 24};
  static int   epsfact[5]  = {24, 6, 2, 1, 1};

  for(i1=0,i2=0,j1=0,j2=0;i1<4&&i2<4&&e1[i1]>0&&e2[i2]>0;  )  
  {  if(e1[i1]==e2[i2]) { i1++; i2++; sgn=1&(sgn+i1-i2);L--; }
     else if(e1[i1]>e2[i2]) c1[j1++]=e1[i1++]; 
     else  c2[j2++]=e2[i2++];
  }   
  for(;i1<4;) c1[j1++]=e1[i1++];
  for(;i2<4;) c2[j2++]=e2[i2++];
  if(sgn) sgn=1; else sgn=-1;
      
  sum=NULL;
  for(k=0;k< epsnumb[L]; k++)
  {  tensor t=newtensor1();
     for(i1 = 0; i1<L&&t; i1++)
      {  char c,cc;
         c = c1[i1];
         cc = c2[epsdata[k].contr[i1]];
         if(c > 0) t->tens[c-1] = cc;
         if(cc> 0) t->tens[cc-1] = c;
         else if(c < 0) multtenspoly(&t,scalarmult(c,cc));
      }      
      multtensint(&t, -epsdata[k].sgn * epsfact[L] * sgn);
      sewtens(&sum,&t);
  }

  return sum;  
}


Etens mult2Etens(Etens S1, Etens S2)
{
  Etens s1,s2,sum,Scp,f1,f2;
  int i;
  tensor sumt;

/*
writeEtens(S1);
printf("\n");
writeEtens(S2);
printf("\n");
*/  
  sum=NULL;
  
  if(!S1||!S2) return NULL;
  
  if(S1->eps[0]==X_MARK) { f1=S1;S1=S1->next;} else f1=NULL;
  if(S2->eps[0]==X_MARK) { f2=S2;S2=S2->next;} else f2=NULL;

/*
if(S1&&S1->eps[0]==X_MARK) printf("1 X_MARK again\n");
if(S2&&S2->eps[0]==X_MARK) printf("2 X_MARK again\n");
*/

  if(f1&&f2) sumt=multtwotens(f1->tcoef,f2->tcoef); else sumt=NULL;

  if(f1) 
  { Scp=copyEtens(S2);
    multEtensTens(&Scp,f1->tcoef);
    addEtens(&sum,Scp);
  } 
  if(f2) 
  { Scp=copyEtens(S1);
    multEtensTens(&Scp,f2->tcoef);
    addEtens(&sum,Scp);
  } 


  for(s1=S1;s1;s1=s1->next) for(s2=S2;s2;s2=s2->next)
  { tensor t2,t1;

/*    char *c1=s1->eps, *c2=s2->eps;

printf("%d %d %d %d     %d %d %d %d\n", c1[0],c1[1],c1[2],c1[3],
                                         c2[0],c2[1],c2[2],c2[3]);
*/

    t1=multeps(s1->eps,s2->eps);
  
    t2=multtwotens(t1,s1->tcoef);
    deltensor(&t1);
    t1=multtwotens(t2,s2->tcoef); 
    deltensor(&t2);
    sewtens(&sumt,&t1);
  }
    
  if(sumt)
  { Etens sum0;
    NewUnit(sum0)
    sum0->next=sum;
    sum0->tcoef=sumt;
    for(i=0;i<4;i++) sum0->eps[i]=X_MARK;
    sum=sum0;
  }       
  return sum;
}
