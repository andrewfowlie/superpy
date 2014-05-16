/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include "syst.h"
#include "polynom.h"
#include "tensor.h"
#include "spinor.h"
#include <string.h>

int spinLength;

static void  kahane(SpinTensor * spnrs)
{  int     m1, m2, nrevol, i, j, l, mI1=maxIndex+1; 
   SpinTensor  s1, s2; 
   char    pos[120],cc,*c,last,ind;   
/*   
printf("\nkahane,input="); writespinor(*spnrs);
printf("\n");
*/      
   s1 = *spnrs; 
   *spnrs=NULL;
   
   while(s1) 
   {  l=s1->l;
      if(l<2){s2=s1; s1=s1->next; s2->next=NULL; addSpin(spnrs,s2); continue;}
      c=s1->g;
       
      bzero(pos,mI1);
      for(m2=0;m2<l;m2++) 
      { ind=c[m2];
        if(ind>0) {m1=pos[ind]; if(m1) {m1--; break;} else pos[ind]=m2+1;}
      }
      if(m2==l)
      { s2=s1; s1=s1->next; s2->next=NULL; addSpin(spnrs,s2); continue;}
/*
printf("Convolution %d %d\n",m1,m2);
*/
      s1->l-=2; 
      last=c[m2-1];
      for(i=m2+1; i<l;i++) c[i-2]=c[i];
      c[l-2]=0; c[l-1]=0;
       
      nrevol=m2-m1-1; 
      if(!nrevol)  multtensint(&s1->tcoef,4);
      else  if(nrevol&1) /*  odd gamma case  */
      { 
         c[m1]=last; 
         for(i=m1+1,j=m2-2; i<j ;i++,j--) {cc=c[i]; c[i]=c[j]; c[j]=cc;}
         multtensint(&(s1->tcoef),-2); 
      }else              /* even gamma case  */
      {
         multtensint(&(s1->tcoef),2);
         NewUnit(s2); 
         s2->next=s1->next;
         s1->next=s2;
         s2->tcoef=copytens(s1->tcoef);
         memcpy(&(s2->g5),&(s1->g5),l);
         s2->g[m1]=last;
         c[m1]=c[m2-2];
         for(i=m1+1,j=m2-3; i<j ;i++,j--) {cc=c[i]; c[i]=c[j]; c[j]=cc;}
         c[m2-2]=last;         
      } 
   }
/*
   printf("kahane,ouput="); writespinor(*spnrs);
   printf("\n"); 
*/
}   


SpinTensor spin1(void)
{ SpinTensor new;
  NewUnit(new);
  new->next=NULL;
  new->tcoef=newtensor1();
  new->g5=0;
  new->l=0;
  return new;
}

void delSpin(SpinTensor s)
{
   if(s)
   {  SpinTensor    m,  mm;
      for(mm = s; mm;)
      {  
        m = mm; mm=m->next;
        deltensor(&(m->tcoef));
        delunit(m);
      }  
   }     
}        

SpinTensor copySpin(SpinTensor s)
{
  SpinTensor  pred, first=NULL;
  
  for(; s;s=s->next)
  {  SpinTensor mm;
     NewUnit(mm);
     if(first) pred->next=mm;  else  first=mm;
     pred=mm;
     memcpy(&(mm->g5),&(s->g5),s->l+2);
     mm->tcoef=copytens(s->tcoef);
  }
  if(first) pred->next=NULL;
     
  return first;
}



void  addSpin(SpinTensor* t1, SpinTensor t2)
{  SpinTensor    m, mm, m1, m2;
   SpinTensor_str ans_str;
   
   if (t2 == NULL) return;
   if (*t1 == NULL) { *t1 = t2; return; }
   m1 = *t1;
   m2 = t2;

   ans_str.next=NULL;
   m = &ans_str;
    
   for(;;)
   { int neq=0;
     if(m1->g5 != m2->g5) { if(m1->g5 > m2->g5) neq= 1; else neq= -1;}
     else if(m1->l != m2->l)   { if(m1->l > m2->l)   neq= 1; else neq= -1;}
     else if(m1->l) neq=memcmp(m1->g,m2->g,m1->l);

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

void multSpinInt(SpinTensor *t , long l)
{ SpinTensor tt;
  if(l==0) { delSpin(*t); *t=NULL;} else 
  for(tt=*t;tt;tt=tt->next) multtensint(&(tt->tcoef),l);
}



void  multSpinPoly(SpinTensor * t,poly p)
{  SpinTensor   tt;
   if (p == NULL)  {delSpin(*t); *t=NULL;} 
   else  for (tt = *t;tt;tt = tt->next) multtenspoly(&(tt->tcoef),p);
   
} 


void  multSpinTens(SpinTensor * spn,tensor  tns)
{
  SpinTensor s,sum,pred,sum0;
  char pos[120];
  int Ni,i,i0,l,m,n;
  char * c, *v;
  tensor t;  
/*
printf("------\nmultSpinTens arguments\n1:");
writespinor(* spn);
printf("\n2:");
writetens(tns);
printf("\n arguments OK\n");
*/


  if(!tns) {delSpin(*spn); *spn=NULL; return;}
  if(!(*spn)) return;

  bzero(pos,maxIndex+1);
  v=tns->tens-1;
  for(i=1;i<=maxIndex;i++)  if(v[i]<0) pos[i]++;
  l=(*spn)->l;
  c=(*spn)->g;
  for(i=0;i<l;i++) { n=c[i]; if(n>0) pos[n]++;}
  v=(*spn)->tcoef->tens-1;
  for(i=1;i<=maxIndex;i++) if(v[i]) pos[i]++;
  for(Ni=0,i=1;i<=maxIndex;i++) if(pos[i]>1) Ni++;  

  for(s=*spn,pred=NULL;s; )
  { t=multtwotens(s->tcoef,tns);
    deltensor(&s->tcoef);
    if(t) {s->tcoef=t;pred=s; s=s->next; } else 
    { SpinTensor s_=s;
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
    SpinTensor s_=s;
    s=s->next;
    s_->next=NULL;

    c=s_->g;
    l=s_->l;
           
    for(i=0;i<l;i++) {n=c[i]; if(n>0 && pos[n]>1) break;}
       
    if(i==l) 
    { if(pred) pred->next=s_; else sum0=s_;
      pred=s_;
      continue;
    }
    i0=i;
    for(t=s_->tcoef; t; )
    { SpinTensor s2;
      int kah=0;
      if(t->next) 
      { NewUnit(s2);
        s2->next=NULL;
        memcpy(&(s2->g5),&(s_->g5),s_->l+2);
      } else s2=s_;
      
      s2->tcoef=t;
      v=t->tens-1;
      t=t->next;
      s2->tcoef->next=NULL;
      c=s2->g;
      
      for(i=i0;i<l;i++) 
      {  n=c[i]; 
         if(n>0)
         { m=v[n];
           if(m) 
           { v[n]=0; c[i]=m;
             if(m>0) {v[m]=0; if(pos[m]>1) kah=1;}
           }
         }
      } 
      if(kah)kahane(&s2);
      addSpin(&sum,s2);
    }
  }
  addSpin(&sum,sum0);    
  *spn=sum;
} 


SpinTensor mult2Spin(SpinTensor S1, SpinTensor S2, int forspur)
{
  SpinTensor s1,s2,sum2,S1cp;
  int l2,i;
  char*c1,*c2;

  sum2=NULL;
  for(s2=S2;s2; )
  { SpinTensor s2_=s2;
    s2=s2->next;
    s2_->next=NULL;
    S1cp=copySpin(S1); 
    multSpinTens(&S1cp,s2_->tcoef); 
    l2=s2_->l;
    c2=s2_->g;

    for(s1=S1cp;s1;s1=s1->next)
    { 
      s1->g5=(s1->g5+s2_->g5)&1;
      c1=s1->g+s1->l;
      c2=s2_->g;
      l2=s2_->l;
      for(i=0;i<l2;i++) c1[i]=c2[i];
      if(s2_->g5&&(s1->l&1)) multtensint(&(s1->tcoef),-1);
      s1->l+=l2;
    }
    kahane(&(S1cp));
    addSpin(&sum2,S1cp);
  }
  return sum2;
}
