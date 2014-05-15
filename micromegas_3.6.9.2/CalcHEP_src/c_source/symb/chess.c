/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include "syst.h"
#include "syst2.h"
#include "physics.h"

#include "chess.h"

vertinfostr	vertinfo[6*maxvert];
int  n_vrt;
int  prgcode[6*maxvert][2];


static int  ncode;


static int  setpower(set ss)
{  int sp, i;
   for(i=set_first(&ss,1),sp=0;i;sp++, i=set_first(&ss,i+1)){;}
   return sp;
}


static void split(set s,set *ss1,set *ss2)
{  int  nsub;
   set   s1_, s2_,s1,s2;
   int  w, w_;
   int  i, j, l, cross, cross_,dim;
   int *nextelem;

   dim=1<<(setpower(s)-1);
   
   nextelem=m_alloc(dim*sizeof(int));

   for(l=set_first(&s,1),w=0,nsub=0;l; )
   {  int l1;
      w += vertinfo[l-1].weight;
      l1=set_first(&s,l+1);
      if(l1) 
      {
         nextelem[nsub]=l;
         for(i=1; i<=nsub; i++)  nextelem[nsub+i]=-nextelem[nsub-i]; 
         nsub = 2 * nsub + 1; 
      } 
      l=l1;
   } 
   
   
   l = nextelem[0]-1; 
   w -= 2 * vertinfo[l].weight; 
   s2=set_constr(l+1,_E);
   
   s1=set_aun(s,s2); 
   cross = 0; 
   for(j=0; j<vertinfo[l].vlnc; j++) if(set_in(vertinfo[l].link[j],s1))cross++; 
   s1_=s1; 
   s2_=s2; 
   w_ = w; cross_ = cross; 
   for (i = 1; i < nsub; i++)
   { 
      l = abs(nextelem[i])-1; 
      if (nextelem[i] > 0) 
      { 
         set_del1(&s1_,l+1);
         for (j = 0; j < vertinfo[l].vlnc; j++) 
         { int ll=vertinfo[l].link[j];
           if(set_in(ll,s1_)) (cross_)++; else if(set_in(ll,s2_)) (cross_)--; 
         }   
         set_add1(&s2_,l+1); 
         w_ -= 2 * vertinfo[l].weight;
      } 
      else 
      { 
         set_del1(&s2_,l+1); 
         for (j = 0; j < vertinfo[l].vlnc; j++) 
         { int ll=vertinfo[l].link[j];
           if(set_in(ll,s2_)) (cross_)++; else if(set_in(ll,s1_)) (cross_)--; 
         }
         set_add1(&s1_,l+1); 
         w_ += 2 * vertinfo[l].weight; 
      } 

      if (MEMORY_OPTIM) 
      { 
         if(abs(w_) < abs(w) || (abs(w_) == abs(w) && cross_ < cross)) 
         {  s1=s1_; s2=s2_; w = w_;  cross = cross_; } 
      } 
      else 
      { 
         if(cross_ < cross || (cross_ == cross && abs(w_) < abs(w)))
         {  s1=s1_; s2=s2_;  w = w_;  cross = cross_; } 
      } 
   }
   free(nextelem);
   *ss1=s1;
   *ss2=s2;
}


static void  programer(set s)
{  set  s1, s2; 
   if(setpower(s) > 1) 
   { 
      split(s,&s1,&s2);
      prgcode[ncode][0] = set_first(&s1,1);
      prgcode[ncode][1] = set_first(&s2,1);
      ncode++; 
      programer(s1); 
      programer(s2);
   } 
} 

void makeprgcode(void)
{  set   ss; 
   ss=set_constr(1,UpTo,n_vrt,_E); 
   ncode = 0;    
   programer(ss); 
} 
