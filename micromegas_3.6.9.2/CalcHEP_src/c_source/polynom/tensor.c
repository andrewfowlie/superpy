/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include "syst.h" 
#include "polynom.h"
#include "tensor.h"
#include <string.h>

int maxIndex,tensLength;


tensor newtensor1(void)
{ tensor t;
  NewUnit(t);
  t->next=NULL;
  t->re=plusone();
  t->im=NULL;
  memset(t->tens,0,maxIndex);
  return t;
}


static tensor  multvectors(char * c1, char * c2)   
{  poly         p, q; 
   int         i, nloop;

   char * cres;
   tensor tres;
   
   NewUnit(tres);

   cres=tres->tens-1;

   memset(tres->tens,0,maxIndex);
   for (i = 1; i <= maxIndex; i++) 
   if(c1[i])  cres[i]= c2[i] ? 123 : 121; 
      else    cres[i]= c2[i] ? 122 : 0;

   nloop = 0; 
   p = plusone(); 
   for (i = 1; i <= maxIndex; i++) if (cres[i] == 123) 
   {  int s1, s2;
      cres[i] = 0; 
      for(s1=i;;) 
      {      
        s1=c1[s1]; if(s1<0) break; if(cres[s1]==123) cres[s1]=0;else break;
        s1=c2[s1]; if(s1<0) break; if(cres[s1]==123) cres[s1]=0;else break;
      }

      if (s1 == i)  ++(nloop); 
      else 
      { 
          for(s2=i;;) 
         {  
           s2=c2[s2]; if(s2<0) break; if(cres[s2]==123) cres[s2]=0;else break;       
           s2=c1[s2]; if(s2<0) break; if(cres[s2]==123) cres[s2]=0;else break;
         }
         if(s1>0) 
         { if (s2 > 0) {cres[s1]=s2; cres[s2]=s1;} else cres[s1]=s2;} 
         else 
         { if(s2>0) cres[s2] = s1; else 
            {  q=scalarmult(s1,s2);
               if(!q) { delpoly(&p); delunit(tres); return NULL;}
               q = multtwopoly(p,q); 
               delpoly(&p); 
               p = q; 
            } 
         } 
      } 
   }
   for (i = 1; i <= maxIndex; i++) 
        if (cres[i] == 121) cres[i] = c1[i]; else 
        if (cres[i] == 122) cres[i] = c2[i];
         
   if(nloop) multpolyint(&p,1 << (2 * nloop)); 
   tres->re = p; 
   tres->im=NULL;
   tres->next = NULL; 
   return tres;
} 

tensor  multtwotens(tensor t1,tensor t2)
{  tensor  ans2,t_1,t_2,t1_,t2_; 

   if(!t1 || !t2) return NULL;

   for(t_1 = t1,t_2 = t2; t_1 && t_2; t_1 = t_1->next,t_2 = t_2->next )
   {  
      if(!t_1)   break;
      if(!t_2) {t_1 = t1;t1 = t2;t2 = t_1; break;}
   }

   for(ans2=NULL,t2_=t2; t2_; t2_=t2_->next)
   {  tensor ans1;
      for(ans1=NULL,t1_=t1; t1_;t1_=t1_->next)  
      {  tensor tres=multvectors(t1_->tens-1,t2_->tens-1);
        if (tres) 
        { 
           multtensComplexpoly(&tres,t1_->re,t1_->im); 
           sewtens(&ans1,&tres); 
        } 
      }
      
      multtensComplexpoly(&ans1,t2_->re,t2_->im); 
      sewtens(&ans2,&ans1);
   } 
   return ans2; 
}  
