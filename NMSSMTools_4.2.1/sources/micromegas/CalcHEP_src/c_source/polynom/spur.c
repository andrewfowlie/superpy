/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include "syst.h"
#include "polynom.h"
#include "tensor.h"
#include "ps_tensor.h"

#include "spinor.h"

int spinLength;

static int   nused;        /* from calcspur     */
static int   used[118];    /* from calcspur     */ 

static tensor calconespur(SpinTensor spnr);

Etens  calcspur(SpinTensor spnr1)
{  

   int     i, i1, i2, i3, i4, j; 
   int     c, cc, sgn; 
   char    ee[4]; 

   SpinTensor spnr;
   Etens     sum0,sum=NULL; 
   tensor    sumt=NULL; 
   
   for(spnr=spnr1;spnr;spnr=spnr->next) 
   { 
      if(1&spnr->l) continue; 
      if (spnr->g5 == 0) 
      {  tensor onesp,onesp2;
         nused = 0; 
         for(j=0; j<spnr->l; j++) used[j]=0; 
         onesp=calconespur(spnr);
         onesp2=multtwotens(onesp,spnr->tcoef);
         deltensor(&onesp);
         sewtens(&sumt,&onesp2);  
      } 
      else 
      {  Etens onesp=NULL;
         if(spnr->l < 4) continue; 
         else 
         { 
            onesp = NULL; 
            for(i1 = 0;    i1<spnr->l - 3; i1++) 
            for(i2 = i1+1; i2<spnr->l - 2; i2++) 
            for(i3 = i2+1; i3<spnr->l - 1; i3++) 
            for(i4 = i3+1; i4<spnr->l;     i4++) 
            { 
               nused = 4; 
               for(j=0; j<spnr->l; j++) used[j]=0; 
               ee[0] = spnr->g[i1]; used[i1] = 1; 
               ee[1] = spnr->g[i2]; used[i2] = 1; 
               ee[2] = spnr->g[i3]; used[i3] = 1; 
               ee[3] = spnr->g[i4]; used[i4] = 1; 
               j = 1; 
               sgn = ((i1 + i2 + i3 + i4) & 1) == 0 ? 1 : -1; 
               do 
               { 
                  c = ee[j-1]; 
                  cc = ee[j]; 
                  if(cc < c) j++; 
                  else if(cc > c) 
                  { 
                     ee[j-1] = cc; 
                     ee[j] = c; 
                     sgn = -sgn; 
                     if(j > 1) j--; else j++; 
                  }  else  sgn = 0; 
               } while (!(j == 4 || sgn == 0)); 
               if (sgn) 
               { Etens onesp2;
                 tensor coef = calconespur(spnr); 
                 if(!coef) continue;                 
                 if(sgn == -1) multtensint(&coef,-1);
                 NewUnit(onesp2);
                 onesp2->tcoef= coef;
                 for(;coef;coef=coef->next)
                 {  poly qq=coef->re;
                    coef->re=coef->im;
                    coef->im=qq;
                 }
                 onesp2->next=NULL;
                 for(i=0; i<4; i++) onesp2->eps[i] = ee[i]; 
                 addEtens(&onesp,onesp2);
               } 
            }
            multEtensTens(&onesp,spnr->tcoef);
            addEtens(&sum,onesp);      
         } 
      } 
   }
   if(sumt)
   { 
      NewUnit(sum0);
      sum0->next=sum;
      sum0->tcoef=sumt;
      for(i=0;i<4;i++) sum0->eps[i]=X_MARK; 
      sum=sum0;
   } 
   return sum;
}  


static tensor  calconespur(SpinTensor spnr) 
{  int    nu, sign, si, sk, i, k; 
   tensor ans, r1, r2; 

   nu = spnr->l; 
   if(nu == nused) return newtensor1();

   ans = NULL; 
   for(k =0; used[k]; k++){;} 
   used[k] = 1; 
   nused += 2; 
   sk = spnr->g[k]; 
   sign = -1; 

   for (i = k + 1; i <spnr->l; i++) 
   if (!used[i]) 
   {  
      sign = -sign;
      si = spnr->g[i]; 
      used[i] = 1; 
      r1 = calconespur(spnr); 
      used[i] = 0; 

      if (!r1) continue;
       
      if (sk < 0 && si < 0)  multtenspoly(&r1,scalarmult(sk,si)); 
      else 
      { 
         r2 = r1; 
         while (r2 != NULL) 
         { 
            if (si > 0) r2->tens[si-1] = sk; 
            if (sk > 0) r2->tens[sk-1] = si; 
            r2 = r2->next; 
         } 
      } 
      if (sign == -1) multtensint(&r1,-1); 
      sewtens(&ans,&r1); 
   }
   used[k] = 0; 
   nused -= 2; 
   return ans; 
} 
