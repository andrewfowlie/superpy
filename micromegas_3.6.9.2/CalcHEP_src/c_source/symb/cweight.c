/*
 Copyright (C) 1997, Alexander Pukhov
*/

#include "physics.h"
#include "process.h"
#include "syst2.h"
#include "colorf.h"
#include "cweight.h"

int NcInfLimit=0;

static int a2k(vampl * g, int nc, int * chains, int * nv, cvertex * vl)
{
  int   i,k, ne, en=0, maxnv=*nv; 
  int   maptar[maxvert][MAXVALENCE];
  int   endl[MAXINOUT];
  int   endc[MAXINOUT];  

  for(i=0;i<MAXINOUT;i++){endl[i]=0;endc[i]=0;}
    
  for(i=0;i<maxvert;i++) for(k=0;k<MAXVALENCE;k++) maptar[i][k]=0;
  *nv=0;
  for(i=0;i<g->size; i++)  if (typev(g->vertlist[i],g->valence[i]) != zv) 
  { 
     if(*nv>= maxnv) return 1;
     vl[*nv].e[0]=0; vl[*nv].e[1]=0;  vl[*nv].e[2]=0;
     vl[*nv].vt = typev(g->vertlist[i],g->valence[i]); 

     for (k=0, ne=0; ne<g->valence[i] ;ne++)
     {  int dim= prtclbase[g->vertlist[i][ne].partcl-1].cdim;
        if(dim!=1)
        {  int l=maptar[i][ne];
           if(!l)
           { 
               l = ++en; 
               if(g->vertlist[i][ne].link.vno!= nullvert)
               {  maptar[i][ne] = l;
                  maptar[g->vertlist[i][ne].link.vno]
                     [g->vertlist[i][ne].link.edno]=l;
               }
               else 
               {  int  k=g->vertlist[i][ne].link.edno;
                  endl[k]=l;
                  endc[k]=prtclbase[g->vertlist[i][ne].partcl-1].cdim; 
               }
            } 
         
            if (vl[*nv].vt != g3 && vl[*nv].vt != g2) 
            switch (dim) 
            {
              case  8:  vl[*nv].e[0]=l;  break; 
              case -3:  vl[*nv].e[1]=l;  break; 
              case  3:  vl[*nv].e[2]=l;  break; 
              default: return 2;
            } else vl[*nv].e[k++] = l;
         } 
      }  
      (*nv)++;   
  }

  for(i=0; i<MAXINOUT;i++) if(endc[i]==8)
  {  
     if(*nv>= maxnv) return 1;
     vl[*nv].vt=qg;
     vl[*nv].e[0]=endl[i]; 
     endl[i]=++en;
     vl[*nv].e[2]=en;
     vl[*nv].e[1]=en+1;
     en++; 
     (*nv)++;
  }
  
  for(i=0;i<nc;i++)
  {  int from=chains[2*i+1];
     int to  =chains[2*i];  

     if(*nv>= maxnv) return 1;

     vl[*nv].vt=tv;
     vl[*nv].e[0]=0; 

     vl[*nv].e[1]=endl[from];
     if(endc[to]==-3) vl[*nv].e[2]=endl[to];
     else             vl[*nv].e[2]=endl[to]+1; 
     (*nv)++;
  }
  return 0;
}  


static int maxNcPower(vcsect* g)
{
  int i,j;
  int power=0;

  for(i=0;i<g->sizel; i++)
  for(j=0;j<g->valence[i];j++) if (g->vertlist[i][j].link.vno >=g->sizel)
  { int np=g->vertlist[i][j].partcl;
    switch(prtclbase[np-1].cdim)
    { case  3:
      case -3: power++; break;   
      case  8: power+=2;break;
    }
  }
  return power/2;
}

static void getLeadingTerm(factor * f, int  maxP, long *num, long *den)
{

   if(maxP<f->len - f->dpow -1) {printf("BUG in my BRAIN \n"); sortie(58);}

   if(maxP>f->len - f->dpow- 1) {*num=0; *den=1;} else
   {  *num=f->nc[f->len-1];
      *den=f->dc;
      while (maxP--) (*num) *=3;
   }
}  


void c_basis_coef(vampl * g,int pow,int nc,int * chains,long * num,long * den)
{ int i,nv;
  cvertex   vl[3*MAXINOUT];
  factor * f;

  if(!pow) return;

  for(i=0; i<pow; i++)
  { 
     nv=3*MAXINOUT;
     a2k(g, nc, chains +2*nc*i, &nv, vl);

     f=colorFactor(nv,vl);
     getLeadingTerm(f,nc,num+i,den+i);

     free(f->nc);
     free(f);
  }
}


void cwtarg(vcsect* g)
{  factor * f;
   cvertex   vl[2*MAXINOUT];
   int nv=2*MAXINOUT;
   t2k2(g,&nv,vl);
   f=colorFactor(nv,vl);

   if(NcInfLimit) getLeadingTerm(f,maxNcPower(g),&(g->clrnum),&(g->clrdenum));
   else fct_num_calc(f,3,&(g->clrnum), &(g->clrdenum));

   free(f->nc);
   free(f);
}  /* CWTarG */


static void  lreduce(long * l1, long  * l2)
{  long    c, i1, i2;

   i1 = *l1; i2 = *l2;
   if(i1<0) i1=-i1; if(i2<0) i2=-i2;

   if (i2 > i1) { c = i1; i1 = i2; i2 = c;}
   while (i2 != 0) { c = i2; i2 =i1%i2; i1 = c; }
   (*l1) /= i1;
   (*l2) /= i1;
}


int generateColorWeights(csdiagram*csdiagr,int cBasisPower,int nC,int*cChains,
     long * cCoefN,long * cCoefD)
{
   vcsect vcs;
   int NcInfLimit_tmp=NcInfLimit;
   int i;

   transfdiagr(csdiagr,&vcs);
   cwtarg(&vcs);
      
   NcInfLimit=1;
      
   if(vcs.clrnum) 
   {  vampl left, right;
      long n=1,d=1;
      long * cCoefNr=malloc(cBasisPower*sizeof(long));
      long * cCoefDr=malloc(cBasisPower*sizeof(long));
 
      decompose(vcs,&left,&right);
      c_basis_coef(&left,cBasisPower,nC,cChains,cCoefN,cCoefD);
      c_basis_coef(&right,cBasisPower,nC,cChains,cCoefNr,cCoefDr);

      for(i=0;i<nC;i++) d*=3;

      for(i=0;i<nin+nout;i++) 
      { vertlink Q=left.outer[i];
        if(8==prtclbase[left.vertlist[Q.vno][Q.edno].partcl-1].cdim)n*=2;
      }

      for(i=0;i<right.size;i++)
      { int n8=0;
        int j;

        for(j=0;j<right.valence[i];j++)
        if(8==prtclbase[right.vertlist[i][j].partcl-1].cdim) n8++;
        if(n8==3) n*=-1;
      }

      for(i=0;i<cBasisPower;i++)
      { 
         cCoefN[i] *= cCoefNr[i]*n *vcs.clrdenum;
         cCoefD[i] *= cCoefDr[i]*d *vcs.clrnum;
         lreduce(cCoefN+i, cCoefD+i);
      }
      for(i=0,n=0,d=1;i<cBasisPower;i++)
      { 
         n=cCoefN[i]*d+n*cCoefD[i];
         d*=cCoefD[i];
         lreduce(&n,&d);
      }
/*      if(n!=vcs.clrnum || d!=vcs.clrdenum) printf("BUG in COLOR (%d %d) (%d %d)\n",
      vcs.clrnum,n, vcs.clrdenum,d  ); 
*/      
      free(cCoefNr); free(cCoefDr);       
   } else for(i=0;i<cBasisPower;i++){cCoefN[i]=0; cCoefD[i]=0;}

   NcInfLimit=NcInfLimit_tmp;
   return vcs.clrnum;
}

static int *used, *perm, *wrt, *pos3, *pos_3;
static int nc_;

static void recurGen(int k)
{
  int i;
  
  if(k==nc_) for(i=0;i<nc_;i++) { *(wrt++)=pos_3[i]; *(wrt++)=pos3[perm[i]];}
  else
  { 
    for(i=0; i< nc_; i++) if(!(used[i] || pos_3[i]==pos3[k]))
    {  used[i]=1;
       perm[i]=k;
       recurGen(k+1);
       used[i]=0;
    }
  }
}


int infCbases(int np, int * cweight, int *nc, int *pow, int ** chains)
{
   int n3=0, n_3=0;
   int i;

   pos3 =(int *)malloc(np*sizeof(int));
   pos_3=(int *)malloc(np*sizeof(int));

   for(i=0;i<np;i++) switch (cweight[i])
   { 
     case -3:  pos_3[n_3++]=i; break;
     case  1:  break;
     case  3:  pos3[n3++]=i; break; 
     case  8:  pos3[n3++]=i; pos_3[n_3++]=i; break;
     default:  return 1; 
   }
   if(n3 !=n_3)  return 2;
    
   nc_=n3;
   if(nc_)
   {  int pow_=1;
      for(i=2;i<=nc_;i++) pow_*=i;

      wrt=(int *)malloc(2*nc_*pow_*sizeof(int));
      *chains=wrt;

      used= (int *) malloc(nc_*sizeof(int)); for(i=0;i<nc_;i++) used[i]=0;
      perm= (int *) malloc(nc_*sizeof(int));   

      recurGen(0);

      free(perm); free(used);
      *pow= (wrt-(*chains))/(2*nc_);
      *chains= (int *) realloc(*chains, 2*nc_*(*pow)*sizeof(int));
      *nc=nc_;
   
   } else { *nc=0; *pow=0; *chains=NULL;}

   free(pos3), free(pos_3);
   return 0;
}
