/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include "syst2.h"
#include "physics.h"
#include "ghosts.h"
#include "process.h"
#include "prepdiag.h"

vertexhlp   vertexes[2 * maxvert];
linkhlp     massindpos[5*maxvert];
fermloopstp fermloops[maxvert];
set         setmassindex,setmassindex0;
int         nloop;
int         fermmap[2*maxvert];
char        inoutmasses[MAXINOUT][VAR_NAME_SIZE];
momsum      momdep[3*maxvert];

int consLow=1;

void  coloringvcs(hlpcsptr currentghst)
{  int i, j;
   for(i=0; i<vcs.sizet; i++) for(j=0; j<vcs.valence[i]; j++)
   {     
       vcs.vertlist[i][j].partcl += currentghst->hlpcs[i][j];
       switch (currentghst->hlpcs[i][j]) 
       {
         case ghostmark:
         case antighostmark:
         case sbosonmark:    vcs.vertlist[i][j].lorentz = 0;
       }
   }
   setmassindex=findmassindex();
                                              
} 

void  attachvertexes(void)
{ int   v, l; 
  arr4byte  vert; 

  for (v = 0; v < vcs.sizet; v++) 
  { 
    for(l=0; l<vcs.valence[v]; l++)  vert[l]=vcs.vertlist[v][l].partcl;      
    for(l=vcs.valence[v]; l<MAXVALENCE; l++)  vert[l]=0;
    vertinlgr(vert,v+1,vertexes[v].subst,&vertexes[v].lgrnptr); 
  }
}

void findReversVert(void)
{
   int v,l,lp,k,vin,lin,fin,fout,fl_len;
   for(v=0;v<vcs.sizet;v++) vertexes[v].r_vert=0;

   for(lp=0;lp<nloop;lp++)
   {  fl_len=fermloops[lp].len;
      for(k=0;k<fl_len;k++)
      {
         l=fermloops[lp].ll[k];
         v=fermloops[lp].vv[k];
         if (k>0) {vin=fermloops[lp].vv[k-1];
                   lin=fermloops[lp].ll[k-1];
                  }
            else  {vin=fermloops[lp].vv[fl_len-1];
                   lin=fermloops[lp].ll[fl_len-1];
                  }
         lin=vcs.vertlist[vin-1][lin-1].nextvert.edno;
         fin=1 ; while (lin !=vertexes[v-1].subst[fin-1]) fin++;
         fout=1; while (l   !=vertexes[v-1].subst[fout-1])fout++;
         vertexes[v-1].r_vert=(fin>fout);
      }
   }
}


static void  findinoutmasses(void)
{ int v,l;
  for(v=0; v<vcs.sizet; v++) for(l=0; l<vcs.valence[v]; l++) 
  { int mom=vcs.vertlist[v][l].moment;
    if(mom > 0 && mom <= nin + nout) 
    strcpy(inoutmasses[mom-1],prtclbase1[vcs.vertlist[v][l].partcl].massidnt); 
  }
} 


set  findmassindex(void)
{ int  v, l; 
  set index=set_constr(_E);
  for(v=0; v<vcs.sizet; v++) for(l=0; l<vcs.valence[v]; l++) 
  { edgeinvert * L=&vcs.vertlist[v][l];
    if(L->lorentz && L->link.vno<v && !photonp(L->partcl) 
    && !gaugep(L->partcl) && prtclbase1[L->partcl].spin==2 ) 
    {     
      set_add1(&index,L->lorentz); 
      massindpos[L->lorentz-1].vrt1= v+1;
      massindpos[L->lorentz-1].ln1 = l+1; 
      massindpos[L->lorentz-1].vrt2= L->link.vno+1; 
      massindpos[L->lorentz-1].ln2 = L->link.edno+1; 
    } 
  }
  return index;
} 


static void  vectorLn(int v,int *l)
{  int s;
   for(;*l>=0;(*l)--)
   { if(PLR_PRTCL&vcs.vertlist[v][*l].prop)continue;
     s=prtclbase1[vcs.vertlist[v][*l].partcl].spin;
     if((s==2||s==4)) return;
   }
}

static void nextFerm(int*v, int*l)
{
  int l1;
  l1=vcs.vertlist[*v][*l].link.edno;
  *v=vcs.vertlist[*v][*l].link.vno;
  for(*l=0;;(*l)++)
  if(*l!=l1 && prtclbase1[vcs.vertlist[*v][*l].partcl].spin&1) break;
}


static void  findfermcycles(void)
{  int  v, v1, l, l1, count;

   for(v=0; v<vcs.sizet; v++) fermmap[v] = 0;

   for(v=0,nloop=0; v<vcs.sizet; v++) if(fermmap[v] == 0)
   {  fermloopstp * FL=&(fermloops[nloop]);
      int vr,lr;
      for(l=vcs.valence[v]-1; l>=0; l--) 
          if(a_fermionp(vcs.vertlist[v][l].partcl)) break;
      if(l< 0) continue;
      count = 0;
      FL->g5 = 0;
      for(v1=v; !fermmap[v1]; nextFerm(&v1,&l),count++)
      {  int nInt=0;
         
         FL->vv[count]=v1+1;
         FL->ll[count]=l+1;
         if(strchr("LR",prtclbase1[vcs.vertlist[v1][l].partcl].hlp)
           || PLR_PRTCL&vcs.vertlist[v1][l].prop) FL->g5=1;
         fermmap[v1]=nloop+1;
         FL->spin[count]=prtclbase1[vcs.vertlist[v1][l].partcl].spin;
         FL->ind[count][0]=vcs.vertlist[v1][l].lorentz;
           lr=vcs.vertlist[v1][l].link.edno;
           vr=vcs.vertlist[v1][l].link.vno;
         FL->ind[count][1]=vcs.vertlist[vr][lr].lorentz; 
               
         for(l1=vcs.valence[v1]-1;l1>=0 ;l1--)  
         {  vectorLn(v1,&l1);
            if(l1>=0 && fermmap[vcs.vertlist[v1][l1].link.vno]==nloop+1) 
                     FL->intln[count][nInt++]=l1+1;
         }
         FL->nint[count]=nInt;         
      } 
      FL->len = count;      
      nloop++;      
   }
}


static void  findinnerverts(void)
{/* int  v , l, nl; 

   for (v = 1; v <= nloop; v++) strcpy(fermloops[v-1].invrt,""); 

   for (v = 1; v <= vcs.sizet; v++) 
      if (fermmap[v-1] == 0)
      { 
         l = vectorslot(v);
         if (l != 0)
         {
            nl = fermmap[vcs.vertlist[v-1][l-1].nextvert.vno-1];
            if (nl != 0)
            {
               for (l = l - 1; l >= 1; l--)
                  if (nl != fermmap[vcs.vertlist[v-1][l-1].nextvert.vno-1])
                     goto label_1;
               sbld(fermloops[nl-1].invrt,
                    "%s%c",fermloops[nl-1].invrt,v); 

label_1: ; 
            } 
         } 
      } 
*/
} 


static void  findsubst(int v,int l,char* subst)
{momsum      frontsubst;
 int        i, j, vv, ll; 

   subst[0] = 0; /* strcpy(subst,""); */
   if (/*(setof(inp,intrp,_E) & vcs.vertlist[v-1][l-1].prop) != setof(_E)*/
   
       vcs.vertlist[v-1][l-1].prop & (IN_PRTCL|OUT_PRTCL)
   
   ) 
      subst[++subst[0]] = vcs.vertlist[v-1][l-1].moment;
      /* sbld(subst,"%s%c",subst,vcs.vertlist[v-1][l-1].moment); */
   else 
      for (i = 1; i <= vcs.valence[v-1]; i++) 
         if (i != l)
         { 
            vv = vcs.vertlist[v-1][i-1].nextvert.vno; 
            ll = vcs.vertlist[v-1][i-1].nextvert.edno; 
            findsubst(vv,ll,frontsubst);
            for (j = 1; j <= frontsubst[0]; j++)
               subst[subst[0] + j] = frontsubst[j];
            subst[0] += frontsubst[0]; 
            /* sbld(subst,"%s%s",subst,frontsubst); */
         } 
}  /*  FindCond  */ 


static void  changesign(char* subst)
{int        i; 

   for (i = 1; i <= /*strlen(subst)*/ subst[0]; i++) 
      subst[i] = - subst[i]; 
} 


static void  optimsubst(int v,int l,char* subst)
{momsum      frontsubst, backsubst; 
 int        i, vv, ll; 

   findsubst(v,l,frontsubst); 
   vv = vcs.vertlist[v-1][l-1].nextvert.vno; 
   ll = vcs.vertlist[v-1][l-1].nextvert.edno; 
   findsubst(vv,ll,backsubst); 
   if (/*strlen(frontsubst) <= strlen(backsubst)*/
       frontsubst[0] <= backsubst[0]) 
      for (i = 0; i <= frontsubst[0]; i++) subst[i] = frontsubst[i];
      /* strcpy(subst,frontsubst); */
   else 
   {  for (i = 0; i <= backsubst[0]; i++) subst[i] = backsubst[i];
      /* strcpy(subst,backsubst); */
      changesign(subst); 
   } 
}   /*  OptimCond  */


static void  standartsubst(int v,int l,char* subst)
{  int  i, vv, ll, ch, flg1 = 0, flg2 = 0; 
 
   if (vcs.vertlist[v-1][l-1].moment == nin + nout) 
   {  subst[0] = 0;
      /* strcpy(subst,""); */
      for (i = 1; i <= nin; i++)
         subst[++subst[0]] = i;          
         /* sbld(subst,"%s%c",subst,i); */
      for (i = nin + 1; i <= nin + nout - 1; i++)
         subst[++subst[0]] = -i;
         /* sbld(subst,"%s%c",subst,-i); */
      return;
   }

   ch = nin + nout; 
   findsubst(v,l,subst); 
   for (i = 1; i <= subst[0]; i++) 
   {  if (subst[i] ==  ch) flg1 = 1;
      if (subst[i] == -ch) flg2 = 1;
   }
   if (flg1 && flg2) 
   { 
      vv = vcs.vertlist[v-1][l-1].nextvert.vno; 
      ll = vcs.vertlist[v-1][l-1].nextvert.edno;
      findsubst(vv,ll,subst); 
      changesign(subst); 
   } 
}   /*  StandartSubst  */ 


static void  findinternalmoments(int indep)
{int     l, v, vln;
 int     m;

   for (m = 1; m <= 3 * maxvert; m++) momdep[m-1][0] = 0;

   for (v = 1; v <= vcs.sizet; v++)
   {
      vln = vcs.valence[v-1];
      for (l = 1; l <= vln; l++)
      {
         m = vcs.vertlist[v-1][l-1].moment;
         if (m > 0)  
         {   if (indep) standartsubst(v,l,momdep[m-1]);
              else       optimsubst(v,l,momdep[m-1]);
         }
      }
   }
}


void  preperdiagram(void)
{
   findinoutmasses();
   findfermcycles();
   setmassindex0=findmassindex();
   findinnerverts();
   findinternalmoments(/*consLow*/   nin+nout<=5);
}
