/*
 Copyright (C) 1997, Alexander Pukhov, e-mail pukhov@theory.npi.msu.su 
*/
#include"chep_crt.h"
#include "syst2.h"
#include "physics.h"
#include "s_files.h"
#include "sortarr.h"
#include "process.h"
#include <unistd.h>


#include "constr.h"

#define maxref   (MAXINOUT - 2)
#define lref     (2 * maxref - 2)

#define indexlink struct indexStruct *
     typedef struct indexStruct {
                 particleNumType  outlist[MAXINOUT];
                 unsigned         num;
                 indexlink    ilink;
                          } indexStruct;
#undef indexlink
     typedef struct indexStruct *indexlink;

     typedef particleNumType listprtcl[lref];

#define elementlink struct element *
     typedef struct element {
                listprtcl    prtcl;
                elementlink  next;
                            } element;
#undef elementlink
     typedef struct element *elementlink; 

static  FILE *       bufres;   
static  int          nprimary, n_two;
static  unsigned     kmenu; 
static  indexlink    head; 
static  int          ndecay; 
static  unsigned     n_diagram; 
static  unsigned     m_diagram; 
static  elementlink  first, zero; 

static  whohow     inclp;
 

typedef struct procListStr
{ struct procListStr *  next;
  short  p[MAXINOUT];
}procListStr;
 
static procListStr*  genList=NULL;

static int errorcode=0;

static void   addprtcl(whohow p_list,int n)
{int i;
   for (i = 0; i < whohowMAX-1; i++)
   {
      if (p_list[i].who == n) { p_list[i].how++; return; } 
      if (!p_list[i].who)  
      { p_list[i].who = n;  p_list[i].how = 1; p_list[i+1].who = 0; return; }
   }
}


static void   delprtcl(whohow p_list,int n)
{int  i=1;
   while (p_list[i-1].who != 0)
   {
      if (p_list[i-1].who == n)
		{
         if (--p_list[i-1].how == 0)
         do {
               p_list[i-1].who = p_list[i].who;
               p_list[i-1].how = p_list[i].how;
					i++;
				}  while (!(p_list[i-1].who == 0));
         return;
      }
      i++;
   } 
} 


static void cleanProcList(void)
{ while(genList){procListStr*run=genList;genList=genList->next;free(run);} }


static int checkGenList(short * out_x)
{ procListStr* run;
  short new[MAXINOUT]; 
  int i,j,k;
 
  if(!n_x) return 1;
  if(nin==1)
  { for(i=0;i<n_x;i++) 
    if(out_x[i]==nprimary || out_x[i]==prtclbase[nprimary-1].anti) return 0;
  }

  new[0]=nprimary;
  if(nin==2) new[1]=prtclbase[n_two-1].anti;
  i=nin; 
  for(k=0;inclp[k].how;k++) for(j=0;j<inclp[k].how;j++) new[i++]=inclp[k].who;
  for(j=0; j<n_x;j++) new[i++]=out_x[j];
  { short * q=new+nin;
     SORTARR(q,nout);
  }

  if(genList)
  for(i=0;i<nin;i++) if(genList->p[i]!=new[i]){cleanProcList(); break;} 
  
  if(!genList) 
  { genList=malloc(sizeof(procListStr));
    genList->next=NULL;
    for(i=0;i<nin+nout; i++) genList->p[i]=new[i];
    return 1;
  }
  
  for(run=genList;run;run=run->next)
  { 
    for(i=0;i<nin+nout;i++) if(run->p[i]!=new[i])break;
    if(i==nin+nout) return 0;
  }
  run=genList;
  genList=malloc(sizeof(procListStr));
  genList->next=run;
  for(i=0;i<nin+nout;i++)genList->p[i]=new[i];
  return 1;
}


/* static  elementlink  ref[127][maxref];   */

typedef  elementlink  ref_arr[maxref];  
static  ref_arr * ref;

static void clearref(void)
{ 
   elementlink  e1, e2; 
   int          i, j; 

   for (j = 1; j <= nparticles; j++) 
   { 
      ref[j-1][0] = NULL; 
      for (i = 2; i <= ndecay - 1; i++) 
      { 
         e1 = ref[j-1][i-1]; 
         ref[j-1][i-1] = NULL; 
         if (e1 != zero) 
             while (e1 != NULL) 
             { 
                e2 = e1; 
                e1 = e1->next; 
                free(e2); 
             } 
      } 
   } 
} 

static int  testin(int l,decayDiagram restmp)
{
   int     nneed, np;
   int     i, j;
   int     copyincl[whohowMAX];
   int     copyins[whohowMAX];
   int     copyout[whohowMAX];    

   nneed = l - n_x;   /*  N_x  number of unknoun inclusev particles
         L   number of decay paticles */

   for(j=0; inclp[j].who; j++)  copyincl[j] = inclp[j].how;
   for(j=0; liminsp[j].who; j++) copyins[j] = liminsp[j].how;
   for(j=0; limout[j].who; j++) copyout[j] =  limout[j].how;
   
   for (i = 1; i <= 2 * (l - 1); i++)
   {
      np = restmp[i-1];
      if (np > 0)
      {  if (pseudop(np)) return 0;
         j = 1;
         while (inclp[j-1].who != 0 && inclp[j-1].who != np) j++;
         if (inclp[j-1].who == 0 || copyincl[j-1] == 0)
         { int k;
           for(k=0;limout[k].who;k++)
           if(limout[k].who==np) 
           { copyout[k]--;
             if(copyout[k]==0) return 0; else break;
           } 
	 }
         else { copyincl[j-1]--; nneed--; }
      }
      if (np < 0)
      {
         np = -np;
         if (np > prtclbase[np-1].anti) np = prtclbase[np-1].anti;
	 for(j=0;liminsp[j].who;j++)
	 {
           if (liminsp[j].who == np)
           { if(copyins[j]>0) { if (copyins[j] == 1) return 0; else copyins[j]--;}
	     else             { if (copyins[j]<-1) copyins[j]++;} 
	   }  
	 }
      }
   }

   if(l==ndecay) 
   {
     for(j=0;liminsp[j].who;j++) if(copyins[j]<-1) return 0;
     int nVV=0,pdg;
     for (i = 1; i <= 2 * (l - 1); i++)
     { np = restmp[i-1];
       if (np < 0)
       {  np = -np-1;
          pdg=prtclbase[np].N;
          if(pdg==23 || abs(pdg)==24)nVV++;               
       } 
     }
     
     if(nVV > ZWmax || nVV < ZWmin) return 0;
   }
   
   return (nneed <= 0);
}

static void  conon(decayDiagram diagram,int k)
{ 
   decayDiagram     m_diagram; 
   int     c; 
   int         l, li, i, length, shift; 

   l = 1; 
   do
   {  
      li = l; 
      c = 0; 
      do 
      {   
         if (++li == k) goto label_1;
         if (diagram[li-1] > 0) c++; else c--; 
      }  while (c != 1);
      for (i = l + 1; i <= li; i++) m_diagram[i - l-1] = diagram[i-1]; 
      length = li - l; 
      if (diagram[li] != 0) 
      {  /*  3 - particle  vertex  */ 
         k -= length; 
         c = 0; 
         do 
         {  /*  until c=1  */ 
            if (diagram[++li-1] > 0) c++; else c--; 
            diagram[li - length-1] = diagram[li-1]; 
         } while (c != 1); 
         li -= length; 
         for (i = 1; i <= length; i++) diagram[li + i-1] = m_diagram[i-1]; 
      } 
      else   
      {  /*  4 - particle vertex  */ 
         li++; 
         c = 0; 
         do /*  until c=1  */ 
            if (diagram[++li-1] > 0) c++; else c--; 
         while (c != 1);

         if (li >= k)
         {  /*  K in Fragment N 2  */  
            shift = length + 1; 
            k -= shift; 
            for (i = l + 1 + shift; i <= li; i++) 
               diagram[i - shift-1] = diagram[i-1]; 
            li += - shift + 1; 
            diagram[li-1] = 0; 
            for (i = 1; i <= length; i++) diagram[li + i-1] = m_diagram[i-1]; 
         } 
         else
         {  /*  K in Fragment  N 3  */  
            for (i = l + 2 + length; i <= li; i++) 
            m_diagram[i - l - 2] = diagram[i-1]; 
            shift = li - l; 
            k -= shift; 
            c = 0; 
            do 
            {  /*  until c=1  */ 
               if (diagram[++li-1] > 0) c++; else c--;
               diagram[li - shift-1] = diagram[li-1]; 
            }  while (c != 1); 
            li += - shift + 1; 
            diagram[li-1] = 0; 
            for (i = 1; i <= shift - 1; i++) 
            diagram[li + i-1] = m_diagram[i-1]; 
         } /* end if */ 
      } 

      label_1: l++;
   }  while (k != l);
}

static void dooutres(decayDiagram res)
{
   int         lmax, m, k, l, i;
   int      cond;
   decayDiagram     copyres;
   decayDiagram     tmplist[MAXINOUT];


   if (nin == 2 )
   {
      lmax = 2 * ndecay - 1;
      m = 0;
      for (k = 2; k <= lmax; k++)
      {
         if (res[k-1] == n_two)
         {
            memcpy(copyres,res,sizeof(copyres));
            conon(copyres,k);
            if (m != 0) for (l = 1; l <= m; l++)
            {
               cond = 1;
               for (i = 1; i < lmax; i++)
                  cond = cond && (copyres[i] == tmplist[l-1][i]);
               if (cond) goto label_1;
            }
            n_diagram++;
            FWRITE1(copyres,bufres);
            m++;
            memcpy(tmplist[m-1],copyres,sizeof(*tmplist));
            label_1: ;
         }
      }
   }
   else
   {
      n_diagram++;
      f_write(res,sizeof(decayDiagram),1,bufres);
   }
}


static void dtc(int* incond,int l,decayDiagram restmp)
{
   int         i;
   decayDiagram    res; 

   *incond = testin(l,restmp);
   if (!(*incond) || l < ndecay) return;

   *incond = 0;
   res[0] = -nprimary;
   for (i = 1; i <= 2 * (ndecay - 1); i++) res[i] = restmp[i-1];
   dooutres(res);
	goto_xy(1,24);  print("%u",n_diagram); refresh_scr();
   if (n_diagram > m_diagram)
      { if ( mess_y_n(35,15,"Continue")) m_diagram += 500;  else  errorcode = -1;}
}


static void decay(int n_part, int l)
{
   decayDiagram     restmp; /* Used in testin,dtc */
   int      existence, equalcond, incond;
   decaylink    pp;
   int     pk, d, c;
   int         m, ni, li, k, jk, ik, i1, i2, i3,
                i2_min, i2_max, i3_min, i3_max,
                i[3], j[3];
   elementlink  el1, el2, el3, elk, newelement, oldelement;
   elementlink  el[3];
  /* * * * * * * * * * * * * * * * */

   if (l == 1)
   {
/*      k = 1;
      while (inclp[k-1].who != 0 && inclp[k-1].who != n_part) k++;
	if (inclp[k-1].who == 0 && missingmass - prtclbase[n_part-1].mass < 1.0E-5)
         ref[n_part-1][0] = zero;
      else
*/      
       ref[n_part-1][0] = first;
      return;
   }
   existence = 0;
   pp=prtclbase[n_part-1].top;

   while (pp != NULL)
   {
      if (pp->part[2] == 0)
        { i3_min = 0; i3_max = 0; }
      else
      {
         if (l == 2) goto label_100;
         i3_min = 1;
         i3_max = l - 2;
      }
      for (i3 = i3_min; i3 <= i3_max; i3++)
      {
         i2_max = l - i3 - 1; 
         if (pp->part[1] == pp->part[2] && i2_max > i3) 
            i2_max = i3; 
         if (pp->part[1] == pp->part[0]) 
            i2_min = (l + 1 - i3) / 2;
         else i2_min = 1; 
         for (i2 = i2_min; i2 <= i2_max; i2++) 
         { 
            i1 = l - i2 - i3; 
            if (i3 == 0) ni = 2; 
            else ni = 3; 
            i[2] = i3; i[1] = i2; i[0] = i1; 

     /*  New REORDER    beg  */ 
            for (k = 1; k <= 3; k++) j[k-1] = k; 
            k = 1;
            while (k < ni) 
            { 
               d = i[j[k-1]-1] - i[j[k]-1]; 
               if(d > 0 || (d == 0 && 
                  pp->part[j[k-1]-1]<pp->part[j[k-1]-1])) 
               { 
                  c = j[k-1]; 
                  j[k-1] = j[k]; 
                  j[k] = c; 
                  if (k == 1) ++(k); else --(k); 
               } 
               else ++(k); 
            } 
     /*  New REORDER     end  */ 
            for (k = 1; k <= ni; k++) 
            { 
               jk = j[k-1]; 
               ik = i[jk-1]; 
               pk = pp->part[jk-1]; 
               if(ref[pk-1][ik-1] == NULL)  decay(pk,ik);
               if(ref[pk-1][ik-1] == zero) goto label_101;
            } 
            if (i3 == 0) el3 = NULL; 
            else el3 = ref[(pp->part[2])-1][i3-1]; 
            do 
            {  /* until El3=NULL */ 
               el2 = ref[(pp->part[1])-1][i2-1]; 
               do 
               {  /*  until El2=NULL  */ 
                  el1 = ref[(pp->part[0])-1][i1-1]; 
                  do 
                  {  /* until El1=NULL */ 
                     el[0] = el1; el[1] = el2; el[2] = el3; 
                     li = 1; 
                     for (k = 1; k <= ni; k++) 
                     { 
                        jk = j[k-1]; 
                        pk = pp->part[jk-1]; 
                        ik = i[jk-1]; 
                        if (ik == 1) restmp[li++-1] = pk; 
                        else 
                        {
                           restmp[li++-1] = -pk; 
                           elk = el[jk-1]; 
                           for (m = 1; m <= 2 * ik - 2; m++) 
                              restmp[li++-1] = elk->prtcl[m-1]; 
                        } 
                        if (ni == 3 && k == 1) restmp[li++-1] = 0; 
                     } 
                     dtc(&incond,l,restmp); 
                     if (errorcode != 0) goto label_102;
                     if (incond) 
                     { 
                        newelement = (elementlink)m_alloc(sizeof(element));
                        for (m = 1; m <= 2 * (l - 1); m++) 
                           newelement->prtcl[m-1] = restmp[m-1]; 
                        if (existence) 
                           oldelement->next = newelement;
                        else 
                        { 
                           ref[n_part-1][l-1] = newelement; 
                           existence = 1; 
                        }
                        oldelement = newelement; 
                     } 
                     equalcond = el1 == el2 ? 1 : 0; 
                     el1 = el1->next; 
                  }  while (!(el1 == NULL || equalcond)); 
                  equalcond = el3 == el2 ? 1 : 0; 
                  el2 = el2->next; 
               }  while (!(el2 == NULL || equalcond)); 
                  if (i3 != 0) el3 = el3->next; 
            }  while (el3 != NULL); 
/*  memory optimisation ?=>	if (l == ndecay && maxavail() < 20000) clearref(); */
            label_101: ; 
         }  /* I2 FOR circl  */ 
      }  /*  I3 FOR circl  */ 
      label_100:   /*  UNTIL pp=Nil; */ 
      pp=pp->next; 
   }                  /*  while pp<>NULL */ 
   label_102: 
   if (existence) oldelement->next = NULL; 
   else 
   if (l < ndecay) ref[n_part-1][l-1] = zero;
} 

static void doindex(void)
{ 
    int         i, j, k, c; 
    unsigned         recno; 
    decayDiagram     res; 
    indexlink    next, old; 
    int      switch_; 
    whohow       p_list; 
/*    int         nout; */

   recno = 0; 
/*   nout = hadr2.how == 0 ? ndecay : ndecay - 1;  */
   old = head;  
   while(FREAD1(res,bufres)==1)   
   {  
      switch_ =(nin==1); /*  hadr2.how == 0 ? 1 : 0; */
      next = (indexlink)m_alloc(sizeof(indexStruct));
      old->ilink = next; 
      old = next; 
      
      nilprtcl(p_list); 
      i = 2; j = 1; 
      do 
      {   
         if (res[i-1] > 0) 
         {  if (switch_) 
            { 
              addprtcl(p_list,res[i-1]); 
              j++; 
            } 
            else  switch_ = 1;
         }    
         i++; 
      }  while (j <= nout); 

      j = 1; 
      while (inclp[j-1].who != 0)
      {  for (i = 1; i <= inclp[j-1].how; i++)
            delprtcl(p_list,inclp[j-1].who);
         j++;
      }
      j = 1; i = 1;
      while (p_list[j-1].who != 0)
      {
         for (k = 1; k <= p_list[j-1].how; k++)
            next->outlist[(i++)-1] = p_list[j-1].who;
         j++;
      }
      if (n_x > 1)
      {
         i = 1;
         do
         {  /*  until i=N_X  */
            c = next->outlist[i-1];
            if (c < next->outlist[i])
            {
               next->outlist[i-1] = next->outlist[i];
               next->outlist[i] = c;
               if (i > 1) i--; else i = 2;
            }
            else  i++;
         }  while (i != n_x);
      }
      next->num = recno++;
   } 
   next->ilink = NULL;
}


static indexlink  mark1, mark2, mark3; /* From sortindex */

static void sorttwoblocks(void)
{ indexlink    mark_1, mark_2;
  int      i, diff;

   mark_2 = mark2->ilink;
   do
   {
      mark_1 = mark1->ilink;
      i = 1;
      do
         diff = mark_1->outlist[i-1] - mark_2->outlist[i-1];
      while (!(++i > n_x || diff != 0));
      if (diff < 0)
      {
         /*  Reoder  */
         mark2->ilink = mark3->ilink;
         mark3->ilink = mark_1;
         mark1->ilink = mark_2;
         /*  Rename  */
         mark_2 = mark2;   /* Temporary */
         mark2 = mark3;
         mark3 = mark_2;
         mark_2 = mark_1;
      } 
      mark1 = mark1->ilink; 
   }  while (mark1 != mark2); 
} 

static void sortindex(int ndiagram)
{int  lblock, i, di, iend; 
   if (n_x == 0) return;
   lblock = 1; 
   while (lblock < ndiagram) 
   { 
      mark3 = head; 
      iend = 0; 
      while (iend + lblock < ndiagram) 
      { 
         mark1 = mark3; 
         mark2 = mark1; 
         for (i = 1; i <= lblock; i++) mark2 = mark2->ilink; 
         iend += lblock; 
         di = ndiagram - iend; 
         if (di > lblock) di = lblock; 
         mark3 = mark2; 
         for (i = 1; i <= di; i++) mark3 = mark3->ilink; 
         sorttwoblocks(); 
         iend += di; 
      } 
      lblock *= 2; 
   } 
} 


static shortstr   recor_;         /* From addbuf */
static long       firstrec, nsubc; /* From addbuf */
static whohow     outprtcls;      /* From addbuf */

static void  addrecordtomenu(void)
{ shortstr     recor; 
  int         i, j, wh, hw, len; 

   /*  Sorting , Will be removed    */ 
   i = 1; 
   while (outprtcls[i].who != 0) 
   { 
      if (outprtcls[i-1].who <= outprtcls[i].who) ++(i); 
      else 
      { 
         wh = outprtcls[i-1].who; hw = outprtcls[i-1].how; 
         outprtcls[i-1] = outprtcls[i]; 
         outprtcls[i].who = wh; 
         outprtcls[i].how = hw; 
         if (i > 1) --(i); else ++(i); 
      } 
   }    
   /*  End Sorting      */ 
   strcpy(recor,recor_); 
   i = 1; 
   while (outprtcls[i-1].who != 0) 
   { 
      for (j = 1; j <= outprtcls[i-1].how; j++)
      { 
        strcat(recor,prtclbase[outprtcls[i-1].who-1].name); 
        strcat(recor,",");
      }
      ++(i); 
   } 

   len = strlen(recor); 
   recor[len-1] = ' '; 
   
   wrt_menu(1,++kmenu,recor,0,0,nsubc,firstrec);
   firstrec += nsubc;
}   /*  AddRecordToMenu  */


static void addbuf(void)
{
    int         i, j;
    int         mem[MAXINOUT];
    indexlink    old, next;
    decayDiagram     res;
    int wrtcode;

   strcpy(recor_,prtclbase[nprimary-1].name);
   if (nin==2)
   {
      strcat(recor_,",");
      strcat(recor_,prtclbase[prtclbase[n_two-1].anti-1].name);
   }
   strcat(recor_," -> ");
   next = head->ilink;
   for (i = 1; i <= n_x; i++) mem[i-1] = next->outlist[i-1];
   nsubc = 0;
   wrtcode=checkGenList(next->outlist); 
   while (next)
   {
      for (i = 0; i < n_x; i++)
      if(next->outlist[i] != mem[i])  
      {  if(wrtcode)
         {  memcpy(outprtcls,inclp,sizeof(outprtcls)); 
            for (j = 0; j < n_x; j++) addprtcl(outprtcls,mem[j]); 
            addrecordtomenu();
         } 
         nsubc = 0; 
         for (j = 0; j < n_x; j++) mem[j] = next->outlist[j]; 
         wrtcode=checkGenList(next->outlist);
         break;
      }
       
      if(wrtcode)
      {  nsubc++;      
         fseek(bufres,sizeof(decayDiagram)*next->num,SEEK_SET);       
         FREAD1(res,bufres);
         {  adiagram result;
            memcpy((void*)result.dgrm0,(void*)res,sizeof(decayDiagram));
            result.delMark=0;
            result.nsub=kmenu+1; 
            FWRITE1(result,diagrp);
         }
      }    
      old = next; 
      next = next->ilink; 
      free(old); 
   }
   if(wrtcode)  
   { memcpy(outprtcls,inclp,sizeof(outprtcls)); 
     for (j = 0; j < n_x; j++) addprtcl(outprtcls,mem[j]); 
     addrecordtomenu();
   } 
} 

static int* prclist(long * power )
{
  int i;
  int num[MAXINOUT];
 
  long N=1,m;
  int l=nin+nout-n_x;
  int * list;

  for(i=0;i<l;i++) N*=hadrons[i].pow;
  list=malloc(l*N*sizeof(int));
   
  for(i=0;i<l;i++) num[i]=0;

/* fill */
  for(m=0;m<N;m++)
  {  int c;
     for(i=0;i<l;i++) list[m*l+i]=hadrons[i].parton[num[i]]; 
     for(c=l-1;c>=0;c--)
     {   num[c]++;  
         if(num[c]==hadrons[c].pow)num[c]=0; else break;
     }
  }

  if(nin==2) /* energy test */
#ifdef QQQ
  for(m=0;m<N;m++)
  { double min=0, mout=0;
    for(i=0;i<nin;i++) min+=prtclbase[list[m*l+i]-1].mass;
    for(i=nin;i<l;i++) mout+=prtclbase[list[m*l+i] -1].mass;
    if(nin==1) {if(min<=mout) list[m*l]=0;} else
    { if(sqrts<=min || sqrts<=mout) list[m*l]=0;}
  }
#endif 
;
  else /* forbid A->A+X */
  for(m=0;m<N;m++)
  for(i=nin;i<l;i++) if(list[m*l]==list[m*l+i]) {list[m*l]=0; break;}
 

/*internal reorder*/

  for(m=0;m<N;m++)
  if(list[m*l])
  { int *  arr=list+m*l+nin;
    int len = l-nin;
    SORTARR(arr,len);
  }

/*delete copies */
  { long  n0=1,m1,k;

    for(i=0;i<nin;i++)n0*=hadrons[i].pow;
    m1=N/n0;
    for(k=0; k<n0;k++)
    for(m=0; m<m1;m++) if(list[(k*m1+m)*l])
    { int * a=list+(k*m1+m)*l;
      int * b=list+(k+1)*m1*l-l;
      for(;b>a;b-=l)
      {  for(i=nin;i<l;i++) if(a[i]!=b[i]) break;
         if(i==l) b[0]=0;
      }
    }
  }

/*
  printf("PRINT4\n");
  for(m=0;m<N;m++) if(list[m*l]) 
  {
    printf("%s ",prtclbase[list[l*m]-1].name);
    if(nin==2) printf("%s ",prtclbase[list[l*m+1]-1].name);
    printf("->");
    for(i=nin;i<l;i++) printf("%s ",prtclbase[list[l*m+i]-1].name);
    if(n_x) printf("%d*x",n_x);
    printf("\n");
  }
*/ 
  *power=N;
  return list;
}



int construct(void)
{
   int          i,j,ndiagram,l;
   char buf_name[STRSIZ];
   int *list;
   long N,m;
   char format[2]={62,71};

   errorcode=0;   
   firstrec=0;
   ref=m_alloc(sizeof(ref_arr)*nparticles);    
   diagrp=fopen(DIAGRP_NAME,"wb");
   sprintf(buf_name,"%stmp%cbuf.res",pathtouser,f_slash);
   menup=fopen(MENUP_NAME,"wb");
   m_diagram = 500;
   kmenu = 0;
   f_write(format,2,1,menup);
   head = (indexlink)m_alloc(sizeof(indexStruct));
   first = (elementlink)m_alloc(sizeof(element)); first->next = NULL;
   zero = (elementlink)m_alloc(sizeof(element));  zero->next = NULL;

   for(i=0; i<nparticles; i++) for(j=0; j<maxref; j++) ref[i][j]=NULL;
   ndecay = nin + nout - 1;
   n_diagram = 0;
   list=prclist(&N);
   l=nin+nout-n_x;
   for(m=0;m<N;m++) if(list[m*l])
   { 
      nprimary = list[m*l];
      nilprtcl(inclp);
      for(i=nin;i<l;i++)  addprtcl(inclp,list[m*l+i]);
      
      if(nin==2) 
      {  n_two=prtclbase[list[m*l+1] -1].anti;
         addprtcl(inclp, n_two);
      }        

      bufres=fopen(buf_name,"wb");
      decay(nprimary,ndecay);
      if (nin == 2) delprtcl(inclp,n_two);
      clearref();
      ndiagram=(ftell(bufres))/sizeof(decayDiagram);
         

      if (ndiagram > 0)
      {
         fclose(bufres); bufres=fopen(buf_name,"rb");
         doindex();
         sortindex(ndiagram);
         addbuf();
      }
      fclose(bufres);
      if (errorcode) break;
   }
   free(list);
   cleanProcList();
   n_diagram=(ftell(diagrp))/sizeof(decayDiagram);
   fclose(diagrp);

   
   unlink(buf_name);
   subproc_f = kmenu;
   fclose(menup);
   free(head);
   free(first);
   free(zero);
   free(ref);
   if (n_diagram == 0) errorcode = -1; else errorcode = 0;
   
   for (i = 18; i <= 24; i++) { goto_xy(1,i); clr_eol(); }   
   return errorcode;
}
