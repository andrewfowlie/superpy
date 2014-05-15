/*
 Copyright (C) 1997, Alexander Kryukov
*/
/**********************************************************/
/*  CopyRight (C) 1990, SCL                               */
/*  Author        A.Kryukov                               */
/*  E-mail        kryukov@theory.npi.msu.su               */
/*  Version       4.61                                    */
/*--------------------------------------------------------*/
/*  Last Rev.     08/01/90                                */
/*                18/02/94    t2k - gluon transfer vert.  */
/*                19/03/99    findl fix tedpole bug       */
/**********************************************************/

#ifdef CALCHEP
#include "syst2.h"
#define malloc m_alloc
#else
#include<stdio.h>
#include<stdlib.h>
#include <memory.h>
#include <string.h>

#endif
#include "syst.h"
#include "colorf.h"

int MAXGLEN, SIZE, MAX_POW;
#define CERRLEV1 0    /*  Run time error level            */
#define CERRLEV2 0    /*  Halt level                      */ 

static int cerror(int n,char* s)
{ /*  - generate error message occur in color package - 08/01/90  */

   fprintf(stderr,"***** %s\n",s);
        if (n > CERRLEV1) fprintf(stderr,"Error  %u in  colorf.c \n",n), sortie(55);
   else if (n > CERRLEV2) sortie(56); 
   return 0;
}  

/*==============     FACTOR SESSION ==================*/

static void fct_init(factor * fct) 
{ fct->len=0; fct->dpow=0; fct->dc=1; 
  MAX_POW=10;
  fct->nc=malloc(MAX_POW*sizeof(long));
}

static void add_fct(factor *f, long sgn, long pow2,int powN,int powNN_1)
{ long dd;
  int i,j;
  if(!sgn) return;
   
  if (pow2>=0) sgn*=f->dc<<pow2  ; else
  {  pow2=1<<(-pow2);
     if(pow2>f->dc) {dd=pow2/f->dc; f->dc=pow2; for(i=0;i<f->len;i++) f->nc[i]*=dd;}
     else if(pow2<f->dc) sgn*=f->dc/pow2;
  }

  if(-powN>f->dpow)
  { dd=-powN-f->dpow; f->dpow=-powN;
    f->len+=dd; 
    if(f->len>MAX_POW) {MAX_POW=2*f->len; f->nc=realloc(f->nc,MAX_POW*sizeof(long));}
    for(i= f->len-1; i>=dd; i--) f->nc[i]=f->nc[i-dd];
    for(i=0;i<dd;i++)f->nc[i]=0;
    powN=0;
  } else powN=f->dpow+powN;

  if(f->len < powN+2*powNN_1+1)
  {  dd=powN+2*powNN_1+1 - f->len;
     f->len=powN+2*powNN_1+1;
     if(f->len>MAX_POW){MAX_POW=2*f->len; f->nc=realloc(f->nc,MAX_POW*sizeof(long));}
     for(i= f->len-dd; i<f->len;i++) f->nc[i]=0;
  }

  if(powNN_1&1) sgn*=-1;
  for(i=0,j=0; i<1+2*powNN_1; i+=2,j++)
  { f->nc[i+powN]+=sgn;
    sgn=  -sgn*(powNN_1-j)/(j+1);
  }

  for(dd=0; !f->nc[dd] && dd<f->len;dd++);
  if(dd)
  {
     for(i=0;i<f->len-dd; i++) f->nc[i]= f->nc[i+dd];
     f->len-=dd;
     f->dpow-=dd;
  }
  if(!f->len) {f->dpow=0; f->dc=1; f->len=1; f->nc[0]=0;}

  dd=f->dc-1;
  for(i=0;i<f->len && dd;i++)  while(f->nc[i]&(long)dd) dd/=2;
  if(dd){ dd++; for(i=0;i<f->len;i++) f->nc[i]/=dd; f->dc/=dd; }
}
/* ================= END OF FACTOR ==================== */

#define CDEBLEV  0    /*  Total debug level               */ 

/*################################################################*/
       
typedef struct cgraph 
   {  struct cgraph *    next;

      int          sgn, pow2, powN, powNN_1; /* factor */

      int          en;            /* Name of next edge */
      int          gl;            /* Number of vertecies (graph length) */
      cvertex       vl[2];         /* array of verticies */
   }  cgraph; 

/* ************************** Cross reference ************************* */ 
/* *                                                                  * */ 
/* *  GevV                                                            * */ 
/* *    +-----> CError                                                * */ 
/* *                                                                  * */ 
/* *  GetEN                                                           * */ 
/* *                                                                  * */ 
/* *  CError                                                          * */ 
/* *                                                                  * */ 
/* *  WrCG                                                            * */ 
/* *                                                                  * */ 
/* ******************************************************************** */ 


#  if (CDEBLEV > DEBLEV) 
static char  vtarr[5][4]  = {"ZV", "TV", "G2", "QG", "G3"};

static void wrcg(cgraph* cg)
{/*  Write C-graph on standard device - 04/01/90 */ 
   int      i; 
   printf("Kr: [(%ld*2^%ld) N^%ld (N^2-1)^%ld \n",
	   cg->sgn,cg->pow2,cg->powN,cg->powNN_1);
   printf("en=%d, gl=%d\n",cg->en,cg->gl);
   for (i = MAXGLEN; i >= 1; i--)
      if (cg->vl[i-1].vt != zv)
      {
	fprintf(stderr,"(%d/%s ",i,vtarr[cg->vl[i-1].vt-1]);
	fprintf(stderr,"%d,%d,%d)",
	 cg->vl[i-1].e[0],cg->vl[i-1].e[1],cg->vl[i-1].e[2]);
      }   /* if */
   fprintf(stderr,"]\n");
} 
#endif

static int getv(cgraph* cg)
{ /* return number of first free vertex in C-graph - 04/01/90 */ 
   int      ok; 
   int      n; 
    
   ok = 0; 
   n = 1; 
   while (n <= MAXGLEN && !ok)
      if (cg->vl[n-1].vt == zv)  ok = 1;  else  n++; 

   if (n > MAXGLEN) cerror(254,"GetV: no free vertex in C-graph"); 
   else  { ++(cg->gl);  return n; }  
   return 0;
}  


#define SINGL   1    /* Colour singlet */ 
#define TRIPL  -3    /* Colour triplet */ 
#define ATRIPL  3    /* Colour antitriplet */ 
#define OCTET   8    /* Colour octet   */ 
#define DEBLEV 10    /*  Debug level   */ 


static void initcg(cgraph* cg)
{ /* Initiate color graph - 04/01/90 */ 

   int  i; 
   cg->sgn=1;  cg->pow2=0;  cg->powN=0; cg->powNN_1=0;

   cg->en = 0; 
   cg->gl = 0; 
   for (i = 0; i < MAXGLEN; i++)
   {  cg->vl[i].vt = zv; 
      cg->vl[i].e[0] = 0; 
      cg->vl[i].e[1] = 0; 
      cg->vl[i].e[2] = 0; 
   } 
}


static int findv(vtype vt,cgraph* cg,int * n)
{ /* return True and number first vertex with type VT in C-graph - 06/01/90 */ 
 
  int      i = 1; 
    
  while (i <= MAXGLEN && cg->vl[i-1].vt != vt) i++;
  if (i > MAXGLEN) return 0;  else  *n = i; 
  return 1; 
} 


static int findl(int e1,int n,cgraph* cg)
{ /* - find number of vertex contane edge e1 - 08/01/90  */ 
                                             /* 19/03/99 */
  int  i = 1, ok=0, e = cg->vl[n-1].e[e1-1]; 

  while (i <= MAXGLEN && !ok)
     if (i != n && cg->vl[i-1].vt != zv && 
          (cg->vl[i-1].e[0] == e || 
           cg->vl[i-1].e[1] == e || 
           cg->vl[i-1].e[2] == e))  ok = 1;
     else if (i == n)
           {  int j;
              for(j=1;j<3 && !ok;j++) if(cg->vl[i-1].e[(e1-1+j)%3] == e) ok=1;
              if (!ok)  i++;
            }
      else 
          i++; 
         
   if (i > MAXGLEN) cerror(253,"FindL: nonconnected edge"); 
   return i; 
} 


static cgraph * addcg(cgraph ** cg) 
{ /*- add C-graph CG to weight structure - 16/10/99 -*/   
   cgraph  *pgl = (cgraph *) malloc(SIZE);
   memcpy(pgl, *cg ,SIZE); 
   pgl->next = *cg; 
   *cg = pgl; 
   return pgl;
}

static void remqg_qg2(int n0,int n1,cgraph* cg)
{  int     n2; 
/* - remove subgraph (see figure) from C-graph - 08/01/90 */
     cg->powNN_1++;                        /*           v1     */
     cg->powN--;                           /*   v2 -->--*--    */
     cg->pow2--;                           /*           :  |   */
   n2 = findl(2,n1,cg);                    /*           :  |   */
   cg->vl[n2-1].e[2] = cg->vl[n0-1].e[2];  /*   v3 --<--*--    */
}                                          /*           v0     */


static void remqg_qg(int n0,int n1, cgraph ** c)
{ /* - remove gluon connected vertex n0 and n1 (see fugure) from first C-graph - 08/01/90  */ 
  int  n2, n5;    
  cgraph * cg=*c;   

#  if (CDEBLEV > DEBLEV) 
     printf(".......RemQG-QG........%u,%u\n",(unsigned int) n0,
             (unsigned int) n1);
#  endif 
   
   cg->vl[n0-1].vt = zv; 
   cg->vl[n1-1].vt = zv; 
   cg->gl -= 2; 
   if (cg->vl[n0-1].e[1] == cg->vl[n1-1].e[2] &&   /*         v1      */
       cg->vl[n0-1].e[2] == cg->vl[n1-1].e[1])     /*    -->--*--     */
   {                                               /*   |     :  |    */
      cg->powNN_1++;                               /*   |     :  |    */
      cg->pow2--;                                  /*    --<--*--     */
   }                                               /*         v0      */
   else  if (cg->vl[n0-1].e[1] == cg->vl[n1-1].e[2])  remqg_qg2(n0,n1,cg); 
   else  if (cg->vl[n0-1].e[2] == cg->vl[n1-1].e[1])  remqg_qg2(n1,n0,cg); 
   else 
   {  
      n2 = findl(2,n0,cg);                         /*         v0        */
      cg->vl[n2-1].e[2] = cg->vl[n1-1].e[2];       /*  v2-->--*-->--v3  */
      n5 = findl(2,n1,cg);                         /*         :         */
      cg->vl[n5-1].e[2] = cg->vl[n0-1].e[2];       /*         :         */
      cg->pow2--;                                  /*  v4--<--*--<--v5  */
                                                   /*         v1        */
      cg=addcg(c); 
      cg->sgn*=-1;
      cg->powN--; 
      cg->vl[n2-1].e[2] =  cg->vl[n0-1].e[2]; 
      cg->vl[n5-1].e[2] =  cg->vl[n1-1].e[2];  
   }  
#  if (CDEBLEV > DEBLEV) 
      wrcg(cg); 
      if (cg->next != NULL) 
         wrcg(cg->next); 
#  endif 
} 


static void rev3g(int en,cvertex* v)
{ /* - reverse 3G vertex such that edge EN will be first - 08/01/90  */  
   if (en != v->e[0]) 
   {
      if (en == v->e[1]) 
      { 
         v->e[1] = v->e[2]; 
         v->e[2] = v->e[0]; 
         v->e[0] = en; 
      } 
      else  if (en == v->e[2]) 
         { 
            v->e[2] = v->e[1]; 
            v->e[1] = v->e[0]; 
            v->e[0] = en; 
         }  
         else  cerror(255,"Rev3G: Invalid select edge");
   }      
}  


static void remqg_3g(int n0,int n1, cgraph  ** c)
/* - remove gluon connected vertex n0 and n1 (see figure)
     from first C-graph - 08/01/90  */ 
{int          n2, n3;   /*         v1        */ 
 int          en;       /*  v2.....*.....v3  */ 
 cgraph * cg = *c;      /*         :         */ 
                        /*         :         */ 
                        /*  v3-->--*-->--v4  */ 
                        /*          v0       */ 
#  if (CDEBLEV > DEBLEV) 
   fprintf(stderr,".......RemQG-3G........%u,%u\n",
      (unsigned int)n0,(unsigned int)n1);
#  endif 
   rev3g(cg->vl[n0-1].e[0],&cg->vl[n1-1]); 
   n2 = findl(2,n1,cg); 
   if (cg->vl[n2-1].vt == g3) 
      rev3g(cg->vl[n1-1].e[1],&cg->vl[n2-1]); 
   cg->vl[n0-1].vt = qg; 
   cg->vl[n0-1].e[0] = cg->vl[n2-1].e[0]; 
   en = cg->vl[n0-1].e[2]; 
   cg->vl[n0-1].e[2] = cg->vl[n1-1].e[0]; 
   n3 = findl(3,n1,cg); 
   if (cg->vl[n3-1].vt == g3) 
      rev3g(cg->vl[n1-1].e[2],&cg->vl[n3-1]); 
   cg->vl[n1-1].vt = qg; 
   cg->vl[n1-1].e[0] = cg->vl[n3-1].e[0]; 
   cg->vl[n1-1].e[1] = cg->vl[n0-1].e[2]; 
   cg->vl[n1-1].e[2] = en; 

   cg=addcg(c);
   cg->sgn*=-1;
   cg->vl[n0-1].e[0] = cg->next->vl[n1-1].e[0]; 
   cg->vl[n1-1].e[0] = cg->next->vl[n0-1].e[0]; 

#  if (CDEBLEV > DEBLEV) 
       wrcg(cg); 
       wrcg(cg->next); 
#  endif 
}  /* RemQG_3G */ 


static int istadpole(int n,cgraph* cg)
{ /* return True if vertex n is teadpole - 08/01/90  */ 
   return
      cg->vl[n-1].e[0] == cg->vl[n-1].e[1] || 
      cg->vl[n-1].e[1] == cg->vl[n-1].e[2] || 
      cg->vl[n-1].e[0] == cg->vl[n-1].e[2]    ; 
}


static void remg(int n0, cgraph ** c)
{ /* - remove gluon issue from vertex n0  from first C-graph - 08/01/90  */ 
   int  n1; 
#  if (CDEBLEV > DEBLEV) 
      printf(".......RemG........%u\n",(unsigned int)n0);
      wrcg(*c); 
#  endif 
   n1 = findl(1,n0,*c); 
   if (istadpole(n0,*c) || istadpole(n1,*c)) 
   { 
      (*c)->sgn=0;
      (*c)->gl = 0; 
      (*c)->vl[n0-1].vt = zv; 
      (*c)->vl[n1-1].vt = zv; 
#     if (CDEBLEV > DEBLEV) 
         wrcg(*c); 
#     endif 
   }
   else  if((*c)->vl[n1-1].vt==qg) remqg_qg(n0,n1,c); else remqg_3g(n0,n1,c); 

#  if (CDEBLEV > DEBLEV) 
     fprintf(stderr,".......end RemG........\n");
#  endif 
} /* RemG */ 


static void exp3g(int n0, cgraph ** c)
{ 
 int       n1, n2, n4, n5; 
 int       e04, e05, e45; 
 cgraph * cg=*c;
 
/* expand 3G vertex (see figure) - 08/01/90  */ 
/* 14/03/99: Check tadpole before expanding */    

#  if (CDEBLEV > DEBLEV) 
      fprintf(stderr,".......Exp3G........\n");
#  endif 

   n1 = findl(1,n0,cg); 
   if (istadpole(n0,cg) || istadpole(n1,cg)) 
   {
      cg->sgn=0; 
      cg->gl = 0; 
      cg->vl[n0-1].vt = zv; 
      cg->vl[n1-1].vt = zv; 
#     if (CDEBLEV > DEBLEV) 
         wrcg(cg); 
#     endif 
      return;
   }
   n1 = findl(1,n0,cg);      /*      v0            v4  v5     */   
   n2 = findl(2,n0,cg);      /*  v1..*...v2    v1..*-<-*..v2  */   
                             /*      :              \ /       */   
   n4 = getv(cg);            /*      :      ->       *v0      */
   e45 = ++(cg->en);         /*      :               :        */
   e04 = ++(cg->en);         /*      v3              v3       */
   cg->vl[n4-1].vt = qg; 
   if (cg->vl[n1-1].vt == g3) 
      rev3g(cg->vl[n0-1].e[0],&cg->vl[n1-1]); 
   cg->vl[n4-1].e[0] = cg->vl[n1-1].e[0]; 
   cg->vl[n4-1].e[1] = e45; 
   cg->vl[n4-1].e[2] = e04; 
   n5 = getv(cg); 
   e05 = ++(cg->en); 
   cg->vl[n5-1].vt = qg; 
   if (cg->vl[n2-1].vt == g3) 
      rev3g(cg->vl[n0-1].e[1],&cg->vl[n2-1]); 
   cg->vl[n5-1].e[0] = cg->vl[n2-1].e[0]; 
   cg->vl[n5-1].e[1] = e05; 
   cg->vl[n5-1].e[2] = e45; 
   rev3g(cg->vl[n0-1].e[2],&cg->vl[n0-1]); 
   cg->vl[n0-1].vt = qg; 
   cg->vl[n0-1].e[1] = e04; 
   cg->vl[n0-1].e[2] = e05; 
   cg->sgn*=-1;  cg->pow2++;

   cg=addcg(c);   /*  Second term  */ 
   cg->vl[n0-1].e[1] = e05; 
   cg->vl[n0-1].e[2] = e04; 
   cg->vl[n4-1].e[1] = e04; 
   cg->vl[n4-1].e[2] = e45; 
   cg->vl[n5-1].e[1] = e45; 
   cg->vl[n5-1].e[2] = e05; 
   cg->sgn*=-1;
   
#  if (CDEBLEV > DEBLEV) 
      wrcg(cg); 
      wrcg(cg->next); 
      fprintf(stderr,".......end Exp3G........\n");
#  endif 
}  /* Exp3G */ 


static void remtv(cgraph * pgl )
{ /* Remove transfered vertex from all C-graphs - 06/01/90     */ 

 int         n, n1; 
 int         vt0;   /* Original type */   
 int         ee;
                                           /*       v0        */ 
                                           /*  -->--*-->--v1  */ 
#if (CDEBLEV > DEBLEV)                     /*                 */ 
  fprintf(stderr,".......RemTV........\n");/*                 */
#endif                                     /*       v0        */ 
                                           /*  .....*.....v1  */ 

   while (findv(tv,pgl,&n) || findv(g2,pgl,&n))
   {
#  if (CDEBLEV > DEBLEV)                                   
         if (pgl != NULL) wrcg(pgl);                      
#  endif                                                   
      vt0 = pgl->vl[n-1].vt;                               
      pgl->vl[n-1].vt = zv;                                
      pgl->gl--;                                           
      if (istadpole(n,pgl))                                
         if (pgl->vl[n-1].e[0] != 0) pgl->powNN_1++;       
         else  pgl->powN++;                                
      else if (pgl->vl[n-1].e[0] != 0) {                   
        n1 = findl(1,n,pgl);                               
        if (pgl->vl[n1-1].vt == g2 &&                      
            pgl->vl[n1-1].e[0] != pgl->vl[n-1].e[0]) {     
          ee = pgl->vl[n1-1].e[0];                         
          pgl->vl[n1-1].e[0] = pgl->vl[n1-1].e[1];         
          pgl->vl[n1-1].e[1] = ee;                         
        }                                                  
        else if (pgl->vl[n1-1].vt == g3)                   
          rev3g(pgl->vl[n-1].e[0],&pgl->vl[n1-1]);         
        pgl->vl[n1-1].e[0] = pgl->vl[n-1].e[1];            
      }                                                    
      else {                                               
        n1 = findl(2,n,pgl);                               
        pgl->vl[n1-1].e[2] = pgl->vl[n-1].e[2];            
      }                                                    
   }                                                       
}   


factor *colorFactor (int nv, cvertex * vl)
{  /* - calculate color weight (two int n,d) - 08/01/90  */
   cgraph    * cg;
   int         n0,i,j;
   factor *f=(factor *)malloc(sizeof(factor));
   cgraph model;

   fct_init(f);

   SIZE= sizeof(model) + (nv)*((char*)&model.vl[1]-(char*)&model.vl[0]);
   MAXGLEN=nv+2;
      
   cg = (cgraph *) malloc(SIZE);   

   initcg(cg);
   cg->gl=nv;
   cg->next=NULL;
   cg->en=0;
   for(i=0;i<nv;i++)
   {  cg->vl[i]=vl[i];
      for(j=0;j<3;j++) if(vl[i].e[j]>cg->en) cg->en=vl[i].e[j];
   }

#  if (CDEBLEV > DEBLEV) 
       wrcg(cg); 
#  endif 
   remtv(cg); 
#  if (CDEBLEV > DEBLEV) 
       wrcg(cg); 
#  endif 
   while(cg)
   {  cgraph  * pgl;
      while (cg->gl && cg->sgn)
         if (findv(qg,cg,&n0)) remg(n0,&(cg));
         else if (findv(g3,cg,&n0)) exp3g(n0,&(cg));
         else cerror(251,"CWTarG: Invalid type of vertex.");

      pgl = cg;
      if(pgl->sgn)  add_fct(f,pgl->sgn,pgl->pow2,pgl->powN,pgl->powNN_1);
      cg = cg->next; 
      free(pgl); 
   }  
   return f;
} 


static void rednd(long * n,long * d,int b)
{ /* - reduce N and D with respect to B - 08/01/90  */
   if (b != 1)
      while (*d != 1 && *n % b == 0 && *d % b == 0)
      {
         *n /= b;
         *d /= b;
      }  /* while */
} /* RedND */


void fct_num_calc(factor * fct2,int Nc, long * n, long *d)
{ 
	int i;
	long p=1;
	factor *fct=fct2;
	
	*n=0;
	for(i=0;i<fct->len;i++) 
	{
		*n+=p*fct->nc[i]; 
		p*=Nc;
	}

	p=1;
	*d=fct->dc;
	if(fct->dpow>0)
	{
		for(i=0;i< fct->dpow;i++) *d*=Nc;
	}
	else
	{
		for(i=0;i<-fct->dpow;i++) *n*=Nc;
	}
			
	rednd(n,d,Nc); 
	rednd(n,d,2);
}

void fct_print(factor *fct, char *s)
{ int i;

  if(!fct->len)
  {
	  strcpy(s,"(0)"); 
	  return;
  }

  sprintf(s,"(");
  for(i=0;i<fct->len;i++)  
	  if(fct->nc[i]) 
	  {  
		  sprintf(s+strlen(s),"%+ld",fct->nc[i]);
		  if(i) sprintf(s+strlen(s),"N^%d",i);
	  }
  sprintf(s+strlen(s),")");

  if(fct->dpow>0) 
	  sprintf(s+strlen(s),"/(%ld*N^%d)",fct->dc,fct->dpow);
  else if(fct->dpow<0) 
	  sprintf(s+strlen(s),"*N^%d/%ld",-fct->dpow,fct->dc);
  else
	  sprintf(s+strlen(s),"/%ld",fct->dc);
}

  /* ************************** Cross reference ************************* */ 
  /* *                                                                  * */ 
  /* *  colorFactor                                                     * */ 
  /* *    +-----> RemTV                                                 * */ 
  /* *    |       +------> FindV                                        * */ 
  /* *    |       +------> isTadpole                                    * */ 
  /* *    |       +------> FindL                                        * */ 
  /* *    |       +------> Rev3G                                        * */ 
  /* *    |       +------> WrCG (Color)                                 * */ 
  /* *    |                                                             * */ 
  /* *    +-----> FindV                                                 * */ 
  /* *    +-----> RemG                                                  * */ 
  /* *    |       +------> WrCG (Color)                                 * */ 
  /* *    |       +------> FindL                                        * */ 
  /* *    |       +------> isTadpole                                    * */ 
  /* *    |       +------> RemQG_QG                                     * */ 
  /* *    |       |        +------> RemQG_QG1                           * */ 
  /* *    |       |        +------> RemQG_QG2                           * */ 
  /* *    |       |        |        +-------> FindL                     * */ 
  /* *    |       |        |                                            * */ 
  /* *    |       |        +------> FindL                               * */ 
  /* *    |       |        +------> AddCG                               * */ 
  /* *    |       |        +------> WrCG (Color)                        * */ 
  /* *    |       |                                                     * */ 
  /* *    |       +------> RemQG_3G                                     * */ 
  /* *    |                +------> Rev3G                               * */ 
  /* *    |                +------> FindL                               * */ 
  /* *    |                +------> AddCG                               * */ 
  /* *    |                +------> WrCG                                * */ 
  /* *    |                                                             * */ 
  /* *    +-----> Exp3G                                                 * */ 
  /* *    |       +------> FindL                                        * */ 
  /* *    |       +------> GetV (Color)                                 * */ 
  /* *    |       +------> GetEN (Color)                                * */ 
  /* *    |       +------> Rev3G                                        * */ 
  /* *    |       +------> AddCG                                        * */ 
  /* *    |                                                             * */ 
  /* *    +-----> CError (Color)                                        * */ 
  /* *    +-----> DispCG                                                * */ 
  /* *            +------> RedND                                        * */ 
  /* *                                                                  * */ 
  /* ******************************************************************** */

vtype typev(vert0 v,int valence)
{int  ng = 0, ne, nq = 0;
 
/* Return color type of vertex - 06/01/90  */ 
      
   for (ne = 0; ne < valence; ne++)
      if (prtclbase[v[ne].partcl-1].cdim != 1) 
      {   if (prtclbase[v[ne].partcl-1].cdim == 8) ng++; else nq++;} 
   switch (ng) 
   {    
      case 0:   return nq == 2 ? tv : zv; 
      case 1:   return qg; 
      case 2:   return g2; 
      case 3:   return g3; 
      default:  return cerror(252,"TypeV: invalid vertex type"); 
   }  /* case */ 
}  /* TypeV */ 


void t2k2(vcsect* g, int * nv, cvertex * vl)
{
  int   i,k, ne, en=0, maxnv=*nv; 
  int   maptar[2*maxvert][MAXVALENCE];
    
/* Transfer Taranov's representation of graph to Kryukov's representation  */ 
/* - 07/01/90  */ 
    
   for(i=0;i<2*maxvert;i++) for(k=0;k<MAXVALENCE;k++) maptar[i][k]=0;
   *nv=0;
   for(i=0;i<g->sizet; i++)  if (typev(g->vertlist[i],g->valence[i]) != zv) 
   { 
      if(*nv>= maxnv) cerror(251,"To many vertices");
      vl[*nv].e[0]=0; vl[*nv].e[1]=0;  vl[*nv].e[2]=0;
      vl[*nv].vt = typev(g->vertlist[i],g->valence[i]); 

      for (k=0, ne=0; ne<g->valence[i] ;ne++)
      {  int dim= prtclbase[g->vertlist[i][ne].partcl-1].cdim;
         if(dim!=1)
         {  int l=maptar[i][ne];
            if(!l)
            { 
               l = ++en; 
               maptar[i][ne] = l;
               maptar[g->vertlist[i][ne].link.vno]
                     [g->vertlist[i][ne].link.edno]=l; 
            } 
         
            if (vl[*nv].vt != g3 && vl[*nv].vt != g2) 
            switch (dim) 
            {
              case  8:  vl[*nv].e[0]=l;  break; 
              case -3:  vl[*nv].e[1]=l;  break; 
              case  3:  vl[*nv].e[2]=l;  break; 
              default: cerror(252,"t2k - invalid particle color");
            } else vl[*nv].e[k++] = l;
         } 
      }  
      (*nv)++;
   } 
}

