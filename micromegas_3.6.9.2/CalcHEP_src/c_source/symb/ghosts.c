/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include "physics.h"
#include "syst.h"
#include "ghosts.h"
#include "process.h"

#define S_PERMIT 8
#define V_PERMIT 16
#define X_PERMIT 32
#define Y_PERMIT 64

j_lgrptr j_lgrarray[2 * maxvert];


#define maxmark 64
#define undefprtclmark 8

static int         sklt[2 * maxvert][MAXVALENCE];
static int         vertmrk[2 * maxvert];
static vcsect      vcs_stat;

void  vertinlgr(int * vert_,int nvert,int * subst,algvertptr* lgr)
{j_lgrptr     q;
 int         i, m;
 arr4byte     s;
 arr4byte vert;

   for(i=0;i<4;i++) vert[i]= vert_[i] ? prtclbase[vert_[i]-1].anti : 0;

   q = j_lgrarray[nvert-1];
   for (i = 1; i <= 4; i++) s[i-1] = i;
   i = 1;
   while (i < 4)
   {
      if (vert[i-1] >= vert[i]) ++(i);
      else
      {
         m = vert[i-1];
         vert[i-1] = vert[i];
         vert[i] = m;
         m = s[i-1];
         s[i-1] = s[i];
         s[i] = m;
         if (i == 1) ++(i);
         else --(i);
      }
   }

   for(q=j_lgrarray[nvert-1]; q != NULL &&  
                           (q->lgrnptr->fields[0] != vert[0] || 
                            q->lgrnptr->fields[1] != vert[1] || 
                            q->lgrnptr->fields[2] != vert[2] || 
                            q->lgrnptr->fields[3] != vert[3]);  q = q->next);

   if(q)
   {
      *lgr = q->lgrnptr;
      for (i = 0; i < 4; i++) subst[(*lgr)->perm[i]-1] = s[i];
   } else *lgr=NULL;
}



static void  fillLgrArray(void)
{int  i, j, m, mm; 
 int  f[4]; 
 algvertptr   lagr; 
 j_lgrptr     j_lgr; 

   for (i = 0; i < 2 * maxvert; i++) j_lgrarray[i] = NULL; 
   for (i = 0; i < vcs_stat.sizet; i++) 
   { 
      for (j = 0; j < vcs_stat.valence[i]; j++) 
         f[j] = prtclbase[vcs_stat.vertlist[i][j].partcl-1].anti; 
      m = 1; 
        /* Sorting */ 
      while (m < vcs_stat.valence[i]) 
         if (f[m-1] >= f[m])  ++(m); 
         else 
         { 
            mm = f[m-1]; 
            f[m-1] = f[m]; 
            f[m] = mm; 
            if (m == 1) ++(m); else --(m); 
         } 
      if (vcs_stat.valence[i] == 3) f[3] = 0; 
      lagr = lgrgn; 
      do 
      { 
         if (ghostmother(lagr->fields[0]) == f[0] && 
             ghostmother(lagr->fields[1]) == f[1] && 
             ghostmother(lagr->fields[2]) == f[2] && 
             ghostmother(lagr->fields[3]) == f[3]) 
         { 
            j_lgr = (j_lgrptr) m_alloc(sizeof(struct j_lgrrec));
            j_lgr->next = j_lgrarray[i]; 
            j_lgr->lgrnptr = lagr; 
            j_lgrarray[i] = j_lgr; 
         } 
         lagr = lagr->next; 
      }  while (lagr != NULL);
   } 
} 



static void     addLine(int i, int j, int n);
 
static void  addVertex(int i,int  n)
{int  j; 
   vertmrk[i-1] = n; 
   for (j = 1; j <= vcs_stat.valence[i-1]; j++)  
    if (sklt[i-1][j-1] == 1) addLine(i,j,n);
}


static void  addLine(int i,int  j,int  n)
{int   i1, j1;

    i1 = vcs_stat.vertlist[i-1][j-1].nextvert.vno;
    j1 = vcs_stat.vertlist[i-1][j-1].nextvert.edno;
    if (vertmrk[i1-1] == 1)
    {
       sklt[i1-1][j1-1] = n;
       sklt[i-1][j-1] = n;
       addVertex(i1,n - 1);
    }
 }


static int  prtclInVert(int  p,int  v)
{j_lgrptr  q; 
 int      i; 

    p=prtclbase[p-1].anti; 
    q = j_lgrarray[v-1]; 
    while (q != NULL) 
    { 
       for (i = 0; i < 4; i++) 
          if (q->lgrnptr->fields[i] == p)  return 1; 
       q = q->next; 
    } 
    return 0; 
 } 


static void  makeSkelet(void)
{int   i, j, p, i1; 
 
   for (i = 0; i < vcs_stat.sizet; i++) 
   { 
      vertmrk[i] = 0; 
      for (j = 0; j < vcs_stat.valence[i]; j++) sklt[i][j] = 0; 
   } 

   for (i = 0; i < vcs_stat.sizet; i++)   /*  mark subgraph vertex with '1'  */ 
      for (j = 0; j < vcs_stat.valence[i]; j++) 
      { 
/*       if(PLR_PRTCL&vcs_stat.vertlist[i][j].prop) continue;*/
         i1 = vcs_stat.vertlist[i][j].nextvert.vno - 1 ; 
         p = vcs_stat.vertlist[i][j].partcl; 
         if (gaugep(p) && i < i1) 
         {   /*  check the  Ghost particles existence in vertex  */ 
            if ((prtclInVert(p + 1,i+1) && 
                 prtclInVert(prtclbase[p + 1-1].anti,i1+1)) || 
                (prtclInVert(p + 2,i+1) && 
                 prtclInVert(prtclbase[p + 2-1].anti,i1+1))) 
            { 
               sklt[i][j] = 1; 
               sklt[i1][vcs_stat.vertlist[i][j].nextvert.edno-1] = 1; 
               vertmrk[i] = 1; 
               vertmrk[i1] = 1; 
            } 
         } 
      } 
   /*  end of mark   */ 
   for (i = 0; i < vcs_stat.sizet; i++) {if (vertmrk[i] == 1) addVertex(i+1,maxmark);}
}  


static void  skipup(int* v,hlpcsptr lpl)
{int   j, i1, j1; 

   j = 1; 
   while (sklt[*v-1][j-1] != vertmrk[*v-1] + 1) ++(j); 
   i1 = vcs_stat.vertlist[*v-1][j-1].nextvert.vno; 
   j1 = vcs_stat.vertlist[*v-1][j-1].nextvert.edno; 
   lpl->hlpcs[*v-1][j-1] = 1; 
   lpl->hlpcs[i1-1][j1-1] = 1; 
   *v = i1; 
} 


static void  mkindeploops(hlpcsptr* indpl)
{int        i, j, ii, jj, i1, j1, v1, v2; 
 hlpcsptr    lpl; 

   *indpl = NULL; 
   for (i = 1; i <= vcs_stat.sizet; i++) 
   for (j = 1; j <= vcs_stat.valence[i-1]; j++) 
   if (sklt[i-1][j-1] == 1)   /*   if (i,j) is't element of OCTOV .... */ 
   { 
      i1 = vcs_stat.vertlist[i-1][j-1].nextvert.vno; 
      j1 = vcs_stat.vertlist[i-1][j-1].nextvert.edno; 
      if (i < i1) 
      { 
         lpl = (hlpcsptr) m_alloc(sizeof(struct hlpcsrec));
         lpl->next = *indpl; 
         *indpl = lpl; 

         for (ii = 1; ii <= vcs_stat.sizet; ii++) 
             for (jj = 1; jj <= vcs_stat.valence[ii-1]; jj++) 
                lpl->hlpcs[ii-1][jj-1] = 0; 

         lpl->hlpcs[i-1][j-1] = 1; 
         lpl->hlpcs[i1-1][j1-1] = 1; 

         v1 = i; v2 = i1; 

         while (vertmrk[v1-1] > vertmrk[v2-1]) skipup(&v2,lpl); 
         while (vertmrk[v2-1] > vertmrk[v1-1]) skipup(&v1,lpl); 
         while (v1 != v2) 
         {  skipup(&v1,lpl); 
            skipup(&v2,lpl); 
         } 
      } 
   }
   
    
}  


static void  mkallloops(hlpcsptr indpl, hlpcsptr* alll)
{hlpcsptr   alll_, tmpl; 
 int       i, j; 

   if (indpl == NULL) 
      *alll = NULL; 
   else 
   { 
      mkallloops(indpl->next,&alll_); 
      *alll = indpl; 
      (*alll)->next = alll_; 
      while (alll_ != NULL) 
      { 
         tmpl = (hlpcsptr) m_alloc(sizeof(struct hlpcsrec));
         tmpl->next = *alll; 
         *alll = tmpl; 
         for (i = 1; i <= vcs_stat.sizet; i++) 
            for (j = 1; j <= vcs_stat.valence[i-1]; j++) 
               (*alll)->hlpcs[i-1][j-1] = 
                  (indpl->hlpcs[i-1][j-1] + alll_->hlpcs[i-1][j-1]) % 2; 
         alll_ = alll_->next; 
      } 
   }   
    
} 


static int  find1(int* i,int* j,hlpcsptr tmp)
{ 
   while(1)
   { 
      if (tmp->hlpcs[*i-1][*j-1] == 1) 
         return 1; 
      else 
         if (++(*j) > vcs_stat.valence[*i-1]) 
         { 
            if (++(*i) > vcs_stat.sizet) return 0; 
            *j = 1; 
         } 
   }
}  


static void  insertorient(int i1,int j1,hlpcsptr tmp)
{int      i2, j2; 
 hlpcsptr  tmpnext; 

   tmp->sgn = -tmp->sgn; 
   tmpnext = (hlpcsptr) m_alloc(sizeof(struct hlpcsrec));
   *tmpnext = *tmp; 
   tmp->next = tmpnext; 

   while (j1 <= vcs_stat.valence[i1-1]) 
   { 
      i2 = vcs_stat.vertlist[i1-1][j1-1].nextvert.vno; 
      j2 = vcs_stat.vertlist[i1-1][j1-1].nextvert.edno; 
      tmp->hlpcs[i1-1][j1-1] = 2; 
      tmp->hlpcs[i2-1][j2-1] = 3; 
      tmpnext->hlpcs[i1-1][j1-1] = 3; 
      tmpnext->hlpcs[i2-1][j2-1] = 2; 
      i1 = i2; 
      j1 = 1; 
      while (tmpnext->hlpcs[i1-1][j1-1] != 1 && 
             j1 <= vcs_stat.valence[i1-1]) ++(j1); 
   } 
}


static void  mkorientedloops(hlpcsptr* alll)
{int        i, j; 
 hlpcsptr    tmp; 

   tmp = *alll; 
   while (tmp != NULL) 
   { 
      tmp->sgn = 1; 
      tmp = tmp->next; 
   } 
   tmp = *alll; 
   while (tmp != NULL) 
   { 
      i = 1; j = 1; 
      while (find1(&i,&j,tmp)) insertorient(i,j,tmp); 
      tmp = tmp->next; 
   } 
   tmp = *alll; 
   while (tmp != NULL) 
   { 
      for (i = 0; i < vcs_stat.sizet; i++) 
      for (j = 0; j < vcs_stat.valence[i]; j++) 
      if (tmp->hlpcs[i][j] != 0) tmp->hlpcs[i][j] --  ;               
      tmp = tmp->next; 
   } 

   tmp = (hlpcsptr) m_alloc(sizeof(struct hlpcsrec)); 
   /*    addition  of origin graph  */ 
   for (i = 1; i <= vcs_stat.sizet; i++) 
      for (j = 1; j <= vcs_stat.valence[i-1]; j++) tmp->hlpcs[i-1][j-1] = 0; 
   tmp->sgn = 1; 
   tmp->next = *alll; 
   *alll = tmp;
          
}   /*  MkOrientedLoops  */ 


static void  insertPermition(hlpcsptr alll)
{ 
 int          i, j, i1, j1; 
 int          first_mark[2 * maxvert][MAXVALENCE];   
 int mark, np;
 
    for (i = 0; i < vcs_stat.sizet; i++) 
    for (j = 0; j < vcs_stat.valence[i]; j++)
    if(vcs_stat.vertlist[i][j].nextvert.vno != nullvert) 
    { 
       i1 = vcs_stat.vertlist[i][j].nextvert.vno-1;
       if (i < i1) 
       { 
          j1 = vcs_stat.vertlist[i][j].nextvert.edno-1;
          np= vcs_stat.vertlist[i][j].partcl;
          mark=0;
          if (gaugep(np) && !zeromass(np) && prtclInVert(np+sbosonmark,i+1) 
          && prtclInVert(prtclbase[np+sbosonmark-1].anti,i1+1) )mark +=S_PERMIT;

          if ( (prtclbase[np-1].spin == 2) /*&& (prtclbase[np-1].cdim !=1)*/ &&
          ( (i >= vcs_stat.sizel) || (i1 < vcs_stat.sizel) ))
          {  if(prtclInVert(prtclbase[np+xbosonmark-1].anti,i1+1)
             && prtclInVert(np+xbosonmark,i+1)) mark += X_PERMIT; 

             if(prtclInVert(prtclbase[np+ybosonmark-1].anti,i1+1 )
             && prtclInVert(np+ybosonmark,i+1)) mark += Y_PERMIT; 
          }
          first_mark[i ][j ]=mark; 
          first_mark[i1][j1]=mark;   
       }    
    } else  first_mark[i ][j ]=0;
          
   while (alll != NULL) 
   { 
      for (i = 0; i < vcs_stat.sizet; i++) 
      for (j = 0; j < vcs_stat.valence[i]; j++) 
      if (alll->hlpcs[i][j] == 0 )  alll->hlpcs[i][j] = first_mark[i][j];   
      alll = alll->next; 
   } 
} 


static void  preliminaryTest(hlpcsptr* alll)
{hlpcsptr c, cpred;
 arr4byte    vert, subst; 
 algvertptr  lgr; 
 int        i, j; 
 int     del; 


   c = *alll; 
   while (c != NULL) 
   { 
      del = 0; 
      
      for(i=0; i < vcs_stat.sizet && !del;i++)
      { 
            for (j = 0; j < vcs_stat.valence[i]; j++) 
               if (c->hlpcs[i][j] >= undefprtclmark)  goto label_1;
            for (j = 0; j < vcs_stat.valence[i]; j++) 
               vert[j] = vcs_stat.vertlist[i][j].partcl + c->hlpcs[i][j]; 
            for (j=vcs_stat.valence[i];j<MAXVALENCE;j++) vert[j] = 0;
            
            vertinlgr(vert,i+1,subst,&lgr); 
            del = (lgr == NULL);
label_1: ;
      }
      if (del) 
      { 
         if (c == *alll) 
         { 
            c = c->next; 
            free(*alll); 
            *alll = c; 
         } 
         else 
         { 
            cpred->next = c->next; 
            free(c); 
            c = cpred->next; 
         } 
      } 
      else 
      { 
         cpred = c; 
         c = c->next; 
      } 
   } 
} 


static void  insertcopy(hlpcsptr c)
{  hlpcsptr cc; 
 
   cc = (hlpcsptr) m_alloc(sizeof(struct hlpcsrec));
   *cc = *c; 
   c->next = cc; 
} 


static int  ins_test(int i,int  j,int  mrk,hlpcsptr  c)
{int      j1; 
 algvertptr   lgr; 
 arr4byte     subst, vert; 
 int         m;
    
   for (j1 = 1; j1 <= vcs_stat.valence[i-1]; j1++) 
   { 
      if (j1 == j) m = mrk; else 
      {   m = c->hlpcs[i-1][j1-1]; 
        if (m >=undefprtclmark) return 1;
      }
      vert[j1-1] = vcs_stat.vertlist[i-1][j1-1].partcl +m;         
   } 
   if (vcs_stat.valence[i-1] == 3) vert[3] = 0; 
   vertinlgr(vert,i,subst,&lgr); 
   return (lgr != NULL); 
} 

static void  insert_v_s_t(hlpcsptr* alll)
{  hlpcsptr    c, cpred; 
   int        i, i1, j, j1; 
   int     del; 
   char  markList[4];
   int lineMark ;
   int l,k;
   
   c = *alll; 
   while (c != NULL) 
   { 
      del = 0; 
      for (i = 1; i <= vcs_stat.sizet; i++) 
         for (j = 1; j <= vcs_stat.valence[i-1]; j++) 
         { 
            i1 = vcs_stat.vertlist[i-1][j-1].nextvert.vno; 
            if (i < i1) 
            {  
               lineMark=c->hlpcs[i-1][j-1];
               if ( lineMark >= undefprtclmark) 
               {  
                  j1 = vcs_stat.vertlist[i-1][j-1].nextvert.edno;                 
                  l=0;
                  if (  ins_test(i,j,vbosonmark,c) 
                      && ins_test(i1,j1,vbosonmark,c)  ) 
                  { markList[l] = vbosonmark; l++;}                   
                  if ( ( S_PERMIT & lineMark ) &&  ins_test(i,j,sbosonmark,c) 
                      && ins_test(i1,j1,sbosonmark,c)  ) 
                  { markList[l] = sbosonmark; l++;}
                  if ( ( X_PERMIT & lineMark ) && ins_test(i,j,xbosonmark,c) 
                      && ins_test(i1,j1,xbosonmark,c)  ) 
                  { markList[l] = xbosonmark; l++;}
                  if ( ( Y_PERMIT & lineMark ) && ins_test(i,j,ybosonmark,c) 
                      && ins_test(i1,j1,ybosonmark,c)  ) 
                  { markList[l] = ybosonmark; l++;}
                  
                  if (l==0 ) { del=1; goto label_1;} else
                  {  for (k=0;k<=l-2;k++)
                     { insertcopy(c);
                        c->next->hlpcs[i -1][j -1] = markList[k];
                        c->next->hlpcs[i1-1][j1-1] = markList[k];  
                     }
                     c->hlpcs[i -1][j -1] = markList[l-1];
                     c->hlpcs[i1-1][j1-1] = markList[l-1];                                          
                  }                   
               } 
            } 
         } 
label_1: 
      if (del) 
      { 
         if (c == *alll) 
         { 
            c = c->next; 
            free(*alll); 
            *alll = c; 
         } 
         else 
         { 
            cpred->next = c->next; 
            free(c); 
            c = cpred->next; 
         } 
      } 
      else 
      { 
         cpred = c; 
         c = c->next; 
      } 
   } 
} 


static void  numerate(hlpcsptr alll)
{unsigned        nn; 
 hlpcsptr    tmp; 

   nn = 0; 
   tmp = alll; 
   while (tmp != NULL) 
   { 
      ++(nn); 
      tmp->num = nn; 
      tmp = tmp->next; 
   } 
   tmp = alll; 
   while (tmp != NULL) 
   { 
      tmp->maxnum = nn; 
      tmp = tmp->next; 
   } 
} 


static void  setLorInd(hlpcsptr alll)
{  
  int i,j,i1,j1,np, ind1, ind2, lorCount;
  hlpcsptr    tmp; 
  lorCount=0;
  for (i=0;i<vcs_stat.sizet;i++)
  for (j=0;j<vcs_stat.valence[i];j++)
  { 
    i1=vcs_stat.vertlist[i][j].link.vno;
    if (i < i1)
    {   
       j1=vcs_stat.vertlist[i][j].link.edno;
       np=vcs_stat.vertlist[i][j].partcl;
       switch(prtclbase1[np].spin)
       { 
         case 0:;
         case 1: ind1=0;ind2=0; break;
         case 2:       
              lorCount++;
              if(PLR_PRTCL & vcs.vertlist[i][j].prop)
              {  
                ind1= lorCount++;
                ind2= lorCount;
                break;
              }      
              if (/* ( prtclbase1[np].cdim != 1 ) && */
                 ( (i >= vcs_stat.sizel) || (i1 < vcs_stat.sizel) ) )
               for(tmp = alll; tmp; tmp = tmp->next) 
               if(tmp->hlpcs[i][j]==xbosonmark||tmp->hlpcs[i][j]==ybosonmark) 
               { lorCount ++; break; }  
              ind1 = lorCount; 
              ind2 = lorCount;          
              break;
         case 3: ind1= ++lorCount;
                 ind2= ++lorCount;
                 break;
         case 4: 
              lorCount+=2;
              ind1=lorCount;
              lorCount+=2;
              ind2=lorCount;       
                 break;
       }
       vcs_stat.vertlist[i][j].lorentz=ind1;
       vcs_stat.vertlist[i1][j1].lorentz=ind2;
    }    
  }     
}

static void  vectorFactor(hlpcsptr alll)
{  
  int i,j,i1;
  hlpcsptr    tmp;
  
  for(tmp = alll; tmp; tmp = tmp->next)
  {    
    for (i=0;i<vcs_stat.sizet;i++)
    for (j=0;j<vcs_stat.valence[i];j++)
    { 
      i1=vcs_stat.vertlist[i][j].nextvert.vno -1;
      if (i < i1)
      { int  np=vcs_stat.vertlist[i][j].partcl-1;
        if(prtclbase[np].spin==2  && tmp->hlpcs[i][j]==vbosonmark) tmp->sgn*=-1;
      }
    }
  }   
}

   

void  generateghosts(vcsect* vcs,hlpcsptr*  alll)
{  hlpcsptr   indpl; 

   vcs_stat=*vcs;
   fillLgrArray();   /*  for every diagram vertex a list       */ 
                     /*  of lagrangian vertex is attached     */ 

   makeSkelet();           /*  this is a set of program       */ 
   mkindeploops(&indpl);   /*  to create list of all posible  */ 
   mkallloops(indpl,alll); /*  ghost and  antiGhost  loops    */ 
   mkorientedloops(alll);  /*                                 */ 
   insertPermition(*alll); 
   preliminaryTest(alll);  /*  remove Ghost maps that have    */ 
   insert_v_s_t(alll); 
   numerate(*alll);
   setLorInd(*alll);
   vectorFactor(*alll);
   *vcs=vcs_stat;
}  /*  GenerateGousts  */  


void  eraseghosts(hlpcsptr gst)
{hlpcsptr     gst_; 
 int         i; 
 j_lgrptr     j_list, j_listmem; 

   while (gst != NULL) 
   { 
      gst_ = gst; 
      gst = gst_->next; 
      free(gst_); 
   } 
   for (i = 0; i < 2 * maxvert; i++) 
   { 
      j_list = j_lgrarray[i]; 
      while (j_list != NULL) 
      { 
         j_listmem = j_list; 
         j_list = j_list->next; 
         free(j_listmem); 
      } 
   } 
} 

void  ghostsForAmplitude(vampl * va, hlpcsptr * ghosts)
{  int i,j;

   vcs_stat.sizel = va->size;
   vcs_stat.sizet = va->size;

   for (i = 0; i < vcs_stat.sizel; i++)
   { vcs_stat.valence[i]=va->valence[i];
     for(j=0;j<MAXVALENCE;j++) vcs_stat.vertlist[i][j] = va->vertlist[i][j];
   }

   *ghosts=(hlpcsptr) m_alloc(sizeof(**ghosts));
   (*ghosts)->next=NULL;  
   for(i=0;i<2*maxvert;i++) for(j=0;j<MAXVALENCE;j++) (*ghosts)->hlpcs[i][j]=0;
   
 
   fillLgrArray();

   insertPermition(*ghosts);

   preliminaryTest(ghosts);  /*  remove Ghost maps that have    */
   insert_v_s_t(ghosts);
}


