/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include"chep_crt.h"
#include "syst2.h"
#include "process.h"
#include "physics.h"
#include "s_files.h"
#include "squar.h"
#include "sortarr.h"

typedef struct groupStr
{
   struct groupStr *  next;
   permut       perm;
   unsigned     left;
   unsigned     res;
}  groupStr;
   
typedef struct groupStr *grouplist;


static void  smpl_linker(permut perm1,permut perm2,permut perm)
{  int  i, j, c;
   permut perma,permb;

   memcpy(perma,perm1,sizeof(perma));
   memcpy(permb,perm2,sizeof(permb));
   for (i = 1; i <= nout; i++)
   {
      c = perma[i-1];
      j = 1;
      while (c != permb[j-1]) j++;
      perm[i-1] = j;
      permb[j-1] = 0;
   }
}


static permut      allptest; 

static int     allpout[MAXINOUT]; 

static void  constrsim(int i,permut stab,grouplist* sim)
{int        j; 
 grouplist    ad_sim; 

   if (i > nout) 
   { 
      ad_sim = (grouplist)m_alloc(sizeof(struct groupStr));
      memcpy(ad_sim->perm,allptest,sizeof(ad_sim->perm)); 
      ad_sim->next = *sim; 
      *sim = ad_sim; 
   } 
   else 
   { 
      for (j = nout; j >= 1; j--) 
      if (allpout[j-1] && (stab[i-1] == stab[j-1])) 
      { 
         allpout[j-1] = 0; 
         allptest[i-1] = j; 
         constrsim(i + 1,stab,sim); 
         allpout[j-1] = 1; 
      } 
   } 
} 

static void  allperm(permut stab,grouplist* sim)
{ int        i; 
  *sim = NULL; 
  for (i = 0; i < nout; i++) allpout[i] = 1; 
  constrsim(1,stab,sim); 
} 


static void  multperm(permut p1,permut p2,permut pres)
{int i=0; permut r;  
   for(; i < nout; i++)  r[i] = p2[p1[i]-1];
   for(i=0; i<nout;i++)  pres[i]=r[i];  
} 

static void  revers1(permut p,permut pres)
{int  i=0; for (; i < nout; i++) pres[p[i]-1] = i+1; } 


static int  eq(decayDiagram branch1,decayDiagram branch2,int l1,int l2)
{int  i; 

   if (l1 == l2) for(i=0; i<l1; i++) if(branch1[i] != branch2[i]) return 0; 
   return (l1 == l2); 
} 

static void  addgen(int nbeg1,int nbeg2,int nl,permlist* gen)
{permut       prm; 
 int         i; 
 permlist     ad_gen; 

   for (i = 1; i <= nout; i++) prm[i-1] = i; 
   for (i = 0; i <= nl - 1; i++) 
   { 
      prm[nbeg1 + i-1] = nbeg2 + i; 
      prm[nbeg2 + i-1] = nbeg1 + i; 
   } 
   ad_gen = (permlist)m_alloc(sizeof(struct permListStr));
   memcpy(ad_gen->perm,prm,sizeof(ad_gen->perm)); 
   ad_gen->next = *gen; 
   *gen = ad_gen; 
} 


static void  onebranch(decayDiagram branch,int beg,int * l,decayDiagram grph)
{int  c, g, k; 

   c = 0; 
   *l = 0; 
   k = 1; 
   while (c != 1) 
   { 
      g = grph[k + beg - 1-1]; 
      branch[k-1] = g; 
      (*l)++; 
      k++; 
      if (g > 0) c++; else c--; 
   } 
} 

static void  dobranches(int root,decayDiagram grph,decayDiagram branch1,
                 decayDiagram branch2,decayDiagram branch3,int * l1,int * l2,int* l3)
{ 
   onebranch(branch1,root + 1,l1,grph); 
   if (grph[root + *l1 + 1-1] != 0) 
   { 
      onebranch(branch2,root + *l1 + 1,l2,grph); 
      *l3 = 0; 
   } 
   else 
   { 
      onebranch(branch2,root + *l1 + 2,l2,grph); 
      onebranch(branch3,root + *l1 + *l2 + 2,l3,grph); 
   } 
} 

static void  dogen(decayDiagram grph,permlist* gen,unsigned* dim)
{int     nn; 
 int     i, ngen; 
 decayDiagram branch1, branch2, branch3; 
 int     nbeg1, nbeg2, nbeg3, l1, l2, l3, nl1, nl2; 

   nn = 1 - nin; 
   *gen = NULL; 
   *dim = 1; 
   for (i = 1; i <= 2 * (nin + nout) - 5; i++) 
   { 
      if (grph[i-1] > 0) nn++; 
      if (grph[i-1] < 0) 
      { 
         dobranches(i,grph,branch1,branch2,branch3,&l1,&l2,&l3); 
         nbeg1 = nn + 1; 
         nl1 = 1 + (l1 - 1) / 2; 
         nbeg2 = nbeg1 + nl1; 
         nl2 = 1 + (l2 - 1) / 2; 
         nbeg3 = nbeg2 + nl2; 
         ngen = 0; 
         if (eq(branch2,branch3,l2,l3)) 
         {  addgen(nbeg2,nbeg3,nl2,gen); 
            ngen++; 
         } 
         if (nn >= 0) 
         { 
            if (eq(branch1,branch2,l1,l2)) 
            {  addgen(nbeg1,nbeg2,nl1,gen); 
               ngen++; 
            } 
            if (eq(branch1,branch3,l1,l3) && ngen < 2) 
            {  addgen(nbeg1,nbeg3,nl1,gen); 
               ngen++; 
            } 
         } 
         if (ngen == 1) *dim *= 2; 
         if (ngen == 2) *dim *= 6; 
      } 
   } 
} 


static int  eqprm(permut perm1,permut perm2)
{ int    i; 

   for (i = 1; i <= nout; i++) 
   if (perm1[i-1] != perm2[i-1]) return 0; 
   return 1; 
} 

static void  doleftident(grouplist group,  ampllist  amp1)
{  grouplist   sim, sim1; 
   permlist    lgen; 
   permut      res; 
   int         c, c1, c2, n; 

   sim = group; 
   n = 1; 
   while (sim != NULL) 
   {  sim->left = n++; 
      sim = sim->next; 
   }  
   lgen = amp1->gen; 
   while (lgen != NULL) 
   { 
      sim = group; 
      while (sim != NULL) 
      { 
         multperm(lgen->perm,sim->perm,res); 
         sim1 = sim; 
         do 
         { 
            if (eqprm(res,sim1->perm)) 
            { 
               c1 = sim1->left; 
               c = sim->left; 
               if (c1 != c) 
               { 
                  if (c1 > c) 
                     c2 = c; 
                  else 
                  {  c2 = c1; 
                     c1 = c; 
                  } 
                  sim1 = group; 
                  while(sim1 != NULL)
                  { 
                     if (sim1->left == c1) sim1->left = c2; 
                     sim1 = sim1->next; 
                  }
               } 
               goto label_1;
            } 
            sim1 = sim1->next; 
         }  while (sim1 != NULL); 
label_1: sim = sim->next; 
      } 
      lgen = lgen->next; 
   } 
} 


static void  doresident(grouplist group, ampllist amp2)
{grouplist   sim, sim1; 
 permlist    rgen; 
 permut      res; 
 unsigned        c,c1,n; 

   sim = group; 
   while (sim != NULL) 
   { 
      sim->res = sim->left; 
      sim = sim->next; 
   } 
   rgen = amp2->gen; 
   while (rgen != NULL) 
   { 
      n = 1; 
      sim = group; 
      while (sim != NULL) 
      {  if (sim->left == n) 
         { 
            multperm(sim->perm,rgen->perm,res); 
            sim1 = sim; 
            do 
            {  if (eqprm(res,sim1->perm)) 
               { 
                  c1 = sim1->res; c = sim->res; 
                  if (c1 != c) 
                  { 
                     sim1 = sim; 
                     while (sim1 != NULL) 
                     { 
                        if (sim1->res == c1) sim1->res = c; 
                        sim1 = sim1->next; 
                     } 
                  } 
                  goto label_1;
               } 
               sim1 = sim1->next; 
            }  while (sim1 != NULL);
         } 
label_1: sim = sim->next; 
         n++; 
      } 
      rgen = rgen->next; 
   } 
} 



static void  clearGroop(grouplist  g)
{ while (g) {   grouplist g_=g; g=g->next; free(g_); }  } 


void clearDiagrams(ampllist a)
{
 
   ampllist a_;
   permlist p,p_;
   while (a) 
   {  
      p=a->gen;
      while (p) 
      { 
         p_ = p; 
         p=p->next; 
         free(p_); 
      }
      a_=a;
      a=a->next; 
      free(a_); 
   } 
}

static ampllist readDiagrams(FILE * diagrp, permut cononout)
{
   int i;
   ampllist  amplitudes,amp1,amp2;
   adiagram buff;
   permut       rvrs;
   permlist     cgen;
   int nsub0;
   long fpos;

   amplitudes = (ampllist)m_alloc(sizeof(struct amplListStr));
   amp1 = amplitudes;
   
   FREAD1(buff,diagrp);

   nsub0=buff.nsub;

   do 
   {  
      fpos=ftell(diagrp);
      if (!buff.delMark)
      {
          memcpy(amp1->dgrm, buff.dgrm0, sizeof(decayDiagram));
          amp1->gen=NULL;
          amp1->dim=1;
          amp2 = amp1;
          amp1 = (ampllist)m_alloc(sizeof(struct amplListStr));
          amp2->next = amp1;
      }
   }   while (1==FREAD1(buff,diagrp) && nsub0==buff.nsub );

   if(!feof(diagrp)) fseek(diagrp,fpos,SEEK_SET);
   
   if (amp1==amplitudes){ free(amp1); return NULL;}
   free(amp1);
   amp2->next = NULL;

   for(amp1 = amplitudes;amp1;amp1 = amp1->next)
   { int buff[MAXINOUT];
      
     InOutPrtclsNumb(amp1->dgrm,buff,0);
     if(amp1== amplitudes) 
     { for(i=0;i<nout;i++) cononout[i]=buff[i+nin];
       SORTARR(cononout,nout);
     }     
     smpl_linker(buff+nin,cononout,amp1->perm);
   }        
     
   for (i = 0; i < nout - 1; i++)
   { if (cononout[i] == cononout[i+1])
     {              
       for(amp1 = amplitudes;amp1;amp1 = amp1->next)
       {
         revers1(amp1->perm,rvrs);
         dogen(amp1->dgrm,&amp1->gen,&amp1->dim);
         for( cgen = amp1->gen;cgen;cgen = cgen->next)
         {
            multperm(rvrs,cgen->perm,cgen->perm);
            multperm(cgen->perm,amp1->perm,cgen->perm);
         }
       }
       i=nout;
     }  
   }
   return amplitudes;
}

ampllist  labeledDiagrams(FILE * diagrp)
{
  permut cononout;
  ampllist ampl, ampl_;
  grouplist group, group_;
  int n;
  
  ampl=readDiagrams(diagrp, cononout);

  allperm(cononout,&group);
  
  for(ampl_=ampl;ampl_;ampl_=ampl_->next)
  {
    doleftident(group,ampl_);
   
    while(ampl_->gen)
    { permlist  p_=ampl_->gen;
       ampl_->gen=(ampl_->gen)->next;
       free(p_);
    }
        
    for(group_=group,n=1; group_; n++,group_=group_->next)
    if(n==group_->left)
    {  permlist p=(permlist)malloc(sizeof(*p));   
       multperm(ampl_->perm,group_->perm,p->perm);       
       p->next=ampl_->gen;
       ampl_->gen=p;
    }       
  }
  
  clearGroop(group);
  return ampl ;
}

static int checkLimQ( particleNumType * diag1 ,particleNumType * diag2)
{ int n,i;
  int L=2*(nin+nout) -3;
  for(n=0; LimQ[n].who; n++)
  {  int p=LimQ[n].who;
     int ap=prtclbase1[p].anti;
     int N=LimQ[n].how;
     for(i=0;i<L;i++) 
     { int np=diag1[i];
       if(np<0 && (-np==p || -np==ap)) N--;
       if(!N) return 1; 
     }
     for(i=0;i<L;i++) 
     { int np=diag2[i];
       if(np<0 && (-np==p || -np==ap)) N--;
       if(!N) return 1; 
     }
     if(N!=1) return 1;   
  }
  return 0;
}

int  squaring(void)
{

 int       i;
 unsigned  n, m;
 shortstr  namesubproc;
 int       testsim;
 ampllist  amplitudes;
 permut    cononout;
 long      constrdiagr, maxdiagr;
 long      firstpos, nsdiagram;
 int       nsubcs;
 ampllist  amp1, amp2;

 grouplist sim, sim1, group;
 permut    rvrs;

 csdiagram    sqres;
 char format[2]={68,77};

   menup=fopen(MENUP_NAME,"rb"); 
   menuq=fopen(MENUQ_NAME,"wb"); 
   
   diagrp=fopen(DIAGRP_NAME,"rb");
   diagrq=fopen(DIAGRQ_NAME,"wb"); 

   maxdiagr = 1000;
   constrdiagr = 0;
   f_write(format,2,1,menuq);    
   nsubcs = 1;
   scrcolor(FGmain,BGmain);
   goto_xy(13,21); print("0      diagrams are constructed ");

   while(!feof(diagrp))
   {
      amplitudes=readDiagrams(diagrp,cononout );

      if (amplitudes)
      {
         nsdiagram = 0;

         firstpos = ftell(diagrq)/sizeof(csdiagram);
         allperm(cononout,&group);

         testsim=(group->next != NULL);

         amp1 = amplitudes;
         while (amp1 != NULL)
         {
            if (testsim) doleftident(group,amp1);
            amp2 = amp1;
            while (amp2 != NULL)
            {  
               if(LimQ[0].who && checkLimQ(amp1->dgrm,amp2->dgrm))
               { amp2=amp2->next; continue;}
                          
               memcpy(sqres.dgrm1,amp1->dgrm,sizeof(sqres.dgrm1));
               memcpy(sqres.dgrm2,amp2->dgrm,sizeof(sqres.dgrm2));
               if(amp1==amp2) sqres.mult=1; else sqres.mult=2;

               if (!testsim)
               {
                  smpl_linker(amp1->perm,amp2->perm,sqres.lnk);
                  sqres.del = 1;
                  sqres.status = 0;
                  sqres.nsub=nsubcs;
                  FWRITE1(sqres,diagrq);
                  nsdiagram++;
               }
               else
               {
                  doresident(group,amp2);
                  sim = group;
                  n = 1;
                  while (sim != NULL)
                  {
                     if (sim->res == n)
                     {
                        if (amp1 == amp2)
                        {
                           sqres.mult = 1;
                           revers1(sim->perm,rvrs);
                           sim1 = sim;
                           while (!eqprm(sim1->perm,rvrs))
                           {
                              sim1 = sim1->next;
                              if (sim1 == NULL) goto label_2;
                           }
                           if (sim1->res == n)
                              sqres.mult = 1;
                           else
                              sqres.mult = 2;
                        }
                        m = 1;
                        sim1 = sim->next;
                        while (sim1 != NULL)
                        {
                           if (sim1->res == n) m++;
                           sim1 = sim1->next;
                        }
                        sqres.del = (amp1->dim * amp2->dim) / m;
                        revers1(amp2->perm,rvrs);
                        multperm(sim->perm,rvrs,rvrs);
                        multperm(amp1->perm,rvrs,sqres.lnk);
                        
                        sqres.status = 0;
                        sqres.nsub=nsubcs;
                        FWRITE1(sqres,diagrq);
                        nsdiagram++;
                     }
label_2:             n++;
                     sim = sim->next;
                  }
               }
               amp2 = amp2->next;
            }
            amp1 = amp1->next;
         }

         proccessName(amplitudes->dgrm,namesubproc); 
	 wrt_menu(2,nsubcs,namesubproc,0,0,nsdiagram,firstpos);
         nsubcs++; 
         
         clearDiagrams(amplitudes);
         clearGroop(group);         
     
         constrdiagr += nsdiagram; 
         goto_xy(13,21); print("%5d",constrdiagr); 
                  
         if ((constrdiagr) > maxdiagr && !feof(diagrp))
          {  if(mess_y_n(3,17," Continue ?")) maxdiagr += 1000;  else goto label_3;}
      }  
   } 
label_3: 
   subproc_sq = nsubcs -1;
   fclose(menup);
   fclose(menuq);
   fclose(diagrp);
   fclose(diagrq);
	for (i = 17; i <= 24; i++) { goto_xy(1,i); clr_eol(); }
   if (constrdiagr == 0)
   { 
      messanykey(5,17," All process are marked as deleted"); 
      return 0; 
   } 
   else 
   { 
/*      be_be();  */
      return 1; 
   } 
} 
