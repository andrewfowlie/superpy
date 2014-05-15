/*
 Copyright (C) 1997, Alexander Pukhov
*/
#include<stdio.h>
#include "process.h"
#include "sortarr.h"
#include "process.h"
#include "diagrams.h"

void  mkverts(decayDiagram diag,vampl* vlist) /*  transformation of diagr representation */
{    
   int i,j,k, particle;
   int newVert, currVert, extNum, restIn; 
   edgeinvert *w;                    
   
   for(i=0;i<maxvert;i++)
   { vlist->valence[i]=3;
      for (j = 0; j < MAXVALENCE; j++) 
      {  w= &( vlist->vertlist[i][j]);
         w->lorentz=0; 
         w->partcl = 0;
         w->prop=0;
         w->moment=0; 
      } 
   }

   w = &vlist->vertlist[0][0];     /*  process 1st edge  */ 
   particle=-diag[0];
   w->partcl = prtclbase[particle-1].anti; 
   w->prop = IN_PRTCL;
   if(polarized(1,particle)) w->prop+=PLR_PRTCL; 
   w->link.vno = nullvert; 
   w->link.edno = 0; 
   vlist->outer[0].vno =  0; 
   vlist->outer[0].edno = 0; 
    
  /*  *** main cycle *********  */ 

   for (i=1, currVert=0, newVert=1, extNum=1,restIn=nin-1; currVert>=0; i++)
   { particle=diag[i];
     if(particle)
     {   if (particle < 0) 
         {
           for(k=0;vlist->vertlist[currVert][k].partcl;k++);
           w= &vlist->vertlist[currVert][k]; 
           w->partcl = -particle;
           w->link.vno = newVert ; 
           w->link.edno = 0; 
      
           w = &vlist->vertlist[newVert][0]; 
           w->partcl = prtclbase[-particle-1].anti; 
           w->link.vno = currVert; 
           w->link.edno = k; 
      
           currVert = newVert++;   /*  next vert to be expanded  */ 
         }
         else  
         {                                                    
           for(k=0;vlist->vertlist[currVert][k].partcl;k++);                           

           vlist->outer[extNum].vno = currVert;
           vlist->outer[extNum].edno = k;

           w = &vlist->vertlist[currVert][k];
           w->partcl = particle;                                               
           if(restIn)
           { restIn--; 
             w->prop=IN_PRTCL;
             if(polarized(2,prtclbase1[w->partcl].anti)) w->prop+=PLR_PRTCL; 
           } else w->prop=OUT_PRTCL;    
           w->link.vno = nullvert;                                             
           w->link.edno = (extNum++);                                              

           while(currVert>=0 && 
                 (vlist->vertlist[currVert][vlist->valence[currVert]-1].partcl) 
                ) currVert--;
         }
     } else  {vlist->valence[currVert]++; }
   }

   vlist->outno = extNum; 
   vlist->size = newVert; 
   
   for(i=0;i<vlist->size;i++) 
   { for(j=vlist->valence[i];j<MAXVALENCE;j++)
     { vlist->vertlist[i][j].link.vno=nullvert;
       vlist->vertlist[i][j].link.edno=0;
       vlist->vertlist[i][j].nextvert=vlist->vertlist[i][j].link;
     }
     for(j=0;j<vlist->valence[i];j++)
     { vlist->vertlist[i][j].nextvert=vlist->vertlist[i][j].link;
        vlist->vertlist[i][j].nextvert.edno++;
        if(vlist->vertlist[i][j].link.vno!=nullvert)
        vlist->vertlist[i][j].nextvert.vno++;
     }          
   }
}   /* ==========mkVerts =========== */ 



void printDiagram(vampl * vlist)
{

int i,j;

for(i=0;i< vlist->size;i++)
{
  for(j=0;j< MAXVALENCE && vlist->vertlist[i][j].partcl ; j++)
  { edgeinvert *v=&vlist->vertlist[i][j];
    printf("vertex=%d, line=%d, prtcl=%s  link(%d,%d) moment=%d,lorentz=%d prop=%d \n",
 i+1,j+1,prtclbase[v->partcl-1].name,
 v->nextvert.vno,v->nextvert.edno,v->moment,v->lorentz, v->prop
           ); 
   }

  printf("\n");
}

 for(i=0;i<vlist->outno;i++)
 
printf("%d=(%d,%d) ",i+1, vlist->outer[i].vno,vlist->outer[i].edno); 

 printf("\n\n");
}


void printCsDiagram(vcsect * vlist)
{

int i,j;

for(i=0;i< vlist->sizet;i++)
{
  for(j=0;j< vlist->valence[i] ; j++)
  { edgeinvert *v=&vlist->vertlist[i][j];
  printf("vertex=%d, line=%d, prtcl=%s  link(%d,%d) moment=%d,lorentz=%d prop=%d \n",
 i+1,j+1,prtclbase[v->partcl-1].name,

 v->nextvert.vno,v->nextvert.edno, v->moment,v->lorentz, v->prop
 
  ); 
  
  }

  printf("\n");
}


 printf("\n\n");
}


                  



static void  mkcsections(csdiagram* diagr,vcsect* vcs)  
/*  constructs vertex-oriented structure for cross section diagramm  */ 
{ 
   int i,j, shift; 
   vertlink   buff[MAXINOUT]; 
   vampl        va1, va2; 


   vcs->symnum   = diagr->mult; 
   vcs->symdenum = diagr->del; 
   
   mkverts(diagr->dgrm1,&va1); 
/*printDiagram(&va1); */
   vcs->sizel = va1.size; 
   for (i = 0; i < va1.size; i++) 
   { vcs->valence[i]=va1.valence[i];
     for(j=0;j<MAXVALENCE;j++) 
     vcs->vertlist[i][j] =  va1.vertlist[i][j]; 
   }

   mkverts(diagr->dgrm2,&va2); 
/*printDiagram(&va2);*/
   for(i=0;i<va2.outno;i++) buff[i]=va2.outer[i]; 
   for(i=nin;i<va2.outno;i++)
   va2.outer[i]=buff[nin+ diagr->lnk[i-nin]-1 ];

/*printDiagram(&va2);*/


   shift=va1.size;
   vcs->sizet = va2.size+shift;
   for (i = 0; i < va2.size; i++) 
   { vcs->valence[i+shift]=va2.valence[i];
     for(j=0;j<MAXVALENCE;j++) 
     vcs->vertlist[i+shift][j]= va2.vertlist[i][j]; 
   }
   
   for(i=vcs->sizel;i<vcs->sizet;i++)
   for(j=0;j<vcs->valence[i];j++)
   {  edgeinvert * v=&vcs->vertlist[i][j];  
      v->partcl = prtclbase[v->partcl-1].anti; 
      if(v->nextvert.vno !=nullvert) v->nextvert.vno +=shift;
      if(v->link.vno !=nullvert) v->link.vno +=shift;
   }


   for(i=0;i< va1.outno;i++)
   { vcs->vertlist[va1.outer[i].vno][va1.outer[i].edno].nextvert.vno=
   va2.outer[i].vno+shift+1;
   vcs->vertlist[va1.outer[i].vno][va1.outer[i].edno].nextvert.edno=
      va2.outer[i].edno+1;

   vcs->vertlist[va2.outer[i].vno+shift][va2.outer[i].edno].nextvert.vno=
   va1.outer[i].vno+1;
   vcs->vertlist[va2.outer[i].vno+shift][va2.outer[i].edno].nextvert.edno=
   va1.outer[i].edno+1;


   vcs->vertlist[va1.outer[i].vno][va1.outer[i].edno].link=va2.outer[i];
   vcs->vertlist[va1.outer[i].vno][va1.outer[i].edno].link.vno +=shift;

   vcs->vertlist[va2.outer[i].vno+shift][va2.outer[i].edno].link=va1.outer[i];
   }
} 

  /* ====================================== mkCSections=======  */ 
  /*    labelling diagramm with momentums       */ 
  /*    This procedures labels cross-section diagramm with               */ 
  /*    momenta and lorentz ind.                                         */ 
  /*    Momenta coded with numbers stored in                             */ 
  /*       vCSect.vertList[vertex,edge].moment                           */ 
  /*    Positive number means outgoing momentum, negative - ingoing.     */ 
  /*    Linear relations are stored in vCSect.linear[vertex]             */ 
  /*    as array of momenta with zero sum. Independent - in 1st place    */ 
  /*    Linear relations must be supplemented with relation, describing  */ 
  /*    dependence of external momenta of the amplitude.                 */ 
  /*    Lorentz ind coded with bytes are stored in                       */ 
  /*     vCSect.vertList[vertex,edge].lorentz                            */ 
  /* ******************* A.Taranov    24.07.89 ------------------------- */ 



static void   assignlor(vcsect* vcs)
{
   int  lorcnt=1,i,j;

   for (i = 0; i < vcs->sizet; i++) 
   for (j = 0; j < vcs->valence[i]; j++) 
   { edgeinvert *  v=&vcs->vertlist[i][j],
     * av=&vcs->vertlist[v->nextvert.vno-1][v->nextvert.edno-1]; 

      if (vectorp(v->partcl)&& (v->nextvert.vno-1 > i) )
      { 
         v->lorentz = lorcnt; 
        av->lorentz = lorcnt; 
        ++(lorcnt); 
      } 
   } 
}   /*  assignLor  */ 


/* ********************main program ************************************** */ 

static void  labelmom(vcsect* vcs1)
{  int  i, j, momcnt, inmom,virtmom, outlst[MAXINOUT]; 
    
   /* ------------------ make list out particles ---------------------- */
   momcnt=0;
   for (i = 0; i < vcs1->sizel; i++)
   for (j = 0; j < vcs1->valence[i]; j++)
   { edgeinvert *  ln=&vcs1->vertlist[i][j];
     if (ln->prop & OUT_PRTCL ) outlst[momcnt++] = ln->partcl;
   }
   /* ---------------- sort List out particles -------------------------- */
     SORTARR(outlst,nout);
   /* ------------------ assign out momentum ----------------------------- */ 
   inmom=1;
   virtmom=nin+nout+1;
   for (i = 0; i < vcs1->sizet; i++) 
   for (j = 0; j < vcs1->valence[i]; j++) 
   { edgeinvert *  v=&vcs1->vertlist[i][j];
     edgeinvert * av=&vcs1->vertlist[v->nextvert.vno-1][v->nextvert.edno-1]; 
     if (i<vcs1->sizel && ( v->prop & OUT_PRTCL))     /* out particles */
     {  int   nprtcl  = v->partcl; 
        int pos=0;
   
        while (outlst[pos] != nprtcl) pos++;
        outlst[pos] = 0;
        v->moment  = pos+nin+1;    
        av->moment = -v->moment; 
     }else if(i<vcs1->sizel && (v->prop & IN_PRTCL))  /* in particles */
     {
        v->moment = -inmom++;   
        av->moment = -v->moment; 
     }else if ( !(v->prop & (IN_PRTCL|OUT_PRTCL))  && i<v->nextvert.vno-1)
     {   v->moment=virtmom++;     
        av->moment=-v->moment;
     }  
     
   }  /* assignOutMom */ 

}   /* ================================================== labelMom======== */ 



void   transfdiagr(csdiagram* diag,vcsect* vcs)
{ 

   mkcsections(diag,vcs);   /*  make cross section diagramm  */ 
   labelmom(vcs);           /*  label it with moments        */ 
   assignlor(vcs);          /*  label it with  Lorentz ind   */ 
#ifdef DEBUG                            
   printCsDiagram(vcs);
#endif
}   /*   transfDiagr  */ 


void InOutPrtclsNumb( decayDiagram a, int * numb, int sort )
{
  int i=1, k=1;

  numb[0]=-a[0]; 
  for (;i<2*(nin+nout)-3;i++) {if(a[i]>0) numb[k++] = a[i];}
  if(nin==2) numb[1]=prtclbase[numb[1]-1].anti;
    
  if(sort) SORTARR((numb+nin) ,nout)

}  

void proccessName(decayDiagram a, char * txt )
{ int k;
  int buff[MAXINOUT];
  InOutPrtclsNumb( a, buff,1);

  txt[0]=0;

  for(k=0;k<nin;k++)
  {  strcat(txt, prtclbase[buff[k]-1].name);
     if(k<nin-1) strcat(txt,",");
  }
  strcat(txt,"->");
   
  for(k=nin;k<nin+nout;k++)
  {  strcat(txt, prtclbase[buff[k]-1].name);
     if(k<nin+nout-1) strcat(txt,",");
  }
}

void decompose( vcsect vcs,  vampl *left,  vampl * right)
{  int i,j;

   left->size=vcs.sizel;
   for(i=0;i<vcs.sizel;i++)       
   {                              
     left->valence[i]=vcs.valence[i];
     for(j=0;j<vcs.valence[i];j++) 
     { left->vertlist[i][j]=vcs.vertlist[i][j];
       if(vcs.vertlist[i][j].link.vno >= vcs.sizel)
       { int m=vcs.vertlist[i][j].moment;
         if(m<0) m=-m;
         m--; 
         left->outer[m].vno=i;
         left->outer[m].edno=j;
         left->vertlist[i][j].link.vno=nullvert;
         left->vertlist[i][j].link.edno=m;
       }
     }
   }                              
      
   right->size=vcs.sizet-vcs.sizel;    
   for(i=vcs.sizel;i<vcs.sizet;i++)
   {  int i_=i-vcs.sizel;
      right->valence[i_]= vcs.valence[i];
      for(j=0;j<vcs.valence[i];j++)
      {  right->vertlist[i_][j]=vcs.vertlist[i][j];
         right->vertlist[i_][j].partcl=
                prtclbase[vcs.vertlist[i][j].partcl-1].anti;
         right->vertlist[i_][j].moment *=-1;
         if(vcs.vertlist[i][j].link.vno < vcs.sizel)
         {  int m=vcs.vertlist[i][j].moment; 
            if(m<0) m=-m;
            m--;
            right->outer[m].vno=i_;
            right->outer[m].edno=j;
            right->vertlist[i_][j].link.vno=nullvert;
            right->vertlist[i_][j].link.edno=m;
         } else  right->vertlist[i_][j].link.vno -=vcs.sizel;
      }
   }
}
