/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>

#include "pdt.h"
#include "files.h"
#include "interface.h"
#include "subproc.h"
#include "chep_crt.h"
#include "strfun.h"
#include "n_calchep_.h"
#include "alphas2.h"
#include "sf_pdt.h"


static pdtList * allPDT=NULL;
static pdtList * curPDT[2]={NULL,NULL};

static pdtStr * pdtData[2][20]={{NULL},{NULL}};
static int  posAux[2]={0,0};
static int  pos[2]={0,0};

static int alphaMode=0;

static double alpha_pdt(double q ){ return interAlpha(q, pdtData[alphaMode-1][0]);}

static void findAllPDT(void)
{  int k,l;
   char * dNames[3];
   char fname[STRSIZ],rootPdt[STRSIZ];
   long pNum1,pNum2;
   
   if(allPDT) return;

   sprintf(rootPdt,"%s/pdTables", pathtocalchep);
   
   dNames[0]="..";
   dNames[1]=".";
   dNames[2]=rootPdt;

   for(k=0 ;k<3;k++)
   {  DIR *dirPtr=opendir(dNames[k]);
      struct dirent * dp;

      if(!dirPtr) continue;
 
      while(dp=readdir(dirPtr))
      { char *c=dp->d_name;
        int l=strlen(c);
        if(l>=4 && strcmp(c+l-4,".pdt")==0) 
        {  sprintf(fname,"%s%c%s",dNames[k],f_slash,c); 
           makePdtList(fname, &allPDT);
        }
      } 
      closedir(dirPtr);
   }
} 


int mc_pdt(int i) {if(curPDT[i-1])return curPDT[i-1]->beamP; else return 0;}

int p_pdt(int * pNum) 
{  
   pdtList * L;
   findAllPDT();
   for(L=allPDT;L;L=L->next) if(checkPartons(pNum,  L)) return 1;
   return 0;       
}

void n_pdt(int i, char *name) 
{  i--;  
   if(curPDT[i]) sprintf(name,"PDT:%s",curPDT[i]->name); 
             else strcpy(name,"PDT:"); 
}

static int updateData(int i,pdtList * p)
{

   if( curPDT[i] && p!=curPDT[i] ) 
   { int * x=curPDT[i]->items;
     for(;*x;x++) { pdtStr**pd= &(pdtData[i][*x-1]);
                    if(*pd) {freePdtData(*pd); free(*pd); *pd=NULL;}
                  }                                     
   }

   if(p && p!=curPDT[i])
   { int * x= p->items;
     int k=0;
     for(;*x;x++) { pdtStr**pd= &(pdtData[i][*x-1]);
                    if(!*pd) { *pd=malloc(sizeof(pdtStr));
                                getPdtData(p->file,*x,*pd);
                             }
                  }                                
   }

   curPDT[i]=p;

   return 0;
}

int r_pdt(int i, char *name)
{  int k; 
   pdtList * p;
   i--;
   if(name!=strstr(name,"PDT") ) return 0;
   findAllPDT();
   for(p=allPDT;p;p=p->next) if(strcmp(name+4,p->name)==0)break; 

   updateData(i,p);

   return 1;
}


int init_pdt(int i,double * be, double * mass) 
{  
   pdtList *list=NULL;
   pdtList *list_=NULL;
   int ret_code=0; 
   int pNum,pNumAux, N1, N2;
   int k;
   i--;
   pinf_int(Nsub,1,NULL,&N1);
   pinf_int(Nsub,2,NULL,&N2);
   if(i) {pNum=N2; N2=N1; N1=pNum;} else pNum=N1;
   
   pNumAux=0;
   
   switch(N1)
   {  case  81: pNum= 1;  if(abs(N2)==2 ||abs(N2)==4) pNumAux=3;  break;
      case  83: pNum= 3;  if(abs(N2)==2 ||abs(N2)==4) pNumAux=1;  break;
      case -81: pNum=-1;  if(abs(N2)==2 ||abs(N2)==4) pNumAux=-3; break;
      case -83: pNum=-3;  if(abs(N2)==2 ||abs(N2)==4) pNumAux=-1; break;
   }     
   for(k=0;curPDT[i]->partons[k];k++) if(curPDT[i]->partons[k]==pNum)
   {  pos[i]=curPDT[i]->items[k]-1;   
       if(pdtData[i][0]->alpha) { alphaMode=i+1; sf_alpha=&(alpha_pdt);}
       if(pdtData[i][pos[i]]->pow1<0) *be=pdtData[i][pos[i]]->pow1+1; else *be=1.;
       *mass=pdtData[i][pos[i]]->mass;
       break;
   }
   if(pNumAux==0) posAux[i]=0; else 
   for(k=0;curPDT[i]->partons[k];k++) if(curPDT[i]->partons[k]==pNumAux)
   { posAux[i]=curPDT[i]->items[k]-1; break;}
   return 1;
}   


#define WIDTH 30
int m_pdt(int i,int*pString)
{
   pdtList *list_; 
   
   int k=0;
   char * strmen; 
   int pNum;

   findAllPDT();
   if(!allPDT) return 0;
   
   i--;
   pinf_int(Nsub,i+1,NULL,&pNum);
 
   for(list_=allPDT;list_;list_=list_->next) if(checkPartons(pString,list_)) k++;
 
   strmen=malloc(2+WIDTH*k);
   strmen[0]=WIDTH;strmen[1]=0; 
   for(list_=allPDT; list_; list_=list_->next) if(checkPartons(pString,list_)) 
      sprintf(strmen+strlen(strmen)," PDT:%-*.*s",WIDTH-5,WIDTH-5,list_->name);
   k=0;
   menu1(80-4-WIDTH,7,"",strmen,"",NULL,&k);
   free(strmen);
   if(k==0) return 0;  
   k--;   
   for(list_=allPDT; k ; list_=list_->next) if(checkPartons(pString,list_))k--;  ;
   updateData(i,list_);
   return 1;
}

double c_pdt(int i, double x, double q)
{  double r;
   i--;
   r=interFunc( x, q , pdtData[i][pos[i]]);
   if(pdtData[i][pos[i]]->pow0)    r*=pow(x,  pdtData[i][pos[i]]->pow0);
   if( pdtData[i][pos[i]]->pow1>0) r*=pow(1-x,pdtData[i][pos[i]]->pow1);
   if(posAux[i])
   {  double  sc2=0.221*0.221;
      double rAux=interFunc( x, q , pdtData[i][posAux[i]]); 
      if(pdtData[i][posAux[i]]->pow0)    rAux*=pow(x,  pdtData[i][posAux[i]]->pow0);
      if( pdtData[i][posAux[i]]->pow1>0) rAux*=pow(1-x,pdtData[i][posAux[i]]->pow1);
      r=r*(1-sc2)+rAux*sc2;
   }

   if(r<0)r=0;
   return r;
}
