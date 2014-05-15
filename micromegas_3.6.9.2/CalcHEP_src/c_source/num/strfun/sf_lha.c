/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>
#include <dirent.h>

#include "files.h"
#include "interface.h"
#include "subproc.h"
#include "chep_crt.h"
#include "strfun.h"
#include "n_calchep_.h"
#include "lha.h"
#include "alphas2.h"
#include "../../../include/VandP.h"

#include "sf_lha.h"
static double xMin[2]={0,0},xMax[2]={0,0},qMin[2]={1,1},qMax[2]={1E10,1E10};

/*static char * pdfName[2]={NULL,NULL}; */

static int nGroup[2]={3,3}, nSet[2]={41,41}, sgn[2]={1,1};

static int parton[15]={83,81,2,1,4,3,5,21,-5,-3,-4,-1,-2,-81,-83};

static int pnum[2]={0,0};


static char *param=NULL;

static double sc2=0.22*0.22;

static char* dataPath=NULL;
static char  fname[2][30]={"",""};
static int   setNum[2]={0,0};


int mc_lha(int i) { return 2212*sgn[i-1];}


static int initLHA(void)
{ 
  int i;  
  char buff[500];
  
  if(dataPath) return 1;
  
  getdatapath(buff,500);
  for(i=499; i>=0 && buff[i]==' '; i--) ;
  if(i<0) return 0;
  buff[i+1]=0;
  dataPath=malloc(strlen(buff)+1);
  strcpy(dataPath,buff);
  return 1;  
}

int p_lha(int * pNum) 
{  
  int i;
  if(!initLHA())return 0;
  for(;*pNum;pNum++) 
  { for(i=0;i<15;i++) if(*pNum==parton[i]) break;
    if(i==15) return 0;
  }
  return 1;
}

void n_lha(int i, char *name) 
{  int i1=i-1;
   if(strlen(fname[i1]))  sprintf(name,"LHA:%s:%d:%d",fname[i1],setNum[i1],sgn[i1]);
   else  strcpy(name,"LHA:"); 
}

int r_lha(int i, char *name)
{ int i1=i-1;
  char txt[50];
  int max;
  if(!initLHA())return 0; 
  if(3!=sscanf(name,"LHA:%[^:]%*c%d:%d",txt,setNum+i1,sgn+i1)) return 0;
  if(abs(sgn[i1])!=1) return 0; 
  initpdfsetbynamem(&i,txt,strlen(txt));
  if(i==0) return 0;
  numberpdfm(&i,&max);
  if(max<setNum[i1] || setNum[i1]<0) return 0;
  initpdfm(&i,setNum+i1,xMin+i1,xMax+i1,qMin+i1,qMax+i1);
  strcpy(fname[i1],txt);
  return 1;
}


int init_lha(int i,double * be, double * mass) 
{  int k;
   int N,N1,N2;
   pinf_int( Nsub,1,NULL,&N1);
   pinf_int( Nsub,2,NULL,&N2);

   if(i==1) N=N1; else N=N2;
   if(abs(N1)>80 &&  abs(N2)>80) {if(N>0) N-=80; else N+=80;}
    
   *mass=0.9383;
   *be=1;
   sf_alpha=&(alpha_lha);   
   for(k=0;k<15;k++) if(parton[k]==N) {pnum[i-1]=N; return 1;}
   pnum[i-1]=0; 
   return 0;
}
 
static int filter (const struct dirent * dp)
{  
   char *ch ;
   
/*   if(dp->d_type == DT_DIR) return 0; */

   ch=strstr(dp->d_name,".LHpdf");
   if(ch && ch[6]==0) return 1;

   ch=strstr(dp->d_name,".LHgrid");
   if(ch && ch[7]==0) return 1;
   
   return 0;
}


static char * lhaMenu(void)
{
  struct  dirent **namelist;
  int N,i,width;
  char * menutxt=NULL;

  if(!dataPath) return NULL;

  N=scandir(dataPath, &namelist,filter ,alphasort);   
  if(N<=0) return NULL;
                              
  for(i=0,width=0;i<N;i++) 
  {  int l= strlen(namelist[i]->d_name);
     if(width<l) width=l;
  }
  menutxt=malloc(N*(width+1)+2);
  menutxt[0]=width+1; menutxt[1]=0;
  for(i=0;i<N;i++)
  {
    sprintf(menutxt+1+(width+1)*i," %-*.*s",width,width,namelist[i]->d_name);
    free(namelist[i]);
  }
  free(namelist);
  menutxt[N*(width+1)+1]=0;
  return menutxt;
}  

                        
int m_lha(int i,int*pString)
{ 
  int size=2+8*500;
  int size_=1;
  void *pscr=NULL;
  void *pscr0=NULL;
  static int n1=0;
  char * strmen=lhaMenu(); 
  int i1=i-1;

  if(!strmen) return 0;
  

  int n0=1,k,l;

  if(n1==0 && strlen(fname[i1]))
  { char *ch=strstr(strmen,fname[i1]);
    if(ch) n1= 1+(ch-strmen)/strmen[0];
  }
  menu1(5,7,"LHAlib menu",strmen,"",&pscr,&n1);
  if(n1)
  { char buff[50];
    sscanf(strmen+1+strmen[0]*(n1-1),"%s",buff);
    if(strcmp(buff,fname[i1]))
    { strcpy(fname[i1],buff);
      setNum[i1]=0;
      initpdfsetbynamem(&i,buff,strlen(buff));
      initpdfm(&i,setNum+i1,xMin+i1,xMax+i1,qMin+i1,qMax+i1);
    }
  }  
  else 
  { fname[i1][0]=0;
    setNum[i1]=0;
    sgn[i1]=0;
    return 0;
  } 


  for(;n0!=0 && n0!=3;)
  {  char buff[50];
     int nMax;
     char strmen0[]="\030"
                    " Set = 0                "   
                    " Proton                 "
                    " OK                     ";

     numberpdfm(&i,&nMax);
     if(nMax>1) improveStr(strmen0,"Set = 0","Set = %d [0,%d]",setNum[i1],nMax);
     else       improveStr(strmen0,"Set = 0","Set = 0 (only)");
     
     if(sgn[i1]<0) improveStr(strmen0,"Proton","%s","antiProton");
     
     menu1(5,10,"",strmen0,"",&pscr0,&n0);
     switch(n0) 
     { 
       case 1: if(nMax>1) 
               { correctInt(50,12,"Enter new value ",setNum+i1,1);
                 if(setNum[i1]<0) setNum[i1]=0;
                 if(setNum[i1]>nMax) setNum[i1]=nMax;
                 initpdfm(&i,setNum+i1,xMin+i1,xMax+i1,qMin+i1,qMax+i1);
               }   
               break;
       case 2: sgn[i1]=-sgn[i1]; break;
       case 3: put_text(&pscr0); break;
     }
  }
  
  put_text(&pscr); 
  free(strmen);
  return 1;
}

double c_lha(int i, double x, double q)
{
  double f[14];
  int p;
  int i1=i-1;
  double z;
  
  p=pnum[i1];

  if(x<xMin[i1]) x=xMin[i1]; else if(x>xMax[i1]) x=xMax[i1];
  if(q<qMin[i1]) q=qMin[i1]; else if(q>qMax[i1]) q=qMax[i1];
  
    
  evolvepdfm(&i,&x,&q,f);  
  if(sgn[i1]<0) p=-p;
/*
  switch(p)
  { case 81:           return (f[7]*(1-sc2) + f[9]*sc2)/x;
    case 83:           return (f[7]*sc2 + f[9]*(1-sc2))/x;
    case 2 :           return f[8]/x;
    case 1 :           return f[7]/x;
    case 3 : case -3 : return f[9]/x;
    case 4 : case -4 : return f[10]/x;
    case 5 : case -5 : return f[11]/x; 
    case 21: case -21: return f[6]/x;
    case -1:           return f[5]/x;
    case -2:           return f[4]/x;
    case -81:          return (f[5]*(1-sc2) + f[9]*sc2)/x;
    case -83:          return (f[5]*sc2 + f[9]*(1-sc2))/x;
  } 
*/
  switch(p)
  { case 81:           z=(f[7]*(1-sc2) + f[9]*sc2)/x; break;
    case 83:           z=(f[7]*sc2 + f[9]*(1-sc2))/x; break;
    case 2 :           z=f[8]/x; break;
    case 1 :           z=f[7]/x; break;
    case 3 : case -3 : z=f[9]/x; break;
    case 4 : case -4 : z=f[10]/x; break;
    case 5 : case -5 : z=f[11]/x; break; 
    case 21: case -21: z=f[6]/x; break;
    case -1:           z=f[5]/x; break;
    case -2:           z=f[4]/x; break;
    case -81:          z=(f[5]*(1-sc2) + f[9]*sc2)/x; break;
    case -83:          z=(f[5]*sc2 + f[9]*(1-sc2))/x; break;
  } 
  if(z<0) return 0;
//  if(z<=0) printf("x=%E q=%E z=%E sc2=%E   \n",x,q,z,sc2); 
  return z;
  
}

