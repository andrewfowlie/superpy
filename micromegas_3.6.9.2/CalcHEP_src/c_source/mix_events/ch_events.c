#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>

#include "ch_events.h"

#define maxFileOpen 100
int nFiles=0;





static void Rot3D(double cn , double fi_1, double psi, int ntot, double *pvect)
{
   double sn=sqrt(fabs(1-cn*cn));
   double E[3]={ cn, sn*sin(fi_1), sn*cos(fi_1)};
   double EPS[3][3]= { {  0,  -E[2], E[1]},
                     { E[2],   0 ,-E[0]},
                     {-E[1], E[0],  0  }  }; 

   double P[3][3]=  {   {1-E[0]*E[0],  -E[0]*E[1], -E[0]*E[2] },
                      { -E[1]*E[0], 1-E[1]*E[1], -E[1]*E[2] },
                      { -E[2]*E[0],  -E[2]*E[1],1-E[2]*E[2] } };
                      
                      
  double ROT[3][3];
  double x1,x2,x3,f3,cn2,sn2;
  
  int i,j,k;
  
  x1=psi;
  x2=M_PI;
  for(;;)
  { x3=0.5*(x1+x2);
    f3=x3-psi-sin(x3);
    if(f3<0) x1=x3;else x2=x3;
    if(fabs(f3)<1.E-4) break;
  }
    
  cn2=cos(x3),sn2=sin(x3);

  for(i=0;i<3;i++) for(j=0;j<3;j++) ROT[i][j]= (cn2-1)*P[i][j] + sn2*EPS[i][j]; 

      
  for(k=0;k<ntot;k++) 
  {  double * X= pvect+3*k;
     double Y[3]={X[0],X[1],X[2]}; 
     for(i=0;i<3;i++)  for(j=0;j<3;j++) Y[i]+=ROT[i][j]*X[j];  
     for(i=0;i<3;i++) X[i]=Y[i];
  }                       

}

static unsigned long _time=1;

eventfile_info * All=NULL;

decay_info * Decays=NULL;



static int closeLast(void)
{
   eventfile_info * a,*a0=NULL;
   decay_info * D;
   
   for(a=All;a;a=a->next) if(a->F) { if(!a0) a0=a; else { if(a0->ltime > a->ltime) a0=a;} }
   for(D=Decays;D;D=D->next) for(a= Decays->List ;a;a=a->next) 
                          if(a->F) { if(!a0) a0=a; else { if(a0->ltime > a->ltime) a0=a;} }
   if(a0){ a0->CurrentEventPos=ftell(a0->F); fclose(a0->F); a0->F=NULL; nFiles--; return 0;}
   return 1;
}

eventfile_info * initEventFile(char* fname)
{
   char buff[200];
   eventfile_info * Finfo;
   FILE*F;
   int i,n;
   int ntot;
 
   F=fopen(fname,"r");
   if(F==NULL) { printf("Can not open event file %s\n",fname); exit(1);} 
   
   Finfo=(eventfile_info *)malloc(sizeof(eventfile_info)); 

   Finfo->F=NULL;
   Finfo->ltime=0;
   Finfo->firstRd=1;     
   for(fscanf(F,"%s",buff);!feof(F); )
   {
      if(strcmp(buff,"#Type")==0) 
      { fscanf(F,"%d -> %d", &Finfo->Nin, &Finfo->Nout);
        ntot=Finfo->Nin+Finfo->Nout;        
      } else if(strcmp(buff,"#Initial_state")==0)
      {
        for(i=0;i<2;i++) { Finfo->inPID[i]=0;Finfo->inMom[i]=0; Finfo->pdf[i][0]=0;Finfo->inPID[i]=0; }
        if(Finfo->Nin==2)
        { 
          for(i=0;i<2;i++){fscanf(F," P%d_3=",&n); fscanf(F,"%lf",Finfo->inMom+n-1); }
          for(i=0;i<2;i++) 
          { fscanf(F," StrFun%d=",&n); 
            if(n>0 && n<=2) if(fscanf(F,"\"%[^\"]\" %d", Finfo->pdf[n-1], Finfo->inPID+n-1)!=2) Finfo->inPID[n-1]=0;
          }
        }  
      }
      else if(strcmp(buff,"#PROCESS")==0)
      {  int np;
         for(np=0;np<ntot;np++)
         { fscanf(F,"%d(%[^)])",Finfo->PIDs+np,Finfo->pName[np]);
           if(np==Finfo->Nin-1)fscanf(F," -> ");
         }  
      }  else if(strcmp(buff,"#MASSES")==0)
      { int i;
         for(i=0;i<ntot;i++) fscanf(F," %lf",&Finfo->pmass[i]);
      } else if(strcmp(buff,"#Cross_section(Width)")==0)
      {  if( fscanf(F," %lf",&Finfo->cs)!=1)
         { fprintf(stderr,"Error: unknown  cross section/width in %s\n",fname);
           exit(5);
         }            
      } else if(strcmp(buff,"#Number_of_events")==0)
      {  fscanf(F," %ld",&Finfo->nEvents);      
      } else if(strcmp(buff,"#Events")==0)
      {  fscanf(F,"%*[^\n]%*c");
         for(i=0;i<Finfo->Nin;i++) if(Finfo->inPID[i]==0) Finfo->inPID[i]=Finfo->PIDs[i];
         
         if(Finfo->cs<0) Finfo->cs*=-1; 
         Finfo->cEvent=1;
         Finfo->FirstEventPos=ftell(F);
         Finfo->CurrentEventPos=Finfo->FirstEventPos;
         Finfo->fileName=malloc(strlen(fname)+1);
         strcpy(Finfo->fileName,fname);
         fclose(F);
         return Finfo;
      }
      fscanf(F,"%s",buff);
   }
   return NULL; 
}

int readEvent(eventfile_info *Finfo, int *Nmom, double * mom, int * clr1, int * clr2, double *Qf,double *alphaQCD, int * w)
{ int n=0;

  if(Finfo->F==NULL)
  {
    if(nFiles>=maxFileOpen) closeLast();
    Finfo->F=fopen(Finfo->fileName,"r");
    if(Finfo->F==NULL) { printf("can not open file %s\n", Finfo->fileName); exit(2);}
    fseek(Finfo->F, Finfo->CurrentEventPos,SEEK_SET); 
    nFiles++;
  }

  Finfo->ltime=_time++;
  for(;;)
  {
     if(Finfo->cEvent > Finfo->nEvents)
     { Finfo->cEvent=1;
       fseek(Finfo->F, Finfo->FirstEventPos, SEEK_SET);
       Finfo->CurrentEventPos=Finfo->FirstEventPos;  
       if(Finfo->Nin==2) 
       {
          fprintf(stderr,"Error: File %s : no more events.\n", Finfo->fileName);
          return 1;
       }  else if(Finfo->firstRd==1)
       { 
         fprintf(stderr,"Warning: File %s : no more events. Reread from the beginning. \n", Finfo->fileName);
         Finfo->firstRd=0;   
       }
     }  
     Finfo->cEvent++;
     fscanf(Finfo->F,"%d",w);     
     for(n=0;1==fscanf(Finfo->F," %lf",mom+n);n++);
     *Nmom=n;
     if(2!=fscanf(Finfo->F,"| %lf %lf",Qf,alphaQCD)) return 2;
     if(Finfo->Nin==1 && Finfo->firstRd==0)  Rot3D(2*(drand48()-0.5) ,2*M_PI*drand48(),M_PI*drand48(),Finfo->Nout ,mom);  
     for(n=0;2==fscanf(Finfo->F," (%d %d)",clr1+n,clr2+n);n++);
     clr1[n]=0; clr2[n]=0;
     return 0;
  }
}

static void cleanList(eventfile_info * list)
{  eventfile_info * Finfo;
   for(;list;)
   {
     if(list->F)fclose(list->F); 
     free(list->fileName);
     Finfo=list;
     list=list->next;
     free(Finfo);
   }
}


void closeevents_(void)
{ decay_info * D;

//  totCS=0;
  nFiles=0; _time=1;
  cleanList(All); All=NULL;
  for(;Decays;)
  {
    cleanList(Decays->List);
    D=Decays;
    Decays=Decays->next;
    free(D);
  }   
}
