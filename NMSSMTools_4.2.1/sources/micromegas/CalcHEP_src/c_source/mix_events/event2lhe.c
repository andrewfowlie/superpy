#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<unistd.h>
#include<time.h>
#include"SLHAplus.h"
#include"../../include/version.h"

#include "ch_events.h"
#include "event2pyth.h"


extern int oneFileInit( char* filename);

#define  LPRUP   1
#define  NPRUP   1
#define  PDFSUP -1
#define  PDFGUP -1
#define  AQEDUP -1
#define  IDWTUP  3
#define  XERRUP  0.
#define  XMAXUP  1.
                
static  char* pdg2name(int pdg)
{ static char buff[40];
  switch(pdg)
  { 
    case  2212: return "p"; 
    case -2112: return "anti-p";
    case   11: return "e-"; 
    case  -11: return "e+";
    case   13: return "mu-"; 
    case  -13: return "mu+";
    case   22: return "gamma";
    default: sprintf(buff,"PDG(%d)",pdg);
             return buff;
  }
}
             
int main(int argc,char ** argv)
{
  int N,NEV,MAXEVENTS,II,J,K,err;
  double cs;
  FILE *F=stdout;
  long posNevents, posSize;
  int SLHA=0;
  int ok;
  double Lumi;
  eventfile_info *Finfo;
   
  if(argc<2){ printf("%s needs one argument: name of event file in CalcHEP format.\n",argv[0]); exit(1); }

  Finfo=initEventFile(argv[1]);
  if(!Finfo) exit(2);
  if( Finfo->cs==0 || Finfo->nEvents==0) 
  {  printf("There are no events in file %s \n",argv[0]);
     exit(3);
  }   

  fprintf(F,"<LesHouchesEvents version=\"1.0\">\n");
  fprintf(F,"<!--\n");
  fprintf(F," CalcHEP event file transformed to the LHE format\n");
  fprintf(F,"-->\n");
  fprintf(F,"<header>\n");
  fprintf(F,"<hepml>\n");
  fprintf(F,"<samples xmlns=\"http://mcdb.cern.ch/hepml/0.2/\"\n");
  fprintf(F,"    xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n");
  fprintf(F,"    xsi:schemaLocation=\"http://mcdb.cern.ch/hepml/0.2/ http://mcdb.cern.ch/hepml/0.2/hepml.xsd\">\n");
  fprintf(F,"    <description>\n");
  fprintf(F,"        <process>\n");
  fprintf(F,"            <beam1>\n");
  fprintf(F,"                <particle KFcode=%d>\"%s\"</particle>\n",Finfo->inPID[0],pdg2name(Finfo->inPID[0]));
  fprintf(F,"                <energy unit=\"GeV\">%f</energy>\n",fabs(Finfo->inMom[0]));
  fprintf(F,"                <pdf name= \"%s\"></pdf>\n",Finfo->pdf[0]);  
  fprintf(F,"            </beam1>\n");
  if(Finfo->Nin==2)
  { 
    fprintf(F,"            <beam2>\n");
    fprintf(F,"                <particle KFcode=%d>\"%s\"</particle>\n",Finfo->inPID[1],pdg2name(Finfo->inPID[1]));
    fprintf(F,"                <energy unit=\"GeV\">%f</energy>\n", fabs(Finfo->inMom[1]));
    fprintf(F,"                <pdf name= \"%s\"></pdf>\n",Finfo->pdf[1]);  
    fprintf(F,"            </beam2>\n");
  } 

  if(Finfo->Nin==2) fprintf(F,"            <crossSection unit=\"pb\">%f</crossSection>\n",Finfo->cs);
  else              fprintf(F,"            <width unit=\"GeV\">%f</width>\n",Finfo->cs);
  
//  fprintf(F,"            <subprocesses>\n");
//?  printProcInfo(F);

//  fprintf(F,"                  <FactorisationScale>\n");
//  fprintf(F,"                      <plain></plain>\n");
//  fprintf(F,"                      <Latex></Latex>\n");
//  fprintf(F,"                  </FactorisationScale>\n");  

//  fprintf(F,"            </subprocesses>\n");
  fprintf(F,"        </process>\n");
  fprintf(F,"    </description>\n");
  fprintf(F,"</samples>\n");
  fprintf(F,"</hepml>\n");

    fprintf(F,"</header>\n");	

    fprintf(F,"<init>\n");
    fprintf(F," %5d %5d %18.11E %18.11E %5d %5d %5d %5d %5d %5d\n",
      Finfo->inPID[0],Finfo->inPID[1],fabs(Finfo->inMom[0]),fabs(Finfo->inMom[1]), 
      PDFGUP,PDFGUP ,PDFSUP ,PDFSUP, IDWTUP,NPRUP);

    fprintf(F," %18.11E %18.11E %18.11E %3d\n",
             Finfo->cs,0.,1., LPRUP); 
    fprintf(F,"</init>\n");

 
  for(N=0;N<Finfo->nEvents ;N++)
  { int Nmom,CC,i;
    double mom[40],Qf,alphaQCD;
    int clr1[10],clr2[10],icol[10][2];
    int w;

    if(readEvent(Finfo, &Nmom, mom, clr1, clr2, &Qf,&alphaQCD, &w)) break;
    fprintf(F,"<event>\n");
    
    fprintf(F,"%2d %4d %15.7E %15.7E %15.7E %15.7E\n",
          Finfo->Nin+Finfo->Nout,  //  NUP
           1,                      //  IDPRUP
           1.,                     //  XWGTUP,
           Qf,                     //  SCALUP,         
          -1.,                     //  AQEDUP,
           alphaQCD                //  AQCDUP
          );
          
    for(II=0;II<Finfo->Nin+Finfo->Nout;II++) for(J=0;J<2;J++) icol[II][J]=0;
    for(i=0,CC=500;clr1[i];i++,CC++)
    { int k1=clr1[i]-1,k2=clr2[i]-1; 
      if(k1<2) icol[k1][0]=CC; else icol[k1][1]=CC;
      if(k2<2) icol[k2][1]=CC; else icol[k2][0]=CC;
    }

          
          
    for(II=0;II< Finfo->Nin+Finfo->Nout;II++)
    { double s, P[4]={0,0,0,0};
      int m1,m2;
      fprintf(F," %8d %4d",  Finfo->PIDs[II] , II<Finfo->Nin? -1:1);
       
      if(II< Finfo->Nin) {m1=0;m2=0;} else { if(Finfo->Nin==2) {m1=1;m2=2;} else {m1=1;m2=1;}} 
     
      fprintf(F," %4d  %4d", m1,m2);
      
      
      for(J=0;J<2;J++) fprintf(F," %4d", icol[II][J]);
      if(II<Finfo->Nin) { if(Finfo->Nin==2) P[2]=mom[II];}
      else for(J=0;J<3;J++) P[J]=mom[3*II-2*Finfo->Nin+J];
      for(J=0;J<3;J++) P[3]+=P[J]*P[J];
      P[3]=sqrt(P[3]+Finfo->pmass[II]*Finfo->pmass[II]); 
      for(J=0;J<4;J++) fprintf(F," %18.11E",P[J] );
      fprintf(F," %18.11E", Finfo->pmass[II]);
      fprintf(F," %11.4E %4.1f\n",0.,9.);
#ifdef NotDone        
      fprintf(F," %8d %4d",  E_.IDUP[II] ,E_.ISTUP[II]);
      for(J=0;J<2;J++) fprintf(F," %4d", E_.MOTHUP[II][J]);
      for(J=0;J<2;J++) fprintf(F," %4d", E_.ICOLUP[II][J]);
      for(J=0;J<4;J++) fprintf(F," %18.11E",E_.PUP[II][J]);
      if(E_.ISTUP[II]==2)
      { 
          s= E_.PUP[II][3]*E_.PUP[II][3] 
          - E_.PUP[II][0]*E_.PUP[II][0] -E_.PUP[II][1]*E_.PUP[II][1]-E_.PUP[II][2]*E_.PUP[II][2];
          if(s<0) s=-sqrt(-s); else s=sqrt(s); 
      } else s=E_.PUP[II][4]; 
      fprintf(F," %18.11E", s);
      fprintf(F," %11.4E %4.1f\n",E_.VTIMUP[II],E_.SPINUP[II]);
#endif

    }	
    fprintf(F,"</event>\n");  
  }       

  fprintf(F,"</LesHouchesEvents>\n");
  if(F!=stdout)    fclose(F);     
  return 0;
}
