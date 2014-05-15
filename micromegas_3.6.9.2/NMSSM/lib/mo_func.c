#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include"../../sources/micromegas.h"
#include"../../sources/micromegas_aux.h"
#include"pmodel.h"
#include"pmodel_aux.h"

int readVarNMSSM(char * fname)
{
  char*vlist[32]={"alfSMZ","McMc","MbMb","Mtp","tb","MG1","MG2","MG3",
  "Ml2","Ml3","Mr2","Mr3","Mq2","Mq3","Mu2","Mu3","Md2","Md3","At","Ab","Al",
  "Au","Ad","Am","mu","Lambda","Kappa","aLambda","aKappa","MZ","Ml","wt"};
  return readVarSpecial(fname,32,vlist);
} 

int readVar_nMSSM(char * fname)
{
  char*vlist[34]={"alfSMZ","McMc","MbMb","Mtp","tb","MG1","MG2","MG3",
  "Ml2","Ml3","Mr2","Mr3","Mq2","Mq3","Mu2","Mu3","Md2","Md3","At","Ab","Al",
  "Au","Ad","Am","mu","Lambda","Kappa","aLambda","aKappa","MZ","Ml","wt","MA","MP"};
  return readVarSpecial(fname,34,vlist);
} 


void o1Contents(FILE * f)
{ double val; 
  int err; 
 
  fprintf(f,"\n~o1 = ");
  err=findVal("Zn11",&val);
  if(err==0) fprintf(f,"%.3f*",val); else  fprintf(f,"???*");  
  fprintf(f,"bino");

  err=findVal("Zn12",&val);
  if(err==0) fprintf(f,"%+.3f*",val); else  fprintf(f,"???*");  
  fprintf(f,"wino");
  
  err=findVal("Zn13",&val);
  if(err==0) fprintf(f,"%+.3f*",val); else  fprintf(f,"???*");  
  fprintf(f,"higgsino1");
  
  err=findVal("Zn14",&val);
  if(err==0) fprintf(f,"%+.3f*",val); else  fprintf(f,"???*");  
  fprintf(f,"higgsino2");
  
  err=findVal("Zn15",&val);
  if(err==0) fprintf(f,"%+.3f*",val); else  fprintf(f,"???*");  
  fprintf(f,"singlino\n");
}


static void FillVal(int mode)
{ char name[10];
  int i,j,k;
  double Pa[2][3]={{0,0,0},{0,0,0}};
  double Q;
  char format[20];
  char fmt[20];     

  char * softName[13]={"MG1","MG2","MG3","Ml2","Ml3","Mr2","Mr3","Mq2","Mq3","Mu2","Mu3","Md2","Md3"};
  int      softId[13]={   1,    2,    3,   32 ,  33 ,  35 ,  36 ,  42 ,  43 ,  45 ,  46 ,  48 ,  49};

  char * nmssmName[4]={"Lambda","Kappa", "aLambda","aKappa"};
  int    nmssmId[4]   ={    1   ,   2   ,   3     ,  4     };

  
  { // Masses
    char* massName[35]={"Mha","Mhb","Mh1","Mh2","Mh3","MHc", "MNE1", "MNE2", "MNE3", "MNE4", "MNE5",    "MC1",  "MC2",  "MSG", "MSne", "MSnm", "MSnl", "MSeL", "MSeR", "MSmL", "MSmR", "MSl1", "MSl2", "MSdL", "MSdR", "MSuL", "MSuR", "MSsL", "MSsR", "MScL", "MScR", "MSb1", "MSb2", "MSt1","MSt2"};
    int   massId[35]  ={ 36  ,  46,  25,   35,   45,     37,1000022,1000023,1000025,1000035, 1000045, 1000024,1000037,1000021,1000012,1000014,1000016,1000011,2000011,1000013,2000013,1000015,2000015,1000001,2000001,1000002,2000002,1000003,2000003,1000004,2000004,1000005,2000005,1000006,2000006};   
    for(i=0;i<35;i++) assignValW(massName[i],slhaVal("MASS",0.,1,massId[i]));
    Q=sqrt(fabs(findValW("MSt1")*findValW("MSt2"))); 
  }
  { // Mixing
    char* Zf[55]={"Zb","Zt","Zl","Zu","Zv"};
    char* Qmix[5]={"SBOTMIX","STOPMIX","STAUMIX","UMIX","VMIX"};
    for(i=1;i<=5;i++) for(j=1;j<=5;j++) { sprintf(name,"Zn%d%d",i,j); assignValW(name,slhaVal("NMNMIX",Q,2,i,j));}
    for(i=1;i<=3;i++) for(j=1;j<=3;j++) { sprintf(name,"Zh%d%d",i,j); assignValW(name,slhaVal("NMHMIX",Q,2,i,j));}
    for(i=1;i<=2;i++) for(j=1;j<=3;j++) { sprintf(name,"Za%d%d",i,j); assignValW(name,slhaVal("NMAMIX",Q,2,i,j));}
    for(k=0;k<5;k++) for(i=1;i<=2;i++) for(j=1;j<=2;j++) { sprintf(name,"%s%d%d",Zf[k],i,j); assignValW(name, slhaVal(Qmix[k],Q,2,i,j));}  
  }   
  //EFFECTIVE_COUPLINGS
  for(i=1;i<=7;i++) { sprintf(name,"la%d",i);  sprintf(fmt,"L%d %%lf",i); assignValW(name,slhaValFormat("EFFECTIVE_COUPLINGS", Q,fmt));}
  for(i=1;i<=8;i++) { sprintf(name,"la%ds",i); sprintf(fmt,"K%d %%lf",i); assignValW(name,slhaValFormat("EFFECTIVE_COUPLINGS", Q,fmt));}
  for(i=1;i<=6;i++) { sprintf(name,"aa%d",i);  sprintf(fmt,"A%d %%lf",i); assignValW(name,slhaValFormat("EFFECTIVE_COUPLINGS", Q,fmt));}
  for(i=1;i<=2;i++) { sprintf(name,"B%d",i);   sprintf(fmt,"B%d %%lf",i); assignValW(name,slhaValFormat("EFFECTIVE_COUPLINGS", Q,fmt));}
  assignValW("X",   slhaValFormat("EFFECTIVE_COUPLINGS", Q, "X %lf"    ));
  assignValW("tB", slhaVal("HMIX",Q,1,2) );
  assignValW("dMb", slhaValFormat("EFFECTIVE_COUPLINGS", Q, "DELMB %lf"));

  if(mode>0)
  {  
    for(i=0;i<13;i++) assignValW(softName[i],slhaVal("MSOFT",Q,1,    softId[i])); 
    for(i=0;i<4;i++) assignValW(nmssmName[i],slhaVal("NMSSMRUN",Q,1,nmssmId[i]));
    assignValW("mu", slhaVal("HMIX",Q,1,1));
    assignValW("Al", slhaVal("Ae",Q,2,3,3));
    assignValW("Ab", slhaVal("Ad",Q,2,3,3));
    assignValW("At", slhaVal("Au",Q,2,3,3));
    assignValW("Am", slhaValExists("Ae",2,2,2)>0 ? slhaVal("Ae",Q,2,2,2):slhaVal("Ae",Q,2,3,3));
    assignValW("Au", slhaValExists("Au",2,2,2)>0 ? slhaVal("Au",0.,2,2,2):slhaVal("Au",0.,2,3,3));
    assignValW("Ad", slhaValExists("Ad",2,2,2)>0 ? slhaVal("Ad",0.,2,2,2):slhaVal("Ad",0.,2,3,3)); 
  }
  
  if(mode==2)
  { assignValW("alfSMZ",slhaVal("SMINPUTS",Q,1,3) );
    assignValW("MbMb",  slhaVal("SMINPUTS",Q,1,5) );
    assignValW("Mtp",   slhaVal("SMINPUTS",Q,1,6) );
    assignValW("Ml",    slhaVal("SMINPUTS",Q,1,7) );
  }   

}


static void CharginoZM(void)
{
  double M2,mu,g2v1,g2v2,offQ,TrX2,detX,D,tU,tV,CU,SU,CV,SV;  
  double Zn_[5][5],NMassM[5][5],M[5];
  double mc[2],Zv_[2][2],Zu_[2][2];
  char name[10];
  int i,j,k;

  for(i=1;i<=5;i++) for(j=1;j<=5;j++) 
  { sprintf(name,"Zn%d%d",i,j); Zn_[i-1][j-1]=findValW(name);}  
  for(i=1;i<=5;i++) { sprintf(name,"MNE%d",i); M[i-1]=findValW(name);} 
  for(i=0;i<5;i++) for(j=0;j<5;j++) for(k=0,NMassM[i][j]=0;k<5;k++)
  NMassM[i][j]+=Zn_[k][i]*M[k]*Zn_[k][j];

   M2=NMassM[1][1]; 
   mu=-NMassM[2][3];
   g2v1= -NMassM[1][3]*sqrt(2.);
   g2v2=  NMassM[1][2]*sqrt(2.);
        
   offQ =g2v1*g2v1+g2v2*g2v2;
   TrX2 =offQ +M2*M2+mu*mu;
   detX =mu*M2 - g2v1*g2v2;
   D=TrX2*TrX2 - 4.*detX*detX;

   tU=(g2v2*g2v2-g2v1*g2v1-M2*M2+mu*mu-sqrt(D))/2./(M2*g2v2+mu*g2v1);
   tV=(g2v1*g2v1-g2v2*g2v2-M2*M2+mu*mu-sqrt(D))/2./(M2*g2v1+mu*g2v2);

   CU=cos(atan(tU));
   SU=sin(atan(tU));
   CV=cos(atan(tV));
   SV=sin(atan(tV));
  
  Zu_[0][0]=CU;
  Zu_[0][1]=SU;
  Zu_[1][0]=-SU;
  Zu_[1][1]=CU;
  Zv_[0][0]=CV;
  Zv_[0][1]=SV;
  Zv_[1][0]=-SV;
  Zv_[1][1]=CV;

  for(i=0;i<2;i++) mc[i]=g2v1*Zu_[i][0]*Zv_[i][1]
                        +g2v2*Zu_[i][1]*Zv_[i][0]
                        +  M2*Zu_[i][0]*Zv_[i][0]
                        +  mu*Zu_[i][1]*Zv_[i][1];
 
 for(i=1;i<=2;i++)
 { sprintf(name,"MC%d",i);
   assignValW(name,mc[i-1]);
   for(j=1;j<=2;j++)
   { 
      sprintf(name,"Zu%d%d",i,j);
      assignValW(name,Zu_[i-1][j-1]);
      sprintf(name,"Zv%d%d",i,j);
      assignValW(name,Zv_[i-1][j-1]);
   }    
 }

}


#define V(N) findValW(#N)

int nmssmEWSB(int mode)
{  int err;

if(mode)
     err=ewsb_n_MSSM(V(tb),V(MG1),V(MG2),V(MG3),V(Ml2),V(Ml3),V(Mr2),V(Mr3),
     V(Mq2),V(Mq3),V(Mu2),V(Mu3),V(Md2),V(Md3),V(At),V(Ab),V(Al),V(Am),V(mu),
     V(Lambda),V(Kappa),V(aLambda),V(aKappa),V(MA),V(MP));
else err=ewsbNMSSM(V(tb),V(MG1),V(MG2),V(MG3),V(Ml2),V(Ml3),V(Mr2),V(Mr3),
     V(Mq2),V(Mq3),V(Mu2),V(Mu3),V(Md2),V(Md3),V(At),V(Ab),V(Al),V(Am),V(mu),
     V(Lambda),V(Kappa),V(aLambda),V(aKappa));


   if(err) return err;
   FillVal(0);
//   CharginoZM();   
   return err; 
}

#undef V

int nmssmSUGRA(double  m0,double mhf, double a0,double tb, double sgn,
double  Lambda, double aLambda, double aKappa)
{  int err;

   err= sugraNMSSM(m0, mhf, a0, tb, sgn, Lambda, aLambda,aKappa);
   if(err==0)
   { 
     FillVal(1);
//     CharginoZM();
   }
   return err;
}

int lesHinput(char * fname)
{  double maxl;
   int err=slhaRead(fname, 0);   
   if(err) return err;
   FillVal(2);
//   CharginoZM();
   return 0;
}


double bsgnlo_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,12);
  *P= slhaVal("LOWEN",0.,1,11);
  return slhaVal("LOWEN",0.,1,1);
}

double deltamd_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,22);
  *P= slhaVal("LOWEN",0.,1,21);
  return slhaVal("LOWEN",0.,1,2);
}

double deltams_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,32);
  *P= slhaVal("LOWEN",0.,1,31);
  return slhaVal("LOWEN",0.,1,3);
}

double bsmumu_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,42);
  *P= slhaVal("LOWEN",0.,1,41);
  return slhaVal("LOWEN",0.,1,4);
}

double btaunu_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,52);
  *P= slhaVal("LOWEN",0.,1,51);
  return slhaVal("LOWEN",0.,1,5);
}

double gmuon_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,62);
  *P= slhaVal("LOWEN",0.,1,61);
  return slhaVal("LOWEN",0.,1,6);
}
    
