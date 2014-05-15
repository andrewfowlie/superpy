#include"suspect_path.h"

/*===================================*/
#include"pmodel.h"
#include"pmodel_aux.h"
#include"../../sources/micromegas_aux.h"
#include<sys/wait.h>
#include<unistd.h>

#define FIN  "suspect2_lha.in"
#define FOUT "suspect2_lha.out"

static int SystemCall(int mode)
{ 
  char buff[2000];
  int err;

  if(!access(FOUT,R_OK)) unlink(FOUT);

  sprintf(buff,"%s/suspect.exe",SUSPECT);
  if(access( buff,X_OK))
  { printf("Executable \n %s\n is not found. Program stops.\n",buff);
    exit(13);
  }  
  
  err=System(buff);   
  if(err>=0)  err=slhaRead(FOUT,4); else cleanSLHAdata();   
  return err;
}

static void addSU_ALGO(char * fname, int * ichoice)
{
   FILE * f=fopen(fname,"a");
   fprintf(f,"Block SU_ALGO # ichoice \n");
   fprintf(f,
"   2    %d  # 2-loop RGE \n"
"   3    %d   # 1: g_1(gut) = g_2(gut) consistently calculated from input\n"
"   4    %d   # RGE accuracy: accurate, but slow\n"
"   6    %d   # MA_pole, MU(EWSB) input\n"
"   7    %d   #  choice for sparticles masses rad. corr. (=/= h)\n"
"   8    %d   # 1 (defaut): EWSB scale=(mt_L*mt_R)^(1/2)\n"
"   9    %d   # Final spectrum accuracy: 1 -> 1%% acc.; 2 -> 0.01 %% acc.\n"
"   10   %d   # Higgs boson masses: One-loop  + dominant DSVZ 2-loop\n"
"   11   %d   # RUNNING DRbar Higgs masses at loop-level at mZ\n",
    ichoice[2],ichoice[3],ichoice[4],ichoice[6],ichoice[7],ichoice[8],ichoice[9],
    ichoice[10],ichoice[11]);
   
    fclose(f);
}


double  suspectSUGRAc(double tb, double gMG1,double gMG2,double gMG3,
             double gAl, double gAt, double gAb, double sgn, double gMHu, double gMHd,
             double gMl2,double gMl3,double gMr2,double gMr3,
             double gMq2,double gMq3,double gMu2,double gMu3,double gMd2,double gMd3)
{ 
   int ichoice[12]={0,0,21, 1, 2, 0, 1, 2, 1, 2, 2,  0};
   
   if(sugraLesH(FIN,tb, gMG1,gMG2,gMG3,
             gAl, gAt, gAb, sgn, gMHu, gMHd,
             gMl2,gMl3,gMr2,gMr3,
             gMq2,gMq3,gMu2,gMu3,gMd2,gMd3))
   {  printf("can not write down %s file\n",FIN); exit(10);}
   addSU_ALGO(FIN,ichoice);
   return SystemCall(1);
}

double  suspectSUGRAnuhc(double tb, double gMG1,double gMG2,double gMG3,
             double gAl, double gAt, double gAb, double gMl2,double gMl3,double gMr2,double gMr3,
             double gMq2,double gMq3,double gMu2,double gMu3,double gMd2,double gMd3,double mu,double MA)
{ 
   int ichoice[12]={0,0,21, 1, 2, 0, 0, 2, 1, 2, 2,  0};
   
   if(sugraHiggsLesH(FIN,tb, gMG1,gMG2,gMG3,
             gAl, gAt, gAb, gMl2,gMl3,gMr2,gMr3,
             gMq2,gMq3,gMu2,gMu3,gMd2,gMd3,mu,MA))
   {  printf("can not write down %s file\n",FIN); exit(10);}
   addSU_ALGO(FIN,ichoice);
   return SystemCall(1);
}




double  suspectAMSBc(double m0,double m32, double tb, double sgn)
{  int ichoice[12]={0,0,21, 1, 2, 0, 1, 2, 1, 2, 2,  0};
  
   if(amsbLesH(FIN,m0,m32,tb, (int)sgn))
   { printf("can not write down %s file\n",FIN); exit(10);}
   addSU_ALGO(FIN,ichoice);    
   return SystemCall(1);
}



double suspectEwsbMSSMc(double tb, double MG1, double MG2, double MG3, double Al, double At, double Ab, 
  double mu, double MH3, double Ml1, double Ml2, double Ml3, double Mr1, double Mr2, double Mr3, 
  double Mq1, double Mq2, double Mq3, double Mu1, double Mu2, double Mu3, 
                                    double Md1, double Md2, double Md3)
{  int ichoice[12]={0,0,21, 1, 2, 0, 0, 2, 1, 2, 2,  0};
   int err;
   
   if(EWSBLesH(FIN,tb,MG1,MG2,MG3,Al,At,Ab,mu,MH3,Ml1,Ml2,Ml3,Mr1,Mr2,Mr3, Mq1,Mq2,Mq3,Mu1,Mu2,Mu3, 
                                     Md1,Md2,Md3)
     ){printf("can not write down %s file\n",FIN); exit(10);}
   addSU_ALGO(FIN,ichoice);
   err= SystemCall(0);
 
   return err;
}


