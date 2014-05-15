#include"../../sources/micromegas.h"
#include"../../sources/micromegas_aux.h"
#include"pmodel.h"


#include<sys/wait.h>
#include<unistd.h>

#define FIN  "nngg.in"
#define FOUT "nngg.out"

int loopGamma(double * csAZ, double *csAA)
{
  double sigmav;
  char buff[2000];
  int err;
  FILE*f;
   
  *csAA=0,*csAZ=0; 

  if(!access(FOUT,R_OK)) unlink(FOUT);
  
  sprintf(buff,"%s/../lib/nngg/lGamma.exe",WORK);
  if(access( buff,X_OK))
  { char buf[2000]; 
    sprintf(buf, "make -C %s/../lib/nngg",WORK);
    system(buf);
  } 
  if(access( buff,X_OK)) 
  {  
    printf("Can not found/compile executable %s\n",buff);
    return 10;
  }  

  f=fopen(FIN,"w");
  if(!f) return 1; 
  fprintf(f, "BLOCK Z3_parameters\n");
  
  fprintf(f, " EE    %f\n",findValW("EE"));
  fprintf(f, " SW    %f\n",findValW("SW"));
  fprintf(f, " MZ    %f\n",findValW("MZ"));
  fprintf(f, " Q     %f\n",findValW("Q"));
  fprintf(f, " wZ    %f\n",findValW("wZ"));
  fprintf(f, " wW    %f\n",findValW("wW"));
  fprintf(f, " Mm    %f\n",findValW("Mm"));
  fprintf(f, " Ml    %f\n",findValW("Ml"));
  fprintf(f, " Mu    %f\n",findValW("Mu"));
  fprintf(f, " Md    %f\n",findValW("Md"));
  fprintf(f, " Mc    %f\n",findValW("Mc"));
  fprintf(f, " Ms    %f\n",findValW("Ms"));
  fprintf(f, " Mtop  %f\n",findValW("Mtop"));
  fprintf(f, " wtop  %f\n",findValW("wtop"));
  fprintf(f, " Mh    %f\n",findValW("Mh"));
  fprintf(f, " la3   %f\n",findValW("la3"));
  fprintf(f, " la2   %f\n",findValW("la2"));
  fprintf(f, " la4   %f\n",findValW("la4"));
  fprintf(f, " laS   %f\n",findValW("laS"));
  fprintf(f, " laS1  %f\n",findValW("laS1"));
  fprintf(f, " laS2  %f\n",findValW("laS2"));
  fprintf(f, " laS21 %f\n",findValW("laS21"));
  fprintf(f, " Mdm1  %f\n",findValW("Mdm1"));
  fprintf(f, " Mdm2  %f\n",findValW("Mdm2"));
  fprintf(f, " sinDm %f\n",findValW("sinDm"));
  fprintf(f, " muppS %f\n",findValW("muppS"));
  fprintf(f, " Mcp   %f\n",findValW("Mcp"));
  fclose(f);
 
  if(!access(FOUT,R_OK)) unlink(FOUT);
  
  sprintf(buff+strlen(buff)," %s %s",FIN,FOUT);
  err=System(buff);   
  
  if(err>=0) 
  {  err=slhaRead(FOUT,1);
     *csAZ=0;
     *csAA=slhaVal("Lgamma",0.,1,2)*2.9979E-26;
  }  

//  if(!access(FOUT,R_OK)) unlink(FOUT);
//  if(!access(FIN,R_OK)) unlink(FIN);
  return err;
}  

extern int  loopgamma_(double * cs1, double *cs2); /* fortran */
int         loopgamma_(double * cs1, double *cs2){ return loopGamma(cs1,cs2); }
 