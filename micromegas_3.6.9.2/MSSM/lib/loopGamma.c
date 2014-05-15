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

  err=slhaWrite(FIN);
  if(err) return err; 
  f=fopen(FIN,"a");
  if(slhaDecayExists(36)<0) slhaDecayPrint("H3",f);
  if(slhaDecayExists(25)<0) slhaDecayPrint("h",f);  
  if(slhaDecayExists(35)<0) slhaDecayPrint("H",f);  
  fclose(f);    
  if(!access(FOUT,R_OK)) unlink(FOUT);
  
  sprintf(buff+strlen(buff)," %s %s",FIN,FOUT);
  err=System(buff);   
  
  if(err>=0) 
  {  err=slhaRead(FOUT,1);
     *csAZ=slhaVal("Lgamma",0.,1,1)*2.9979E-26;
     *csAA=slhaVal("Lgamma",0.,1,2)*2.9979E-26;
  }  

//  if(!access(FOUT,R_OK)) unlink(FOUT);
//  if(!access(FIN,R_OK)) unlink(FIN);
  return err;
}  

extern int  loopgamma_(double * cs1, double *cs2); /* fortran */
int         loopgamma_(double * cs1, double *cs2){ return loopGamma(cs1,cs2); }
 