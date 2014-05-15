#include"../../sources/micromegas.h"
#include"../../sources/micromegas_aux.h"
#include"pmodel.h"
#include"lpath.h"

#include<sys/wait.h>
#include<unistd.h>

#define FIN  "nngg.in"
#define FOUT "nngg.out"

int loopGamma(double * csAA, double *csAZ)
{
  double sigmav;
  char buff[2000];
  int err;
  FILE*f;
  int GI=0;
   
  *csAA=0,*csAZ=0; 

  if(!access(FOUT,R_OK)) unlink(FOUT);
  
  sprintf(buff, LPATH "/nngg/lGamma.exe");
  if(access( buff,X_OK))
  { char buf[2000]; 
    sprintf(buf, "make -C " LPATH "/nngg");
    system(buf);
  } 
  if(access( buff,X_OK)) 
  {  
    printf("Can not found/compile executable %s\n",buff);
    return 10;
  }  

  if(GI)  sprintf(buff+strlen(buff)," GI ");
  err=System(buff);   
  
  if(err>=0) 
  {  err=slhaRead(FOUT,1);
     if(err) return err;
     *csAZ=slhaVal("Lgamma",0.,1,1)*2.9979E-26;
     *csAA=slhaVal("Lgamma",0.,1,2)*2.9979E-26;
  }  

//  if(!access(FOUT,R_OK)) unlink(FOUT);
//  if(!access(FIN,R_OK)) unlink(FIN);
  return err;
}  

extern int  loopgamma_(double * cs1, double *cs2); /* fortran */
int  loopgamma_(double * cs1, double *cs2){ return loopGamma(cs1,cs2); }
 