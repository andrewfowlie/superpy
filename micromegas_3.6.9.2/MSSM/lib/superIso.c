#include "../../sources/micromegas.h"
#include "../../sources/micromegas_aux.h"
#include "pmodel.h"

int callSuperIsoSLHA(void)
{  char* command;
   int err;
   char * path=getenv("superIso");
   if(!path) { printf("Define environment variable 'superIso' to specify path to the package\n"); return -1;}
   command=malloc(strlen(path)  +100);
   sprintf(command,"%s/slha.x",path); 
   if(access(command,X_OK))
   { int sysTimeLimTmp=sysTimeLim;
     sysTimeLim=0;
     System("cd %s; make; make slha",path);
     sysTimeLim=sysTimeLimTmp;
   }

   if(access(command,X_OK)) 
   { printf("Can not compile superIso\n");
     free(command);
     return -2;
   }
   slhaWrite("slhaForSuperIso");
   System("%s  slhaForSuperIso >/dev/null  ",command);
   err=slhaRead("output.flha",1);
   unlink("slhaForSuperIso");
   free(command);
   return err;
}

int callsuperisoslha_(void){ return callSuperIsoSLHA();}