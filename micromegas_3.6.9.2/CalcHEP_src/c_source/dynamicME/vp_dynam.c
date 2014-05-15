#include<math.h>
#include<stdio.h>
#include<unistd.h>
#include<sys/stat.h>
#include<sys/types.h>

#include <dlfcn.h>
#include <sys/wait.h> 
             

#include"num_in.h"
#include"num_out.h"
#include"VandP.h"
#include"dynamic_cs.h"
#include"rootDir.h"

static void * VandP=NULL;

int nModelParticles;
ModelPrtclsStr*ModelPrtcls;
int nModelVars;
int nModelFunc;
char**varNames;
REAL *varValues;
static  int (*calcMainFunc_)(void);


int calcMainFunc(void) { return calcMainFunc_();} 

int getDynamicVP(void)
{  int err;
   char *SO=malloc(strlen(compDir)+30);
   sprintf(SO,"%s/so_generated/VandP.so",compDir); 
   if(access(SO,X_OK&R_OK)|| checkMtime(SO) )  // [re]compile VandP.so 
   {
     char *command=malloc(300+strlen(rootDir)+strlen(compDir) ); 
     sprintf(command,"CALCHEP=%s;"
                     "cd %s;"
                     "$CALCHEP/bin/make_VandP models 1 6\n"
                     " . $CALCHEP/FlagsForSh;"
                     " . ./EXTLIBsh;"
                     "$CC $CFLAGS $SHARED -o so_generated/VandP.so VandP.c $CALCHEP/include/VandPgate.c $CALCHEP/lib/dummy.a $EXTLIB $CALCHEP/lib/libSLHAplus.a -lm"
                    ,rootDir,compDir);
     err=system(command);
     free(command);
     if(WIFSIGNALED(err) ||WEXITSTATUS(err) )
     { 
        printf("Can not compile Constraints\n"); { free(SO);return 1;}
     }
     if(VandP) dlclose(VandP); VandP=NULL;
   }
   
   if(VandP==NULL) VandP=dlopen(SO,RTLD_NOW);
   if(!VandP) {  printf("dlopen problem with %s\n",SO);
   messanykey(10,10,dlerror()); free(SO); return 2;}
   { 
     int * nModelVars_;
     int * nModelFunc_;
     int *nModelParticles_;
     char ***varNames_; 
     REAL **varValues_; 
     ModelPrtclsStr **ModelPrtcls_;
   
     nModelVars_     =dlsym(VandP,"nModelVars");
     nModelFunc_     =dlsym(VandP,"nModelFunc");
     nModelParticles_=dlsym(VandP,"nModelParticles" );
     varNames_       =dlsym(VandP,"varNames");
     varValues_      =dlsym(VandP,"varValues");
     ModelPrtcls_    =dlsym(VandP,"ModelPrtcls");
     calcMainFunc_   =dlsym(VandP,"calcMainFunc");

     if(!nModelParticles_||!ModelPrtcls_||!nModelVars_||!nModelFunc_
      ||!varNames_ ||!varValues_||!calcMainFunc_) 
     { messanykey(10,10,"Can not find symbols"); free(SO); return 3;}
     nModelVars=*nModelVars_;
     nModelFunc=*nModelFunc_;     
     nModelParticles=*nModelParticles_;
     varNames=*varNames_;
     ModelPrtcls=*ModelPrtcls_;
     varValues=*varValues_;
   }
   free(SO);   
   return 0;
}


int setModel(char * modelDisp , int nModel )
{  int err,newD;
   int size=100;
   struct stat buf;
   for(;;)    
   {  compDir=realloc(compDir,size+20);
      if(getcwd(compDir,size)) break; else size*=2;
   }
   strcat(compDir,"/aux");
   libDir=malloc(strlen(compDir)+20);
   sprintf(libDir,"%s/so_generated",compDir);
   modelNum=nModel;
   calchepDir= rootDir;
   if(modelDisp[0]=='/')
   {  modelDir=realloc(modelDir,strlen(modelDisp)+1);
      strcpy(modelDir,modelDisp);
   }
   else
   { modelDir=realloc(modelDir,size+10+strlen(modelDisp));
     getcwd(modelDir,size);
     strcat(modelDir,"/");
     strcat(modelDir,modelDisp);
   }
   newD=prepareWorkPlace();
   if(newD) { if(stat(libDir,&buf)) mkdir(libDir,00755);}
   else 
   { if(checkWorkPlace())
     { char*command=malloc(strlen(libDir)+20);
       sprintf(command,"rm -f %s/*.so",libDir);
       system(command);
       free(command);
     }  
   }
   delAllLib();
   if(getDynamicVP()) return 2; 
   cleanDecayTable();
   return 0;
}


int assignVal(char * name, double val)
{
  int i; 
  for(i=0;i<nModelVars;i++)
  { 
    if(strcmp(name,varNames[i])) continue;
    varValues[i]=val; return 0;
  }
  return 2;
}

int assignValW(char*name, double  val)
{
  if(assignVal(name,val)) { printf(" %s not found\n",  name); return 1;}
  else return 0;
}
  


int readVar(char *fname)
{
  double val;   
  char name[80];
  int n;
  FILE * f=fopen(fname,"r");
  if(f==NULL) return -1;
  
  for(n=1;;n++)
  { if(fscanf(f,"%s",name)!=1) { n=0; break;}
    if(name[0]=='#') { fscanf(f,"%*[^\n]"); continue;}
    if(fscanf(f,"%lf",&val)!=1) break;
    fscanf(f,"%*[^\n]");
    { int err=assignVal(name,val);
      if(err==1) break;
    }
  }
  fclose(f);
  return n;                                         
}
 