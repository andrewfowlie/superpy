#include <stdlib.h>


#ifdef __hpux
#include<dl.h>
#else
#include <dlfcn.h>
#endif

#include <unistd.h>
#include <ctype.h>
#include <string.h>

#include <sys/utsname.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


#include "dynamic_cs.h"

#include "../../include/V_and_P.h"

REAL * varAddress(char *name)
{int i;
 for(i=0;i<nModelVars+nModelFunc;i++) if(!strcmp(name,varNames[i]))return varValues+i;
 return NULL;
}
  
char  * libDir=NULL;
char  * modelDir=NULL;
char  * compDir=NULL;
char  * calchepDir=NULL;
int   modelNum=0;


int  prepareWorkPlace(void)
{  char * command;
   struct stat buf;
   int err,len;

   if(!compDir) return  -1;   

   err=stat(compDir,&buf);

   len=strlen(compDir)+500;
   if(err && modelDir) len+=strlen(modelDir);
   command=malloc(len);  
   if(err) 
   {  char * dir[3]={"tmp","results","models"};
      int i;
      if(mkdir(compDir, 00755)) return -3; 
      for(i=0;i<3;i++) 
      { 
         sprintf(command,"%s/%s",compDir,dir[i]);
         mkdir(command,00755);
      } 
     if(modelDir && modelNum)
     {
      sprintf(command, "dName=%s\n"          
       "for FILE in vars func prtcls lgrng extlib\n do\n" 
       "  if(test ! -r $dName/models/\"$FILE\"1.mdl) then\n" 
       "   cp %s/\"$FILE\"%d.mdl $dName/models/\"$FILE\"1.mdl\n"
       "  fi\n"
       "done\n",compDir, modelDir,modelNum);
       system(command);
     }      
   }
   else
   {   
     sprintf(command, "dName=%s\n"
     "for FILE in $dName/tmp/* $dName/results/*\n"
     "do\n"
     " if(test ! -d $FILE) then\n"
     "   rm -f $FILE\n"
     " fi\n" 
     "done\n",compDir);     
     system(command);
   }
   free(command);
   return err;
}

/*#define TEST*/
int cleanWorkPlace(void)
{ 
#ifndef TEST
   char * command;
   command=malloc(strlen(compDir)+100);
   sprintf(command,"rm -fr %s",compDir);
   system(command);
   free(command);
#endif    
   return 0; 
}

static int  checkMtime(char * fname)
{ int i,L;
  time_t tt;
  struct stat buff;
  char * mf[4]={"vars","func","prtcls","lgrng"};
  char *mfname;
  if(modelDir==NULL) return 0;
  stat(fname,&buff); tt=buff.st_mtime;
  L=strlen(modelDir)+20;
  mfname=malloc(strlen(modelDir)+20); 
  
  for(i=0;i<4;i++)
  { sprintf(mfname,"%s/%s%d.mdl",modelDir,mf[i],modelNum);
    stat(mfname,&buff);
    if(buff.st_mtime > tt) break;
  }
  free(mfname);
  if(i<4) {unlink(fname); return 1;}
  return 0;
}

static void* newSymbol(void*handle,char *name)
{
#ifdef __hpux
void * addr;
     if(shl_findsym((shl_t*)&handle,name,TYPE_UNDEFINED,&addr)) return NULL;
       else return addr;
#else
      return dlsym(handle, name);
#endif
}

static void * dLoad(char * libName)
{
void *q;

if(access(libName,R_OK)) return NULL;

#ifdef __hpux
   return  shl_load(libName,0,0L);
#else
   q= dlopen(libName, RTLD_NOW);
   if(!q) printf("%s\n",dlerror()); 
   return q;
#endif
}


static void dClose(void * handle)
{
#ifdef __hpux
       shl_unload(handle);
#else
       dlclose(handle);
#endif
}

extern double aWidth(char *);

static numout* loadLib(void* handle, char * lib)
{ numout * cc=malloc(sizeof(numout));
  char name[100];
  if(!handle) {free(cc); return NULL;}   
  cc->handle=handle;
  sprintf(name,"interface_%s",lib);
  cc->interface=newSymbol(handle, name);
  if(!cc->interface || cc->interface->nprc==0){free(cc); return NULL;}
  else
  {  int i;
     cc->init=0;
     cc->Q=NULL, cc->SC=NULL;
     cc->link=malloc(sizeof(double*)*(1+cc->interface->nvar));
     cc->link[0]=NULL;
     for(i=1;i<=cc->interface->nvar;i++) 
     { char *name=cc->interface->varName[i];
       cc->link[i]=varAddress(name);
if(cc->link==NULL) printf("No link for %s\n",name);       
       if(strcmp(name,"Q")==0) cc->Q=cc->interface->va+i; 
       else if(strcmp(name,"SC")==0) cc->SC=cc->interface->va+i;
     }  
     *(cc->interface->aWidth)=&aWidth;
  }
  return cc;
}


typedef struct  procRec 
{ struct procRec  * next;
  char * libname;
  numout * cc;
}  procRec;   

static  procRec* allProc=NULL;

#include<stdlib.h>
#include<sys/wait.h>

numout*getMEcode(int twidth,int Gauge, char*Process, char*excludeVirtual, 
                          char*excludeOut,char*lib)
{
   char *proclibf,*command;
   void * handle=NULL;
   int new=0;
   numout * cc;
   procRec*test;
   int Len;
   char * lib_;

   lib_=malloc(strlen(lib)+4);
   
   if(Gauge) sprintf(lib_,"%s_u",lib);    else  strcpy(lib_,lib); 
     
   for(test=allProc;test; test=test->next)
   { if(strcmp(lib_,test->libname)==0) return test->cc;}
   
   Len=strlen(compDir)+strlen(lib)+strlen(libDir)+300;
   proclibf=malloc(Len);

   if(Process) Len+=strlen(Process);
   if(excludeVirtual) Len+=strlen(excludeVirtual);
   if(excludeOut) Len+=strlen(excludeOut);
   command=malloc(Len);
 

   sprintf(proclibf,"%s/%s.so",libDir,lib_);

   if(access(proclibf,R_OK)==0 && checkMtime(proclibf)==0) handle=dLoad(proclibf);
   if(!handle)
   {  int i;
   
      for(i=0;Process[i]==' ';i++); if(Process[i]==0)
      {    
        free(command); free(proclibf); free(lib_); 
        return NULL;    
      }      
      if(!handle)
      {
        char options[20];
        char GaugeCh[4];
        int ret;  
        int delWorkDir;
        
        if(twidth) strcpy(options,"5[[{[{}");else strcpy(options,"");
        if(Gauge) strcpy(GaugeCh,"U"); else strcpy(GaugeCh,"F");
 
        delWorkDir=prepareWorkPlace();      

        sprintf(command,"cd %s; %s/bin/newProcess %s %s \"%s\" %s \"%s\"",
                       compDir, calchepDir, lib_, libDir,options,GaugeCh,Process);
  
        if(excludeVirtual) sprintf(command+strlen(command)," \"%s\"",excludeVirtual);
        else  sprintf(command+strlen(command)," \"\"");            
        if(excludeOut) sprintf(command+strlen(command)," \"%s\"",excludeOut);       
        ret=system(command);
      
        if(ret<0 || WIFSIGNALED(ret)>0 ) exit(10);
        if(delWorkDir )cleanWorkPlace();
        
        if(ret==0) handle=dLoad(proclibf); else 
        { printf(" Can not compile %s \n", Process);
          free(command); free(proclibf); free(lib_);
          return NULL;
        } 
        if(!handle)
        { printf(" Can not load the compiled library %s \n",proclibf);
           free(command); free(proclibf); free(lib_);
          return NULL;
        }         
        new=1;   
      }
   }
   cc=loadLib(handle,lib_);
   if(!cc && new) dClose(handle);
   if(cc)
   {  test=(procRec*)malloc(sizeof(procRec));
      test->next=allProc; allProc=test;
      test->libname=(char*) malloc(strlen(lib_)+1);
      strcpy(test->libname,lib_);
      test->cc=cc;
   } else if(new) dClose(handle);  
    free(command); free(proclibf); free(lib_);
    return cc; 
}

txtList  makeProcList(char ** InNames, char** OutNames, int nx)
{ 
  char lname[20],buff[200];
  char * fnameG;
  FILE *f;
  txtList List=NULL;
  char process[50],lib[50];
  char ** c;
  
  process[0]=0;
  sprintf(lib,"pList_");
  
  for(c=InNames;*c;c++)
  {  if(c!=InNames)strcat(process,",");
/*printf("%s\n",*c);  */
     strcat(process,*c);
     pname2lib(*c,lname);
     strcat(lib,lname);
  }
  strcat(process,"->");
  strcat(lib,"_");
  for(c=OutNames;*c;c++)
  {  if(c!=OutNames)strcat(process,",");
     strcat(process,*c);
     pname2lib(*c,lname);
     strcat(lib,lname);     
  }

  if(nx)
  { if(OutNames[0]) strcat(process,",");   
    sprintf(process+strlen(process),"%d*x",nx);
    sprintf(lib+strlen(lib),"_%dx",nx);
  }   
/*  
printf("lib=%s\n",lib);
printf("process=%s\n",process);
*/
  fnameG=malloc(strlen(libDir)+50);
  sprintf(fnameG,"%s/%s",libDir,lib);
    
  if(access(fnameG, R_OK)|| checkMtime(fnameG))
  {  char * command=malloc(strlen(compDir) + strlen(calchepDir) + strlen(libDir) + 200);
     int delWorkDir=prepareWorkPlace();
     sprintf(command,"cd %s;"
       " %s/bin/s_calchep -blind \"{{%s{{{[[{0\";"
       " mv results/list_prc.txt %s/%s",
       compDir,calchepDir,process, libDir,lib);
     system(command);
     free(command);
     if(delWorkDir) cleanWorkPlace();
  }
  f=fopen(fnameG,"r");
  free(fnameG);
  if(!f) return NULL;
  for(; 1==fscanf(f,"%[^\n]\n",buff); )
  { txtList l=malloc(sizeof(txtListStr));
    l->next=List; List=l;
    l->txt=malloc(strlen(buff)+1);
    strcpy(l->txt,buff);
  }
  fclose(f);
  return List;
}


txtList  makeDecayList(char * pname, int nx)
{ 
  char fnameL[50],lname[20],buff[200];
  char * fnameG;
  FILE *f;
  txtList List=NULL;

printf("asks 2->%d for %s\n",nx,pname);
  pname2lib(pname,lname);
  sprintf(fnameL,"dList_%s_%dx",lname,nx);
  fnameG=malloc(strlen(libDir)+50);
  sprintf(fnameG,"%s/%s",libDir,fnameL);
  if(access(fnameG, R_OK)|| checkMtime(fnameG))
  {  char * command=malloc(strlen(compDir) + strlen(calchepDir)+strlen(libDir) +200);
     int delWorkDir=prepareWorkPlace();
     sprintf(command,"cd %s;"
       " %s/bin/s_calchep -blind \"{{%s->%d*x{{{[[{0\";"
       " mv results/list_prc.txt %s/%s",
       compDir,calchepDir,pname,nx,libDir,fnameL);
     system(command);
     free(command);
     if(delWorkDir) cleanWorkPlace();
  }

  f=fopen(fnameG,"r");
  free(fnameG);
  if(!f) return NULL;
  for(; 1==fscanf(f,"%[^\n]\n",buff); )
  { txtList l=malloc(sizeof(txtListStr));
    l->next=List; List=l;
    l->txt=malloc(strlen(buff)+1);
    strcpy(l->txt,buff);
  }
  fclose(f);
  return List;
}


void cleanTxtList(txtList L)
{ 
  while(L) {txtList l=L; free(L->txt); l=L;  L=L->next; free(l);}  
}

void printTxtList(txtList L, FILE *f)
{ for(;L;L=L->next) fprintf(f,"%s\n",L->txt);}


int qNumbers(char*pname, int *spin2, int * charge3, int * cdim)
{
  int n=pTabPos(pname);
  int sign=1,pdg;
  if(!n) { if(pname[0]=='#' && sscanf(pname+1,"%d",&pdg)==1) return pdg; else return 0;}
  if(n<0){ n=-n; sign=-1;}
  if(spin2)   *spin2  =ModelPrtcls[n-1].spin2;
  if(charge3) *charge3=sign*ModelPrtcls[n-1].q3;
  pdg=sign*ModelPrtcls[n-1].NPDG;
  if(cdim)    
  { *cdim   =ModelPrtcls[n-1].cdim; 
    if(sign==-1 &&(*cdim==3 || *cdim==-3)) (*cdim)*=-1;
  }
  return pdg;
}

void delAllLib(void)
{
  procRec* curProc=allProc;
  while(curProc)
  {  procRec*tmp=curProc;
     free(curProc->libname);
     free(curProc->cc->link);
     dClose(curProc->cc->handle);     
     free(curProc->cc);
     curProc=curProc->next;
     free(tmp);
  }
  allProc=NULL;
}          