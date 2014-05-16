#include"SLHAplus.h"
#define  _XOPEN_SOURCE 500
#define  _XOPEN_SOURCE_EXTENDED
#include<unistd.h>

unsigned sysTimeLim=0;
unsigned sysTimeQuant=10;

              
/*
   contains routines adapted for CalCHEP which allow to 
   prepare input SLHA file and to launch external program generation of
   for SLHA output:
             System
             openAppend
             aPrintF
*/
extern int FError;
#define  STRSIZ 2000

/* Easy extention of standard 'system' routine which allows 
   construction  of command line from severa sources 
*/ 

int System( char * format, ...) 
{  
   va_list args;
   char command[STRSIZ];
   int err;
      
   va_start(args,format);
   vsprintf(command,format,args);
   va_end(args);
   fflush(NULL);
   if(!sysTimeLim)  { 
//   printf("call system\n"); 
                      err=system(command);}
   else
   {   
     int id=fork();
     if(id==-1) err=system(command);
     else if(id)  /*parent*/
     {  long time;
        for(time=0;time<sysTimeLim;time+=sysTimeQuant)
        { usleep(sysTimeQuant*1000);
          if(waitpid(id,&err,WNOHANG)) break;
        }
        if(time>=sysTimeLim) 
         { int err,grid; 
           printf("Time limit exeeds for command:\"%s\"\n",command);
           grid=getpgid(id);   
           err=kill(-grid,SIGKILL); FError=1;
           waitpid(id,&err, WUNTRACED);  
           return -2;
        }  
     }else /*child*/
     { setpgid(0, 0);
       err=system(command);
       if(WIFEXITED(err))
       { err=WEXITSTATUS(err);      
         if(err==255) err=254; 
         exit(err);
       } exit(255);  
     }      
   }
   if(err<0||WIFSIGNALED(err) ) err==-1; else
   {
      err=WEXITSTATUS(err);
      if(err==255) err=-1;
   }
   if(err<0)     
   {  FError=1;
      printf("Can not execute \"%s\"\n", command);
   }  
   return err;  
} 

/* File name  for SLHA input */
static char * FName=NULL;

/* Creates  SLHA input file. If file exist then it is cleaned  */ 

int openAppend(char * fileName)
{ FILE*f;
  if(access(fileName,F_OK)==0) unlink(fileName);
  FName=realloc(FName, strlen(fileName)+1);
  sprintf(FName,"%s",fileName);
  f=fopen(fileName,"a");
  if(!f) { FError=1; return -1;}
  fclose(f);  
  return 0;
} 

/* Open SLHA input file, Add record to the end of file. Close the file. */

int aPrintF(char * format, ...)
{  int err;
   FILE* f; 
   va_list args;
   char text[2*STRSIZ];
   if(!FName) {FError=1; return 0;}  
   f=fopen(FName,"a");
   if(f==NULL) {FError=1;return -1;}   
   va_start(args,format);
   vsprintf(text,format,args);
   va_end(args);
   err=fprintf(f,"%s",text);
   fclose(f);
   return err;
}
