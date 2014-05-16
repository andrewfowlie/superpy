#include <dlfcn.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "phys_val.h"
#include "VandP.h"
#include "parser.h"
#include "num_out.h"
#include "rootDir.h"
#include "qcdScale.h"

static int NX;
static FILE*fout;
static double (*scaleCC_)(REAL*,double (*calcPV)(char,char*,double*), double *,double*,double*) =NULL;
static void* scaleLib=NULL;


static void*  rd_num_(char* s)
{
   char*p;
   char key, plist[20];
   int i;
   
   if(isdigit(*s))
   { p=malloc(strlen(s)+10); 
     sprintf(p,"%s",s);
     return p;
   }  
   else 
   { p=malloc(20);
     sprintf(p,"X[%d]",NX++);

     if(checkPhysVal(s,&key,plist)) 
     { fprintf(fout,"  %s=calcPV('%c',\"",p,key);
       for(i=0;plist[i];i++) fprintf(fout,"\\%d",plist[i]);
       fprintf(fout,"\",pvect);\n");   
        return p;
     }

     for(i=0;i<nModelVars+nModelFunc;i++) if(strcmp(varNames[i],s)==0)
     { fprintf(fout,"%s = modelVal[%d];  /* %s  */ \n",p,i, varNames[i]); 
       return p;
     }    
   }  
   return NULL;
}

static void*  rd_num_R(char* s)
{ char*p;
  if(strcmp(s,"Qfact")==0){ p=malloc(20); strcpy(p,"Qfactorization"); return p;}
  return rd_num_(s); 
}

static void*  act_num_(char* ch,int n, void**args)
{  char*  p=NULL,*p1=NULL, *p2=NULL,*p3=NULL;
   p1= (char*)args[0];
   if(n>=2) p2=(char*)args[1];
   if(n>=3) p3=(char*)args[2];
   
   switch (ch[0])
   { 
     case '-': 
     case '+': 
     case '*': p=malloc(strlen(p1)+strlen(p2)+10);
               sprintf(p,"(%s)%c(%s)",p1,ch[0],p2);
               return p;
     case '/':p=malloc(strlen(p1)+strlen(p2)+15);
              sprintf(p,"(%s)/(double)(%s)",p1,p2);  
              return p;           
     
     case '^': p=malloc(strlen(p1)+strlen(p2)+10);
               sprintf(p,"pow(%s,%s)",p1,p2);
               return p;
     case '.': rderrcode=typemismatch; 
               return NULL; 
     default : switch(n) 
               {      
                   case 1: 
                      if(!strcmp(ch,"sqrt")||!strcmp(ch,"sin")||!strcmp(ch,"cos")||!strcmp(ch,"tan")
                       ||!strcmp(ch,"asin")||!strcmp(ch,"acos")||!strcmp(ch,"atan")||!strcmp(ch,"exp")
                       ||!strcmp(ch,"log")||!strcmp(ch,"fabs"))
                       { 
                         p=malloc(strlen(p1)+10);
                         sprintf(p,"%s(%s)",ch,p1);
                       }   
                       return p;
                   case 2:
                      if(!strcmp(ch,"atan2"))
                      { p=malloc(strlen(p1)+strlen(p2)+10);
                        sprintf(p,"atan2(%s,%s)",p1,p2);
                        return p;   
                      }
                   case 3:   
                      if(!strcmp(ch,"if"))
                      { p=malloc(strlen(p1)+strlen(p2)+strlen(p3)+10);
                        sprintf(p,"(%s>0 ? %s:%s)",p1,p2,p3);
                        return p;
                      }  
               }
    }                    

   if(strcmp(ch,"min")==0 || strcmp(ch,"max")==0)
   { int l=0,i;
     char *q;
     for(i=0;i<n;i++) l+=strlen((char*)args[i])+8;
     p=malloc(l);
     q=malloc(l);
     strcpy(p,(char*)args[0]);
     for(i=1;i<n;i++)
     { strcpy(q,p);
       sprintf(p,"%s(%s,%s)",ch,q,(char*)args[i]); 
     }
     free(q);
     return p;
   }
   return NULL;
}

int initScales(char*sF,char*sR,char * mess)
{  char *ch;
   char *command;
    
   if(scaleLib) dlclose(scaleLib); 
   scaleCC_=NULL;
   
   fout=fopen("scale.c","w");
   if(!fout){ if(mess) sprintf(mess,"can't open file scale.c for writing");  return -1;}
   int pos;
   NX=0;
   fprintf(fout,"#include<math.h>\n");  
   fprintf(fout,"#include\"%s/include/nType.h\"\n",rootDir);
   fprintf(fout,"#define min(x,y) (x<y? x:y)\n");
   fprintf(fout,"#define max(x,y) (x>y? x:y)\n"); 
   fprintf(fout,"extern void ScaleCC(REAL*, double (*calcPV)(char,char*,double*), double*,double*,double*);\n");
   fprintf(fout,"void ScaleCC(REAL*modelVal, double (*calcPV)(char,char*,double*), double *pvect,double *QF,double *QR)\n");
   fprintf(fout,"{ double ");
   pos=ftell(fout);
   fprintf(fout,"         \n");
   fprintf(fout," double Qfactorization;\n");   
   ch=readExpression(sF,rd_num_,act_num_,free);
   
   if(!ch) 
   { fclose(fout); 
     if(mess)
     { 
       sprintf(mess,"Error in Qfact scale definition: position %d :\n%s", 
                    rderrpos, errmesstxt(rderrcode));
     }
     return rderrpos;    
   } 
   fprintf(fout," *QF= %s;\n",ch);
   fprintf(fout," Qfactorization=%s;\n",ch);  
   ch=readExpression(sR,rd_num_R,act_num_,free);
   
   if(!ch) 
   { fclose(fout); 
     fout=NULL;
     if(mess)
     { 
       sprintf(mess,"Error in Qren scale definition: position %d :\n%s", 
                    rderrpos, errmesstxt(rderrcode));
     }
     return rderrpos;    
   } 
   fprintf(fout," *QR= %s;\n",ch);
      
   fprintf(fout,"}\n");
   fseek(fout,pos,SEEK_SET);
   fprintf(fout,"X[%d];",NX+1);
   fclose(fout);
   free(ch);
   
    command=malloc(strlen(rootDir)+100);
    sprintf(command,". %s/FlagsForSh; $CC $CFLAGS $SHARED -o scale.so scale.c",rootDir);
    system(command);
    free(command);
    scaleLib=dlopen("./scale.so",RTLD_NOW);
    if(!scaleLib)
    {
       if(mess) sprintf(mess,"Can't load shared library for QCD scale"); 
       return -3;
    } 
    scaleCC_=dlsym(scaleLib,"ScaleCC");
    if(!scaleCC_)
    {  
       if(mess) sprintf(mess,"Problem in library scale.so ");
       return -4;
    }
    return 0;
}


void Scale(double*pv,double *qF, double *qR)
{ double q=91.187; 

  if(scaleCC_)  {(*scaleCC_)(varValues, calcPhysVal, pv,qF,qR ); if(*qF<0) *qF*=-1; if(*qF<1) *qF=1; if(*qR<0) *qR*=-1; if(*qR<1) *qR=1;} 
  else {*qR=q; *qR=q;}
}
