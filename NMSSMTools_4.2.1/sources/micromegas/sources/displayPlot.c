
#include"micromegas.h"
#include"micromegas_aux.h"
#include"../CalcHEP_src/c_source/chep_crt/include/crt.h"
//#define OLD_VERSION
#ifdef OLD_VERSION

void displayPlot(char * title, char*xName, char*yName, double xMin, double xMax, int dim, double *f,double *ff)
{ 
  int i;
  char *buff;
  if(ff) i=2;else i=1;
 
  buff=malloc(100+ strlen(title)+strlen(xName)+strlen(yName)+strlen(calchepDir) +dim*(11*i+2));  
 
  sprintf(buff,"echo \"#type %d\n#title %s\n#xName %s\n#yName %s\n#xMin %.3E\n#xMax %.3E\n#xDim %d\n",
  i-1,title,xName,yName, xMin,xMax,dim);
  if(ff) for(i=0;i<dim;i++) sprintf(buff+strlen(buff),"%.2E %.2E\n",f[i],ff[i]);
  else   for(i=0;i<dim;i++) sprintf(buff+strlen(buff),"%.3E\n",f[i]); 
  
  sprintf(buff + strlen(buff),"\" | %s/bin/plot_view &",calchepDir);

  system(buff);
  free(buff);
}

void  killPlots(void){;}

#else 

#include <sys/types.h>
#include <unistd.h>
#include <signal.h>              
#include <sys/wait.h>
 
static int First=1;
  
static void disconnect(int N) { setsid();}
  

extern void   plot_1(double xMin, double xMax, int dim,
                       double *f, double *ff,char* upstr, 
                       char* xstr, char* ystr);
extern int blind;

  
static int newPID=0;
static int pidList[100];

extern char pathtocalchep[], pathtohelp[];
  
void displayPlot(char * title, char*xName, char*yName, double xMin, double xMax, int dim, double *f,double *ff)
{ int pid;
  
  if(First) { First=0;   signal(SIGUSR1, disconnect);}  
  
  pid=fork();
  if(pid==0) 
  {  int err;
     blind=0; 
     err=start1("micrOMEGAs Plot",NULL ,"calchep.ini",NULL);
     if(err) 
     { printf("Can not display plot because micromegas is compiled without X11\n");
       exit(0);
     }
     sprintf(pathtocalchep,"%s/",calchepDir);
     sprintf(pathtohelp,"%s/help/",pathtocalchep);
     clearTypeAhead();  
     plot_1(xMin,xMax,dim,f,ff,title,xName,yName);
     finish();
     exit(0);
  } else pidList[newPID++]=pid;
}


void  killPlots(void)
{
  int  C,i;
  
  for(i=0;i<newPID;i++) if(waitpid(pidList[i],NULL,WNOHANG)) break;

  if(newPID && i==newPID)
  {
    printf("Kill all plots (Y/N)? "); 
    C=getchar();
    if(C=='y'|| C=='Y')  kill(0,SIGKILL); else kill(0,SIGUSR1);
  }
  newPID=0;
}

#endif

void  killplots_(void) { killPlots();}

void displayFunc(double (*F)(double), double x1  ,double x2, char * mess)
{
  int i;
  double f[100];
  
  for(i=0;i<100;i++) f[i]=F(x1+i*(x2-x1)/99.);
  displayPlot(mess,"x","F(x)",x1,x2,100,f,NULL);
}  

void displayFunc10(double (*F)(double), double x1  ,double x2, char * mess)
{
  int i;
  double f[100];
  
  for(i=0;i<100;i++) f[i]=F(pow(10,x1+i*(x2-x1)/99.));
  displayPlot(mess,"x","F(10^x)",x1,x2,100,f,NULL);
}  
