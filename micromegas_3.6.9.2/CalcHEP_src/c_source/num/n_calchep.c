/*
 Copyright (C) 1997, Alexander Pukhov, e-mail pukhov@theory.npi.msu.su 
*/

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <signal.h>

#include "chep_crt.h"
#include "files.h" 
#include "rw_sess.h"
#include "interface.h"
/*#include "num_in.h" */
#include "n_calchep_.h"
#include "viewdir.h"
#include "rootDir.h"    
#include "../../include/num_out.h"
#include "../../include/version.h"
#include "read_func.h"
#include "rd_num.h"
#include "parser.h"
#include "dynamic_cs.h"
#include "param.h"

#ifdef OLD
double aWidth(char * pName)
{  int dim;
   txtList LL;
   if(!compDir)
   { int size=100;
     for(;;)
     {  compDir=realloc(compDir,size+20); 
        if(getcwd(compDir,size)) break; else size*=2;
     }
     strcat(compDir,"/aux");
     libDir=malloc(strlen(compDir)+20);
     sprintf(libDir,"%s/so_generated",compDir); 
     modelNum=1;
     calchepDir=getenv("CALCHEP");
     if(!calchepDir) calchepDir=interface_ext.CALCHEP;
     ForceUG=interface_ext.forceUG;
   }


   { double r=  pWidth(pName,&LL,&dim);
/*     printf("r=%E\n",r); */
     return r;
   }  
}
#endif

static void f5_key_prog(int x)
{
  int kmenu=1;
  void * pscr=NULL;

  while(kmenu) 
  {  
    char strmen[]="\040"
                  " Virtual   W decays          OF5"
                  " Virtual   Z decays          OF6";
        
    if(VWdecay)  improveStr(strmen,"OF5","ON ");
        else     improveStr(strmen,"OF5","OFF"); 
    if(VZdecay)  improveStr(strmen,"OF6","ON ");
        else     improveStr(strmen,"OF6","OFF");         
                                          
    menu1(20,18,"Switches",strmen,"s_switch_*",&pscr,&kmenu);
    if(kmenu==1) 
    {  VWdecay=!VWdecay; if(checkParam()>0){VWdecay=!VWdecay;checkParam();};}  
    if(kmenu==2) 
    {  VZdecay=!VZdecay; if(checkParam()>0){VZdecay=!VZdecay;checkParam();};}          
  }                     
}



                                                                    
static void f6_key_prog (int x){  viewDir("."); }
                                                                        

static void f10_key_prog (int x)
{
    if( mess_y_n(15,15," Quit session? ")) 
    {
        w_sess__(NULL);
        finish();
        sortie(0);
    }
}

static void f8_key_prog(int x)
{  
  static char FUNC[75]="2*2  % Press ESC to finish, F1 for help, ENTER to calculate";  
  int npos=1;
  void * pscr;
  get_text(1,20,80,21,&pscr);
  scrcolor(Red,White);   
  goto_xy(1,20); print("CALC : ");
  goto_xy(1,21); print("result=");
  scrcolor(Black,White);
  for(;;)
  { double res;
    int err; 
    int key;
    goto_xy(9,20); key=str_redact(FUNC,npos,70); 
    if(key==KB_ESC) break;
    else if(key==KB_F1)  show_help("n_calc");
    else if(key==KB_ENTER)
    {  goto_xy(8,21); 
       err=calcExpression(FUNC,rd_num,&res);
       goto_xy(9,21);
       if(err) {print("Erorr: %s",errmesstxt(err)); npos=rderrpos;}
       else  print("%E                     ",res); 
    }
  }
  put_text(&pscr);
}

static void f9_key_prog(int x)
{
  FILE*f;
  char fname[200];
  sprintf(fname,"%s%cCITE",pathtocalchep,f_slash);
  f=fopen(fname,"r");
  showtext (1, 1, 80,1,"",f);
  fclose(f);
}

static void xw_error(void) {sortie(80);}
  
int main(int argc,char** argv)
{
  int n;
  char icon[200];
  char title[30];

  blind=0;
  for( n=1;n<argc;n++) 
  {
    if(strcmp(argv[n],"--version")==0)  {  printf("%s\n",VERSION_); exit(0); }
    if(strcmp(argv[n],"-blind")==0&& n<argc-1 )
    { blind=1;     
      inkeyString=argv[++n];
    } else  if (strcmp(argv[n],"+blind")==0 ) blind=2;
  }

  if(!writeLockFile(".lock"))
  { fprintf(stderr,"locked by other n_calchep. See .lock\n");
    exit(100);
  }
                 
  setenv("CALCHEP",rootDir,0);
  sprintf(pathtocalchep,"%s%c",rootDir,f_slash);
  sprintf(pathtohelp,"%shelp%c",pathtocalchep,f_slash);
  sprintf(icon,"%s/include/icon",pathtocalchep);
  sprintf(title,"CalcHEP_%s/num", VERSION);

  f3_key[2]=f5_key_prog;   f3_mess[2]="Options";
  f3_key[3]=f6_key_prog;   f3_mess[3]="Results";
  f3_key[5]=f8_key_prog;   f3_mess[5]="Calc";
  f3_key[6]=f9_key_prog;   f3_mess[6]="Ref";
  f3_key[7]=f10_key_prog;  f3_mess[7]="Quit";
  
  { int size=100;
     for(;;)
     {  compDir=realloc(compDir,size+20);
        if(getcwd(compDir,size)) break; else size*=2;
     }
     strcat(compDir,"/aux");
     libDir=malloc(strlen(compDir)+20);
     sprintf(libDir,"%s/so_generated",compDir);
     modelNum=1;
     calchepDir=getenv("CALCHEP");
     if(!calchepDir) calchepDir=interface_ext.CALCHEP;
     ForceUG=interface_ext.forceUG;
  }
  
  
  

/* **  initialization of the session */
  link_process(PtrInterface_ext);  
  
  start1(title,icon,"calchep.ini;../calchep.ini",&xw_error); 
  r_sess__(NULL); 

  goto_xy(10,10); print("Calculation of constraints.  Please, be patient.");
  escpressed();
  n_comphep();
  finish();  
  sortie(0);
}
