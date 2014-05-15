
//#include"num0.inc"

#include "histogram.h"
#include "chep_crt.h"
#include "viewdir.h"
#include "files.h" 
#include "rootDir.h"
#include "../../include/version.h"
#include "interface.h"

static char **proclist;
static int nProc; 
static int topProc;
static int height=20;

static void f6_key_prog (int x){  viewDir("."); }


static void drawList(void)
{
  int i; 
  scrcolor(Black,White);
  for(i=0; i<height;i++)
  { goto_xy(3,i+4); 
    if(i+topProc<nProc)print(" %-40.40s",proclist[i+topProc]); 
          else print(" %-40.40s"," ");     
  }
  goto_xy(45,4);
  if(topProc>0)
  { scrcolor(Red,White); print("F7-up");}   else {scrcolor(Black,LightGray); print("     ");}
  goto_xy(45,2+height+1); 
  if(topProc+height<nProc) 
  { scrcolor(Red,White);   print("F8-down"); } else{ scrcolor(Black,LightGray);  print("       ");}

}

static void f7_key_prog (int x)
{
    if(topProc>0)
    { topProc-=height;
      drawList();
    }   

}

static void f8_key_prog (int x)
{
    if(topProc+height<nProc)
    { topProc+=height;
      drawList();
    }   
}



static void f10_key_prog (int x)
{ if( mess_y_n(15,15," Quit session? ")) { finish(); exit(0);} }    


    
 
int main(int argc, char ** argv)
{ 
   FILE*f;
   int nf=1;
   char *Process;
   int i;
   
   f3_key[3]=f6_key_prog;   f3_mess[3]="Results";
   f3_key[4]=f7_key_prog;   f3_mess[4]="Up"; 
   f3_key[5]=f8_key_prog;   f3_mess[5]="Down";
   f3_key[7]=f10_key_prog;  f3_mess[7]="Quit";
   blind=0;

   if(argc<2) 
   { printf(" This  routine is intended to display distributions produced\n" 
            " in  calcHEP numerical sessions (files dist_#).\n"
            " The name of file should submited to this routine as parameter\n");
     return 1;
   }

   if(strcmp(argv[1],"-blind")==0)
   { blind=1;     
     nf=3; 
     inkeyString=argv[2];
   } else  if (strcmp(argv[1],"+blind")==0 )
   { blind=2;
     nf=2;
   }

  if(nf!=argc-1) 
  { printf(" This  routine is intended to display distributions produced\n" 
           " in  calcHEP numerical sessions (files dist_#).\n"
           " The name of file should submited to this routine as parameter\n");
    return 1;
  }

  sprintf(pathtohelp,"%s%chelp%c",rootDir,f_slash,f_slash); 
  f=fopen(argv[nf],"r");
   
  if(!f)
  {  printf("Can not open the '%s' file\n",argv[nf]);
     exit(1);
  }   

  if(rdr_hist2(f,&Process))
  {  printf("Wrong  format for %s \n",argv[nf]);
         exit(1);
  }else     
  {  char*c1,*c2;
     c1=strstr(Process,",");
     c2=strstr(Process,"->");
     if(c1>c2) nin_int=1; else nin_int=2;
  }
  fclose(f);

  { char txt[80];
    sprintf(txt," %s Viewer for  distributions ",VERSION_);   
    start1(txt,pathtocalchep,"calchep.ini;../calchep.ini",NULL);
  }    
    
  goto_xy(35,2);  
  scrcolor(Blue,LightGray);     print(" File:"); 
  scrcolor(Red,LightGray);      print(" %s     ",argv[nf]); 
  goto_xy(10,2);    
  scrcolor(Blue,BGmain);        print(" Processes:");

  { char *ch1,*ch2;
    for(ch1=Process,nProc=0;; ch1=strchr(ch1,';'), nProc++ )if(ch1)ch1++;else break;
    proclist=malloc(nProc*sizeof(char*));
    for(i=0,ch1=Process;i<nProc;i++,ch1=strchr(ch1,';'))
    {  if(i) ch1++;
       ch2=strchr(ch1,';');
       if(ch2) proclist[i]=malloc(ch2-ch1+1); else proclist[i]=malloc(strlen(ch1)+1);
       sscanf(ch1,"%[^;]",proclist[i]);
    }
  }
       
  topProc=0; 
  drawList();
  
  showHist(54,5,Process);
  
  finish();

  return 0;
}
