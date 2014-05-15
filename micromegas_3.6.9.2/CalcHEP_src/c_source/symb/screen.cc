/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include"chep_crt.h"
#include "syst2.h"
#include "s_files.h"
#include "sos.h"
#include "showgrph.h"
#include "showgrph2.h"
#include "red_out.h"
#include "r_code.h"
#include "file_scr.h"
#include "read_mdl.h"
#include "crt_util.h"
#include "physics.h"
#include "screen.h"
#include "prepdiag.h"
#include "c_out.h"
#include "process.h"
#include "cweight.h"
#include "procvar.h"
#include "viewdir.h"
#include "dynamic_cs.h"
#define tcol Green
#define mpos 7
#define graphpos 8

int    menulevel;

static void diag_stat(int type,int*n_sub,int*n_del,int*n_calc,int*n_rest)
{
  int ndel,ncalc,nrest;
  long nn;
  char buff[50];
  FILE * tmp, **menu;

  *n_sub=0, *n_del=0, *n_calc=0, *n_rest=0;

  if(type==1) { menu=&menup; strcpy(buff,MENUP_NAME);}
  else        { menu=&menuq; strcpy(buff,MENUQ_NAME);}

  tmp=*menu;
  *menu=fopen(buff,"rb"); 
  if(*menu) 
  {
     while (rd_menu(type, *n_sub+1 ,buff, &ndel,&ncalc,&nrest,&nn)) 
     {  
        (*n_sub)++;
        *n_del +=ndel;
        *n_calc+=ncalc;  
        *n_rest+=nrest;
     }
     fclose(*menu);
  }
  *menu=tmp; 
  return; 
}


void writeModelFiles(int l)
{
  int i;
  char fName[STRSIZ];
  for(i=0;i<5;i++)
  {
    sprintf(fName,"%smodels%c%s%d.mdl",pathtouser,f_slash,mdFls[i],l);
    writetable( &modelTab[i],fName);
  }
}



void  editModel( int  edit)
{ int  n,i,j;
  void *  pscr = NULL;
  int edited=0,renamed=0;
  int ok=1;
  char  tabMenu[STRSIZ], tabName[80];
  char menuName[30]; 
  char * tabhelp[5]={ "s_mdl_1", "s_mdl_2", "s_mdl_3", "s_mdl_4","s_mdl_5"};

  n = 1;
 
cont:
  do
  { 
    if(edit) strcpy(menuName,"Edit  model"); else strcpy(menuName,"View   model");
    strcpy(tabMenu,"\017");
    for (i=0;i<5;i++)
    { strcpy(tabName,modelTab[i].headln);
      trim(tabName);
      for(j=strlen(tabName);j<14;j++) tabName[j]=' ';
      tabName[14]=0;
      strcat(tabMenu," "); strcat(tabMenu,tabName);
    }
    if(edit)  strcat(tabMenu,"      RENAME   "
                             "   CHECK MODEL ");
    if(edit) menu1(56,7,"",tabMenu,"s_mdl_e",&pscr,&n);
       else  menu1(17,15,"",tabMenu,"s_mdl_v",&pscr,&n); 
    switch(n)
    { case 0:
        if( edited || renamed )
        {
          if(mess_y_n(15,19," Save corrections ?") )
          {
            if (ok||loadModel(1,0) ) writeModelFiles(n_model); else  goto cont;
            if(renamed) 
            { int i,size=modelmenu[0]; 
              for(i=0;modelTab[0].mdlName[i]&&i<size-1;i++)
              modelmenu[(n_model-1)*size +2+i]=modelTab[0].mdlName[i];
              for(;i<size-1;i++) modelmenu[(n_model-1)*size +2+i]=' ';
              modelinfo();
            }       
          }else  readModelFiles("./models",n_model);   
        }
        return;
    
      case 1: case 2: case 3: case 4: case 5:    
        edited=edittable(1,1,&modelTab[n-1],1,tabhelp[n-1],!edit) || edited;
	if(edited) ok=0;
    	break; 
      case 6:
        if(correctStr (10,17,"New name:",modelTab[0].mdlName,29,1))
        { trim(modelTab[0].mdlName);
          for(i=1;i<5;i++)strcpy(modelTab[i].mdlName,modelTab[0].mdlName);
          renamed=1;
        }         
        break;
      case 7:
        if (!loadModel(1,0) ) goto cont; 
        else {ok=1; messanykey(10,15,"The model is Ok");}
        break;
    }       
  }while(n);

}


void         showheap(void)
{ /* goto_xy(60,24);  print("Memory= %lu",usedmemory);*/ }


void         menuhelp(void)
{
  scrcolor(Red,BGmain);
  goto_xy(23,4);print("Abstract");
  scrcolor(FGmain,BGmain);
  goto_xy(1,6);
  print("     CalcHEP package is created for calculation of \n"); 
  print("  decay and high energy collision processes of     \n"); 
  print("  elementary particles in the lowest order (tree)  \n");
  print("  approximation. The main idea put into the CalcHEP\n");
  print("  was to make available passing from the lagrangian\n"); 
  print("  to the final distributions effectively with the  \n");
  print("  high level of automatization.\n");
  print("     Use F2 key to get information about interface \n");   
  print("  facilities and F1 - as online help.              \n");
 
  scrcolor(Black,BGmain);
  chepbox(1,5,53,15);
  scrcolor(FGmain,BGmain);
}

char * currentModelName(void)
{ static char name[40];
  int size=modelmenu[0];
  if(n_model<=0) name[0]=0; else 
    sprintf(name,"%*.*s",size,size,modelmenu+(n_model*size-size+1));
  trim(name);
  return name;
}

void modelinfo(void)
{ 
  goto_xy(5,1);
  scrcolor(Red,BGmain); print("   Model:  ");
  scrcolor(FGmain,BGmain);
  print(currentModelName());
  if(forceUG) print("/Unitary Gauge/");
  else        print("                  ");
}


void         processinfo(void)
{
   goto_xy(5,3);
   scrcolor(Red,BGmain); print(" Process:  ");
   scrcolor(FGmain,BGmain);
   print("%s",processch);
}


void         diagramsinfo(void)
{  int n_sub, n_del, n_calc, n_rest;

   diag_stat(1,&n_sub,&n_del,&n_calc,&n_rest);
   if(!n_sub) return;
   goto_xy(15,5);
   scrcolor(Red,BGmain);
   print(" Feynman diagrams \n");
   scrcolor(FGmain,BGmain);
   print("      diagrams in      subprocesses  are constructed.\n");
   print("      diagrams  are deleted.");
   scrcolor(Blue,BGmain);
   goto_xy(1,6);
   print("%d",n_del+n_calc+n_rest);
   goto_xy(20,6);
   print("%d",n_sub);
   goto_xy(1,7);
   print("%u",n_del);
   while (where_y() < 7) print(" ");
   
}


void  sq_diagramsinfo(void)
{
   int n_sub, n_del, n_calc, n_rest;


   diag_stat(2,&n_sub,&n_del,&n_calc,&n_rest);
   if(!n_sub) return;

   goto_xy(15,9);
   scrcolor(Red,BGmain);
   print(" Squared diagrams \n");
   scrcolor(FGmain,BGmain);
   print("      diagrams in      subprocesses  are constructed.\n");
   print("      diagrams  are deleted.\n");
   print("      diagrams  are calculated.");
   scrcolor(Blue,BGmain);
   goto_xy(1,10);
   print("%d",n_del+n_calc+n_rest);
   goto_xy(20,10);
   print("%d",n_sub);
   goto_xy(1,11);
   print("%d",n_del);
   while (where_x() < 7) print(" ");
   goto_xy(1,12);
   print("%d",n_calc);
   goto_xy(1,13);
   
}

static void f7_key_prog(int x)
{ 
  static char delstr[5]="{D}";
  inkeyString=delstr;   
}

static void f8_key_prog(int x)
{
  static char delstr[5]="}R{";
  inkeyString=delstr;   
}


static void menu_f (int col, int row,char* label, char* f_name, char * help,
                                             void* hscr,int * kk)
{
FILE * f;
int nline,i;

char * menustr;
char ch1,ch2;

f= fopen(f_name,"r");
if (f==NULL) return;

fread(&ch1,1,1,f);
fread(&ch2,1,1,f);

fseek(f,0,SEEK_END);
nline=ftell(f)/ch2;

menustr=(char *) malloc(2+nline*ch1);
fseek(f,2,SEEK_SET);

menustr[0]=ch1;
for (i=1;i<=nline;i++) 
{ fread(menustr + 1+(i-1)*ch1,ch1,1,f); 
  fseek(f,ch2-ch1,SEEK_CUR);
}
  fclose(f);
  menustr[1+nline*ch1]=0;  
  menu1 (col, row,label,menustr, help, hscr, kk);

}


void  sqdiagrmenu(void)
{  void * pscr = NULL ;

   if (subproc_sq == 1) { nsub = 1; return;}
   menu_f(5,15,"NN      Subprocess                                 Del   Calc  Rest ",
       MENUQ_NAME,"s_sq_proc",&pscr,&nsub);
  if(nsub) put_text(&pscr);

}


void  viewsqdiagr(void)
{
   nsub = 1;
   do
   {
      sqdiagrmenu();
      if (nsub != 0)
      { 
//      if(nin+nout<7)  
        showgraphs(2); 
//        else  messanykey(10,15,"The editor does not work if the number of legs exceed 6.");
      }
      sq_diagramsinfo();
   }  while (!(nsub == 0 || subproc_sq == 1));   /*  Esc  */
   sq_diagramsinfo();
}


void  viewfeyndiag(int del_mode)
{
   void * pscr = NULL;
   int upr= del_mode? 1:-1;
 
   if(del_mode){ f3_key[4]=f7_key_prog; f3_key[5]=f8_key_prog;}       
   nsub = 1;
   do
   {
      if (subproc_f == 1)  nsub = 1; else
      {
         menu_f(9,11,"NN        Subprocess                               Del   Rest ",
          MENUP_NAME,"s_proc",&pscr,&nsub);
      }
      if (nsub > 0) {/*showgraphs2(upr);*/ showgraphs(upr);  if(del_mode) diagramsinfo();}         
   }  while (!(nsub == 0 || subproc_f == 1));

   if(del_mode){ f3_key[4]=NULL; f3_key[5]=NULL;} 
}



int  viewresults(void)
{  
   int  k,kmenu;
   void *  pscr  = NULL;
 
   shortstr  newname;
   int dirStat=checkDir("results");

   if(dirStat==0){messanykey(10,15,"directory RESULTS is empty"); return 1;}

   kmenu = 1;   
label_1:

   menu1(10,10,"","\010"
      " View   "
      " Delete "
      " Rename ","s_res",&pscr,&kmenu);
      
   switch (kmenu)
   {
     case 0:
       return 0;         
     case 1: 
      viewDir("results");    
      break;

     case 2:
      if(dirStat==2)
      { char mess[]="Can not clean dir 'results' because it contains the LOCK file";
        if(blind) { printf("%s\n",mess); sortie(102);} 
        else { messanykey(3,13,mess); break;}
      } 
      if ( mess_y_n( 6,13," Delete files ") ) 
      {  struct dirent **namelist;
         int n,i;
         n = scandir("./results", &namelist, NULL, NULL);
         for(i=0; i<n;i++)
         { 
           char buff[100];
           if(strcmp(namelist[i]->d_name,"aux"))
           { 
             sprintf(buff,"results/%s",namelist[i]->d_name);  
             unlink(buff);
           }  
         }
         free(namelist);
      }
      put_text(&pscr);
      return 1;
     case 3:
      strcpy(newname," ");
      while(1)
      {  void * pscr3;
         get_text(1,maxRow(),maxCol(),maxRow(),&pscr3); 
         goto_xy(1,maxRow());
         print("Enter new name: ");
	 k = str_redact(newname,1,30);
	 if (k == KB_ESC)
	 {   goto_xy(1,24);
             clr_eol();
             goto label_1;
         }
	 if (k == KB_ENTER)
         {
            trim(newname);
            if(rename("results",newname)==0)
            {  char command[200];
               mkdir("results",-1);
               sprintf(command," cp -rp  %s/aux results",newname);
               system(command);
               
               put_text(&pscr);
               put_text(&pscr3);
               return 1;
            }
             else  messanykey(10,15," Can't rename the directory");
         }
         put_text(&pscr3);   
      }
   } 
   goto label_1;
}


void f3_key_prog(int x)
{
  int i;
  for(i=0; i<8 && f3_key_prog != f3_key[i] ; i++);
  if(i<8) f3_key[i]=NULL;  /* LOCK */

    editModel(0);

  if(i<8) f3_key[i]=f3_key_prog;  /* UNLOCK */

}



void f4_key_prog(int x)
{ int nsubtmp;
  int i;
  for(i=0; i<8 && f4_key_prog != f3_key[i] ; i++);
  if(i<8) f3_key[i]=NULL;  /* LOCK */
     nsubtmp = nsub;
     viewfeyndiag(0);
     nsub = nsubtmp;
  if(i<8) f3_key[i]=f4_key_prog;  /* UNLOCK */
}

void f5_key_prog(int x)
{
  int kmenu=1;
  void * pscr=NULL;

  int nfun;
  for(nfun=0; nfun<8 && f5_key_prog != f3_key[nfun] ; nfun++);
  if(nfun<8) f3_key[nfun]=NULL;
       
 
  
  while(kmenu) 
  {  
    char strmen[]="\040"
/*                  " Symbolic conservation low   OF1"   */
                  " Number of QCD colors =      Nc "
                  " Diagrams in C-output        OF3"
                  " Widths in t-channels        OF4"
                  " Virtual W/Z decays          OF5";
/*                   
    if(consLow) improveStr(strmen,"OF1","ON ");
       else     improveStr(strmen,"OF1","OFF");
*/
    if(NcInfLimit) improveStr(strmen,"Nc","Inf");
       else        improveStr(strmen,"Nc","3");
/*
    if(noCChain) improveStr(strmen,"OF2","OFF");
       else      improveStr(strmen,"OF2","ON ");
*/
    if(noPict) improveStr(strmen,"OF3","OFF");
        else      improveStr(strmen,"OF3","ON ");
    if(tWidths) improveStr(strmen,"OF4","ON ");
        else      improveStr(strmen,"OF4","OFF"); 
        
    if(VVdecay)  improveStr(strmen,"OF5","ON ");
        else     improveStr(strmen,"OF5","OFF");    
                                          
                  
    menu1(20,18,"Switches",strmen,"s_switch_*",&pscr,&kmenu);
    switch (kmenu)
    {
    
//      case 1: consLow=!consLow;       break;
      case 1: NcInfLimit=!NcInfLimit; break;
      case 2: noPict=!noPict;         break;
      case 3: tWidths=!tWidths;         break;     
/*      case 5: noCChain=!noCChain;     break; */
    }
    
  }
  if(nfun<8) f3_key[nfun]=f5_key_prog;
                     
}


void f6_key_prog(int x)
{
  int i;
  for(i=0; i<8 && f6_key_prog != f3_key[i] ; i++);
  if(i<8) f3_key[i]=NULL;  
  viewresults();
  if(i<8) f3_key[i]=f6_key_prog;
}

void f10_key_prog(int x)
{
  if ( mess_y_n(56,maxRow()-5," Quit session?  "))
  {
     saveent(menulevel);
     finish();
     sortie(0);
  }
}

void f9_key_prog(int x)
{
  FILE*f;
  char fname[200]; 
  sprintf(fname,"%s%cCITE",pathtocomphep,f_slash);
  f=fopen(fname,"r");
  showtext (1, 1, 80,1,"",f);
  fclose(f);
} 
