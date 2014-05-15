#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>

#include "chep_crt.h"
#include "syst2.h"
#include "physics.h"
#include "screen.h"
#include "lbl.h"
#include "constr.h"
#include "read_mdl.h"
#include "sos.h"
#include "process.h"
#include "s_files.h"
#include "squar.h"
#include "r_code.h"
#include "c_out.h"
#include "red_out.h"
#include "math_out.h"
#include "form_out.h"
#include "symbolic.h"
#include "m_utils.h"
#include "amplitudes.h"
#include "viewdir.h"

#include "saveres.h"
#include "out_serv.h"
#include "fcompare.h"
#include "numcheck.h"
#include "../../include/version.h"
#include "dynamic_cs.h"

static int errorcode=0;

static void xw_error(void) { sortie(80);}
int    menulevel;

static  void checkAuxDir( int nmod)
{  char  n1[50],n2[50];
   char *fList[5]= {"vars","prtcls","extlib","func","lgrng",};   
   int i;
   for(i=0;i<5;i++)
   { sprintf(n1,"models/%s%d.mdl",fList[i],nmod); 
     sprintf(n2,"results/aux/models/%s1.mdl",fList[i]);
     if(fcompare(n1,n2)) break;
   }   
   if(i!=5)
   {
     delAllLib();
     system("if(test -d results/aux) then\n"
            "rm  -r results/aux\n"
            "fi\n"
            "mkdir results/aux; cd results/aux;\n"
            " mkdir models tmp results so_generated");
     for(i=0;i<5;i++)
     { char command[100];
       sprintf(command,"cp models/%s%d.mdl results/aux/models/%s1.mdl", fList[i],nmod,fList[i]);
       system(command); 
     } 
   }
}  


int main(int argc,char** argv)
{   

/*===================================
0 - First start.
1 - Model menu.
2 - Enter process menu.
3 - Feynman diagrams menu; squaring.
4 - Squared diagram menu; symbolic calculation.                        
5 - Write results; new process.
10 -Restart symbolic calculations
==========================================================*/

/* 0-Start; 1-Restart; 2-Heap Error,3-Edit Model,4-UserBreak */

  void *pscr1=NULL,*pscr2=NULL,*pscr3=NULL,*pscr4=NULL,*pscr5=NULL;
  int   k1=1,k2=1,k3=1,k4=1,k5=1;

  int n;
  int pid=0;
  char LOCKtxt[]="Directory 'results/' contains the .lock file created by n_calchep.\n"
                 "To continue you may a)Close the n_calchep session, or b)Rename 'results/',\n";


  blind=0;
  strcpy(pathtocalchep,argv[0]);
  for(n=strlen(pathtocalchep)-1; n && pathtocalchep[n]!=f_slash; n--);
  pathtocalchep[n-3]=0;                                         
     
   for ( n=1;n<argc;n++) 
   { if (strcmp(argv[n],"-blind")==0 && n<argc-1 )
     {  blind=1;
        inkeyString=argv[++n];
     }
     if (strcmp(argv[n],"+blind")==0 )  blind=2;                                     
     if (strcmp(argv[n],"--version")==0 ) { printf("%s\n", VERSION_); exit(0);}
   }    

   if(!writeLockFile(".lock")) 
   { fprintf(stderr,"locked by other s_calchep. See .lock file\n");
      exit(100);
   }
   strcpy(pathtouser,"");  
   sprintf(pathtohelp,"%shelp%c",pathtocalchep,f_slash);
   outputDir="results/";        
   { char * icon=(char *) malloc(strlen(pathtocalchep)+20);
     char  title[30];
     int err;
     sprintf(icon,"%s/include/icon",pathtocalchep);
     sprintf(title,"CalcHEP_%s/symb", VERSION);
     err=start1(title,icon,"calchep.ini",&xw_error);
     if(err && blind==0)
     { printf("Error:You have launched interactive session for version compiled  without X11 library.\n");
       printf(" Presumably X11 development package is not installed in your computer.\n");
       printf(" In this case directory /usr/include/X11/ is empty.\n");       
       printf("Options: a) Use blind session; b) Update Linux and recompile CalcHEP \n");
       printf("Name of needed package\n");
       printf("     libX11-devel      Fedora/Scientific Linux/CYGWIN/Darwin(MAC)\n");
       printf("     libX11-dev        Ubuntu/Debian\n");
       printf("     xorg-x11-devel    SUSE\n");

       exit(66);
     }      
     free(icon);
   }
   fillModelMenu();
   
   f3_key[0]=f3_key_prog;   f3_mess[0]="Model"; 
   f3_key[1]=f4_key_prog;   f3_mess[1]="Diagrams"; 
   f3_key[2]=f5_key_prog;   f3_mess[2]="Switches";
   f3_key[3]=f6_key_prog;   f3_mess[3]="Results"; 
                            f3_mess[4]="Del"; 
                            f3_mess[5]="UnDel";
   f3_key[6]=f9_key_prog;   f3_mess[6]="Ref";    
   f3_key[7]=f10_key_prog;  f3_mess[7]="Quit";

   restoreent(&menulevel);

   if(!blind && menulevel<2) cheplabel(); 

   switch (menulevel)
   {
      case 10: 
      case 6:
      case 5: 
      case 4:
      case 3: readModelFiles("./models",n_model); 
              modelinfo();
              loadModel(0,forceUG);
              processinfo();
              diagramsinfo();
      case 2: k1=n_model;

      case 1: break;
   }

   switch (menulevel)
   {
      case 2:  menuhelp();
               goto label_20;
      case 3:  goto label_31;
      case 4:  goto label_40;
      case 5:  goto label_50;
      case 6:  
      case 10: goto restart2;
   }
        
label_10:   /*   Menu2(ModelMenu): */
   f3_key[0]=NULL; /*models*/ 
   f3_key[1]=NULL; /*diagrams*/
   menulevel = 1;
   forceUG=0;
   menuhelp();
   for(;;)
   { 
      showheap();
      k1=n_model;
      menu1(56,4,"",modelmenu,"s_1",&pscr1,&k1);
      n_model=k1;
      if(n_model == 0)
      {
	if( mess_y_n(56,4,"Quit session")) {  saveent(menulevel); goto exi; }         
      }
      else  if(n_model == maxmodel+1)
      {
         clrbox(1,4,55,18);
         makenewmodel();
         menuhelp();
      }
      else if (n_model > 0)
      { 
	put_text(&pscr1);
	goto label_20;
      }
   }

label_20:   /*  Menu3:Enter Process  */
   f3_key[0]=NULL; 
   f3_key[1]=NULL; 

   menulevel = 2;
   saveent(menulevel);
   
   if(readModelFiles("./models",n_model)) 
   { 
     if(blind) sortie(133);
     goto label_10;
   }   
   modelinfo();
   k2 = 1;
   
   do
   {  char strmen[]="\026"
        " Enter Process        "
        " Force Unit.Gauge= OFF"
        " Edit model           "
        " Numerical Evaluation "
        "======================"
        " Delete model         ";

      if(forceUG)improveStr(strmen,"OFF","ON");
      menu1(56,4,"",strmen,"s_2_*",&pscr2,&k2);
      switch (k2)
      {
         case 0:  goto_xy(1,1); clr_eol(); goto label_10;
         case 2:  forceUG=!forceUG;   modelinfo(); break;
	 case 3:  editModel(1); break;
	 case 4:  numcheck();   
	 case 5:  break;
         case 6: 
	    if ( deletemodel(n_model) )
            {
               goto_xy(1,1);
               clr_eol();
               n_model=1;
               fillModelMenu();
               goto label_10;
            }
            else   readModelFiles("./models",n_model);
      }
   }  while (k2 != 1);

   loadModel(0,forceUG);
   
label_21:

   f3_key[0]=NULL; 
   f3_key[1]=NULL; 
   
   menulevel=2;
   checkAuxDir(n_model);
   errorcode=enter();   /*  Enter a process  */
   newCodes=0;
   showheap();
   if (errorcode)   /*  'Esc' pressed  */ { menuhelp(); goto label_20;}
   errorcode=construct();          /*  unSquared diagrams  */
   if (errorcode) 
   {  if(blind)
      {  printf("Processes of this type are absent\n"); sortie(111); } 
      else 
      { messanykey(5,22,"Processes of this type are absent");  
        clrbox(1,19,80,24); 
        goto label_21;
      }
   }
   else if(!blind)
   { int dirStat=checkDir("results"); 
     if(dirStat!=0)
     {  messanykey( 10,10,"There are files in directory 'results/'.\n"  
                          "To continue you has to clean or rename this directory.");
        viewresults();
        if(checkDir("results")!=0)  goto label_21;                   
     }
     clr_scr(FGmain,BGmain);
     modelinfo();
     processinfo();
     diagramsinfo();
     goto label_31;	
   }
   
label_30: /*  Menu4: Squaring,...*/
   clr_scr(FGmain,BGmain);
   modelinfo();
   processinfo();
   diagramsinfo();
label_31: 

   f3_key[0]=f3_key_prog; 
   f3_key[1]=NULL; 

   menulevel=3;  
   k3 = 1;
   do
   {
      menu1(56,4,"","\026"
         " View diagrams        "
/*         " Amplitude calculation" */
         " Square diagrams      "
         " Write down processes ","s_squa_*",&pscr3,&k3);
      switch (k3)
      {
         case 0: clrbox(1,2,55,11); menuhelp(); goto label_20;
         case 1: viewfeyndiag(1);   break;
         case 3: { FILE*f=fopen("results/list_prc.txt","w");
                   int k,ndel,ncalc,nrest;
                   char process[100];
                   long recpos; 
                   menup=fopen(MENUP_NAME,"r");
                   for(k=1;;k++) 
                   { int err=rd_menu(1,k,process,&ndel,&ncalc,&nrest,&recpos);
                     if(!err) break;
                     trim(process);
                     fprintf(f,"%s\n",process);
                   } 
                   fclose(f);
                   fclose(menup);
                    messanykey(20,14,"See file 'results/list_prc.txt'");
 
                 } break; 

/*       case 2: messanykey(10,10,"Not implemented yet"); Amplitudes(); */
      }
   }  while (k3 != 2);      

   if (!squaring()) goto label_30;  /*  process is absent  */

   clear_tmp();

   saveent(menulevel);
   restoreent(&menulevel);  

label_40:   /*  Menu5: View squared diagrams.....   */

   f3_key[0]=f3_key_prog; 
   f3_key[1]=f4_key_prog; 

   menulevel=4;
   clr_scr(FGmain,BGmain);
   modelinfo();
   processinfo();
   diagramsinfo();
   sq_diagramsinfo(); /*   ????????   */

   k4=1;
   saveent(menulevel);
   pscr4=NULL;
   for(;;)   
   {  int res;
      menu1(56,4,"","\026"
         " View squared diagrams"
         " Symbolic calculations"
         " Make&Launch n_calchep"        
         " Make n_calchep       "
         " REDUCE program       "
         ,"s_calc_*",&pscr4,&k4);

      res=checkDir("results");
      if(res==1)
      {
        int  n_calchep_id=setLockFile("results/.lock");
        if(n_calchep_id) unLockFile(n_calchep_id); else res=2;                          
      }
      switch (k4)
      {  case 0:
            if (mess_y_n(50,3,"Return to previous menu?"))goto label_30;
         break;

         case 1:
            viewsqdiagr();  break;

         case 2:     /*  Compute all diagrams   */
restart2:
            f3_key[0]=f3_key_prog; 
            f3_key[1]=f4_key_prog; 

            menulevel=4;
            calcallproc();
            sq_diagramsinfo();

            if(!continuetest()) break; 
            showheap();  
            put_text(&pscr4);
         break;
         case 3:
            { static char keystr[12]="]{{[{";
              if(res==2) messanykey(3,10,LOCKtxt); else inkeyString=keystr;
            }
            break;
         case 4:
            if(res==2) messanykey(3,10,LOCKtxt); else
            { char command[100];
              saveent(menulevel);
              chdir("results");
              makeVandP(0,"../models",n_model,2,pathtocalchep);
              finish();
              sortie(22);
            }
         case 5:  mk_reduceprograms();  break;
      }
      if(k4==2 && continuetest()) { put_text(&pscr4); break; }
   }

label_50:
   k5=1;
   pscr5=NULL;
   menulevel=5;
   saveent(menulevel);
   
   for(;;)  
   {  int n_calchep_id;   
      menu1(56,4,"","\026"
         " C code               "
         "     C-compiler       "
         "     Edit Linker      "
         " REDUCE code          "
         " MATHEMATICA code     "
         " FORM code            "
         " Enter new process    " ,"s_out_*",&pscr5,&k5);
         
      if((k5==1||k5==2)&&pid) 
      { int epid=waitpid(pid,NULL,WNOHANG);
        if(epid) pid=0; else
        { 
          messanykey(3,10,"This option is frozen while n_calchep runs or is compiled");
          continue;
        }
      }  

      switch (k5)
      {  case 0: goto label_40;  break;
         case 1:
           c_prog(); newCodes=0; saveent(menulevel);
           break;
         case 2:
           if(newCodes) { c_prog(); newCodes=0; saveent(menulevel);}
           fflush(NULL); 
           pid=fork();
           if(pid==0)
           { if(chdir("results")==0)
             { 
                n_calchep_id=setLockFile(".lock"); 
                if(n_calchep_id)
                {  char*command=malloc(100+strlen(pathtocalchep));
                   makeVandP(0,"../models",n_model,2,pathtocalchep); 
                   sprintf(command,"xterm -e %s/sbin/make__n_calchep %d", pathtocalchep,n_model);                                 
                   system(command);
                   sprintf(command," ./n_calchep"); 
                   system(command); 
                   chdir("../");
                   unLockFile(n_calchep_id);
                   free(command);
                }   
             }
             exit(0);
           }
           break;
         case 3: if(edittable(1,1,&modelTab[4],1," ",0))
                 {  char fName[STRSIZ];
                    sprintf(fName,"%smodels%c%s%d.mdl",pathtouser,f_slash,mdFls[4],n_model);
                    writetable( &modelTab[4],fName); 
                 }
                 break;
         case 4: makeReduceOutput(); break;
         case 5: makeMathOutput();   break; 
	 case 6: makeFormOutput();   break;
   
         case 7:
            put_text(&pscr5);
            clrbox(1,2,55,11);
            menuhelp();
            goto label_20;
      }
   }

exi:
   finish();
   sortie(0);
   return 0;
}
