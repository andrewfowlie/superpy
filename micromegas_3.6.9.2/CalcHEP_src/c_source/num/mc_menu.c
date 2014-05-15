/*
 Copyright (C) 1997, Alexander Pukhov, e-mail pukhov@theory.npi.msu.su
*/

#include <unistd.h>
#include "chep_crt.h"
#include "runVegas.h"
#include "alphas2.h"
#include "kininpt.h"
#include "interface.h"
#include "cut.h"
#include "regul.h"
#include "edittab.h"
#include "rw_sess.h"
#include "strfun.h"
#include "width_12.h"
#include "cs_22.h"
#include "histogram.h"
#include "q_kin.h"
#include "subproc.h"
#include "n_calchep_.h"
#include "plot.h"
#include "param.h"
#include "vegas.h"
#include "mc_menu.h"
#include "4_vector.h"
#include "spectrum.h"
#include "num_in.h"
#include "comp.h"

static int sub_men__(void)
{
    int  n, npr, mode=0;
    char * strmen;
    void * pscr = NULL;
    int width,width_;
    int nprc_old=Nsub;
    
    if(nprc_int==1 && !blind ) return 0;
 
    width=0;
    for(npr = 1; npr <= nprc_int; ++npr)
    for(n=1;n<=nin_int+nout_int;n++) 
    width=MAX(width,strlen(pinf_int(npr,n,NULL,NULL)));
    
    width++;
    width_=6+width*(nin_int+nout_int);
       
    strmen=malloc(2+nprc_int*width_);

    for(n=1;n<=width_*nprc_int;n++) strmen[n]=' ';
    strmen[0]= width_;
    strmen[1+nprc_int*width_]=0;

    for(npr = 1; npr <= nprc_int; ++npr)
    { 
      for(n=1;n<=nin_int;n++) 
      { char *s=pinf_int(npr, n,NULL,NULL);
        memcpy(strmen+(npr-1)*width_+2+(n-1)*width ,s,strlen(s));
      }                

      memcpy(strmen+(npr-1)*width_+2+(nin_int)*width ," -> ",4);

      for(n=nin_int+1;n<=nin_int+nout_int;n++) 
      { char* s=pinf_int(npr, n,NULL,NULL);
          memcpy(strmen+(npr-1)*width_+6+(n-1)*width ,s,strlen(s));
      }                
    }         
    menu1(76-width_,6,"",strmen,NULL,&pscr,&mode);
    free(strmen);
           
    if(mode)put_text(&pscr);
    if(mode && mode!=nprc_old) 
    { 
       Nsub = mode;
       wrtprc_();
       if(nin_int==2) initStrFun(0);
       return 1;
    } 
    return 0;
}



static void f7_prog(int mode)
{ int pos=1;
  void *pscr=NULL;
  if(mode>2) 
  { messanykey(10,15," Highlight the corresponding\n"
                     "structure function");
    return;
  }     

  if(!sf_num[mode-1]) return;
  
  f3_key[4]=NULL;

  for(;;)
  {  static double xMin=0.0, xMax=1.0, scale = 91.187;
     static int nPoints=100;
     double f[150];
     
     char strmen[]="\030 "
                  " x-Min = XXX            "
                  " x-Max = YYY            "
                  " Npoints = NNN          "
                  " QCD-scale= QQQ         "
                  " Display plot x*F(x)    "
                  " Display plot F(x)      ";
     
     improveStr(strmen,"XXX","%.3f",xMin);
     improveStr(strmen,"YYY","%.3f",xMax);
     improveStr(strmen,"NNN","%d",nPoints);
     improveStr(strmen,"QQQ","%.1fGeV",scale); 

     menu1(54,14,"",strmen,"n_alpha_view",&pscr,&pos);

     switch(pos)
     {  case 0: f3_key[4]=f7_prog; 
                return;
        case 1:
          correctDouble(55,18,"xMin = ",&xMin,1);                  
          break; 
        case 2:
          correctDouble(55,18,"xMax = ",&xMax,1);                  
          break; 
        case 3:
          correctInt(50,18,"nPoints = ",&nPoints,1);
          break; 
        case 4: 
          correctDouble(50,18,"QCD-scale = ",&scale,1);
          break;
        case 5: case 6:
         if(xMin>=0 && xMax>xMin && xMax<=1 
            && nPoints>=3 && nPoints<=150 && scale>0.5)
         {
           void * screen;
           double dx=(xMax-xMin)/(2*nPoints);
           double be=sf_be[mode-1];
           int i;
           get_text(1,1,maxCol(),maxRow(),&screen);
           for(i=0;i<nPoints;i++)
           {  double x=xMin+(i+0.5)*(xMax-xMin)/nPoints;
              f[i]=strfun_(mode,x,scale);
              if(pos==5) f[i]*=x;
              if(be!=1.) f[i]*=be*pow(1.-x,be-1.);
           }
           {
              char p_name[20], mess[STRSIZ];
              strcpy(p_name,pinf_int(Nsub,mode,NULL,NULL));
              if(pos==5) strcat(p_name,"(x)*x"); else strcat(p_name,"(x)");
              strFunName(mode,mess);
              trim(mess);
              sprintf(mess+strlen(mess)," [QCD-scale = %.1f GeV]",scale); 
              plot_1(xMin+dx,xMax-dx,nPoints,f,NULL,mess,"x",p_name);  
           }
           put_text(&screen);
         } else messanykey(16,5," Correct input is \n"
                                "  0<=xMin<xMax<=1,\n"
                                " QCD-scale > 0.5 GeV\n"
                                " 3<=nPoints<=150");    
         break;    
      }
   }
}


static int  in_setting(void)
{
  int mode=1;
  void * pscr=NULL;    
  double Pcm; 
  void (*f7_tmp)(int)=f3_key[4];
  char * f7_mess_tmp= f3_mess[4];
  char sf_txt[STRSIZ];
  REAL mass[2];
  int i;
  int returnCode=0;

  if(nin_int == 1) return returnCode;

  for(i=0;i<2;i++) 
  if(sf_num[i])mass[i]=sf_mass[i];else pinf_int(Nsub,i+1,mass+i,NULL);
   
/* **   menu loop */

  for(;;)
  {    
    char strmen[] = "*"
    " S.F.1: First_structure_function                     "
    " S.F.2: Second_structure_function                    "
    " First  particle momentum[GeV] = PPP1                "
    " Second particle momentum[GeV] = PPP2                "
    " FirstPol                                            "
    " SecondPol                                           ";

    Pcm=va_int[0];

    strmen[0]=strlen(strmen)/6;         
    if(is_polarized(1,Nsub))
    improveStr(strmen,"FirstPol", "Helicity of  first particle   %.3G",Helicity[0]); 
    else improveStr(strmen,"FirstPol", "First  particle unpolarized");
    if(is_polarized(2,Nsub))
    improveStr(strmen,"SecondPol","Helicity of second particle   %.3G",Helicity[1]); 
    else improveStr(strmen,"SecondPol", "Second particle unpolarized");

    strFunName(1,sf_txt); improveStr(strmen,"First","%-45.45s", sf_txt);
    strFunName(2,sf_txt); improveStr(strmen,"Second","%-45.45s",sf_txt);
    improveStr(strmen,"PPP1","%-10.4G",inP1);
    improveStr(strmen,"PPP2","%-10.4G",inP2);

    f3_key[4]=f7_prog;   f3_mess[4]="Plot";
    menu1(25,7,"",strmen,"n_in_*",&pscr,&mode);
    f3_key[4]= f7_tmp;  f3_mess[4]=f7_mess_tmp;

    switch(mode)
    { case 0:
         for(i=0;i<2;i++) 
         if(sf_num[i])mass[i]=sf_mass[i];
         else pinf_int(Nsub,i+1,mass+i,NULL);         

         if((mass[0]==0 && inP1<=0)|| (mass[1]==0 && inP2<=0)) 
           messanykey(10,10,"For massless particle\nmomentum should be positive\n");
         else 
         {
           initStrFun(0);                                                                                                    
           return returnCode;
         } break;
      case 1:
      case 2: if(sf_menu(mode))
              { initStrFun(mode);            
                returnCode=returnCode|3;
              } 
              break;
      case 3: correctDouble(50,12,"Enter new value ",&inP1,1);
              returnCode=returnCode|1; break;
      case 4: correctDouble(50,12,"Enter new value ",&inP2,1); 
              returnCode=returnCode|1; break;
      case 5: if(is_polarized(1,Nsub))
              {  double buf=Helicity[0]; 
                 int spin2; 
                 char txt[60];
                  (*pinfAux_int)(Nsub,1, &spin2,NULL,NULL);
                  sprintf(txt, "Enter new value [%.1f,%.1f] :", -(spin2/2.),(spin2/2.));
                  correctDouble(40,12,txt,&buf,1);
                  if(fabs(2*buf)>spin2) 
                  { messanykey( 10,10,"Helicity out of limits");
                    if(blind) exit(111);
                  }else 
                  { Helicity[0]=buf;
                    returnCode=returnCode|1;
                  }   
              }   break;
      case 6: if(is_polarized(2,Nsub))
              {  double buf=Helicity[1];
                 int spin2; 
                 char txt[60];
                  (*pinfAux_int)(Nsub,2, &spin2,NULL,NULL);
                  sprintf(txt, "Enter new value [%.1f,%.1f] :", -(spin2/2.),(spin2/2.));
                  correctDouble(40,12,txt,&buf,1);
                  if(fabs(2*buf)>spin2) 
                  { messanykey( 10,10,"Helicity out of limits");
                    if(blind) exit(111);
                  }else 
                  { Helicity[1]=buf;
                    returnCode=returnCode|1;
                  }   
              }   break;                            
    }
  }
} /* in_setting */




static  int w_men__(void)
{
    int key =1;
    void * pscr=NULL;     
L1:
   {
      char strmen[] ="\030"
                      " BreitWigner range XXX  "
                      " T-channel widths   YYY "
                      " GI in t-channel    ZZZ "
                      " GI in s-channel    WWW ";   

      improveStr(strmen,"XXX","%.1f",*BWrange_int);
      improveStr(strmen,"YYY", (*twidth_int)? "ON" : "OFF");
      improveStr(strmen,"ZZZ",(*gtwidth_int)? "ON" : "OFF");
      improveStr(strmen,"WWW",(*gswidth_int)? "ON" : "OFF");

      menu1(54,6,"",strmen,"n_width_*",&pscr ,&key);
   }

   switch(key)
   { case 0:  return 0;
     case 1:  correctDouble(40,18,"Breit-Wigner range = ",BWrange_int ,1);
              if(*BWrange_int<0) *BWrange_int=-*BWrange_int;
              if(*BWrange_int>999.9) *BWrange_int=999.9;
              break;
     case 2:  *twidth_int = ! *twidth_int; break;
     case 3:  *gtwidth_int =! *gtwidth_int; break; 
     case 4:  *gswidth_int =! *gswidth_int; break;
   }    
    goto L1;
} /* w_men__ */

static int checkEnergy(void)
{  int i;
   REAL ms,m_;

   for(i=nin_int+1,ms=0; i<=nin_int+nout_int;i++)
   {  
       pinf_int(Nsub,i,&m_,NULL);
       ms+=m_;
   }
 
   if(nin_int==1)
   { 
      pinf_int(Nsub,1,&m_,NULL);
      if(m_<=ms) return 1;
   } else 
   {                                                         
     REAL  S,m1,m2;
     pinf_int(Nsub,1,&m1,NULL);
     pinf_int(Nsub,2,&m2,NULL); 
     incomkin(m1,m2,inP1,inP2,&S,NULL,NULL);
     if(S <= ms) return 1;
   }
   return 0;

}


static void infor(void)
{
  scrcolor(FGmain,BGmain);  
  clrbox(1,1,53,4);
  goto_xy(4,3); scrcolor(Red,BGmain);    print("(sub)Process: ");
  scrcolor(FGmain,BGmain); print("%s",Process);
  goto_xy(4,4); scrcolor(Red,BGmain);    print("Monte Carlo session: ");
  scrcolor(Black,BGmain);  print("%d",nSess);
  if(integral.old) 
  { print("(continue)");     
    if(integral.n_it>0)
    {
      goto_xy(1,7); scrcolor(Blue, BGmain);
      print(" #IT %s Error[%%]  nCall    Eff.  chi^2", nin_int == 2? "Cross section[pb]":"   Width[GeV]    ");
      goto_xy(1,8);scrcolor(FGmain, BGmain);
      integral.In=integral.s1/integral.n_it;               
      integral.dI=sqrt(integral.s0)/integral.n_it;                        
      if(integral.n_it<=1 || integral.s0==0 ) integral.khi2=0; else           
      integral.khi2=(integral.s2-integral.s1*integral.s1/integral.n_it)*integral.n_it/(integral.n_it-1)/fabs(integral.s0);
      print(" < >   %12.4E %10.2E %8d %8.8s %8.1G" ,
      integral.In, fabs(integral.In)? 100*integral.dI/(double)fabs(integral.In):0., integral.nCallTot, effInfo(),integral.khi2);
    }      
  }  
  else 
  {   print("(begin)");
      scrcolor(FGmain,BGmain);                 
      clrbox(1,8,53,maxRow()-2);      
  }
}


static int mappingMenu(void)
{
  int r=0;
  int mode=1;
  void * pscr=NULL;
  char menutxt[]="\030"
                 " Kinematics             "
                 " Regularization         ";
  for(;;)
  {
     menu1(54,6,"",menutxt,"n_map_*",&pscr, &mode);
     switch(mode)               
     { case  0: return r;
       case  1: r=r|entkin_();  break;
       case  2: do r=r|(2*edittable(1,4,&regTab,1,"n_reg",0));while(fillRegArray());
                break;    
      }
   }
}

static void f10_key_prog_for22(int x)
{
    if( mess_y_n(15,15," Quit session? ")) 
    {
        finish();
        sortie(0);
    }
}



int monte_carlo_menu(void)
{
   static int r=0;
   int mode=1;
   void * pscr=NULL;
   void * pscr_mem=NULL;
   void (*quit)(int)=f3_key[7];
   char menutxt[]="\030"
                  " Subprocess             "
                  " IN state               "
                  " Model parameters       "
                  " Constraints            "
                  " QCD coupling           "
                  " Breit-Wigner           "
	          " Aliases                "
	          " Cuts                   "
	          " Phase space mapping    "
                  " Monte Carlo simulation "
                  " Easy                   ";
                  
   if(nout_int!=2  ) menutxt[menutxt[0]*10+1]=0;
   if(nin_int==1)  improveStr(menutxt,"Easy", "Total width"); 
           else    improveStr(menutxt,"Easy", "1D integration");
 
   get_text(1,10,80,24,&pscr_mem);
   wrtprc_();
   for(;;)
   {  
      infor();
      f3_key[7]=quit;
      menu1(54,4,"",menutxt,"n_mc_*",&pscr, &mode);
      if(mode==1||mode==2||mode==3||mode==5||mode==7)f3_key[7]=NULL;

      switch (mode)
      { 
        case  0: put_text(&pscr_mem); return 0;
        case  1: r=r|3*sub_men__(); break;
        case  2: r=r|in_setting(); break;
        case  3: r=r|change_parameter(54,7,0);  break;
        case  4: { int modeC=1;
                   for(;;)
                   { char menuC[]="\030"
                     " All Constraints        " 
                     " Masses,Widths,Branching"; 
                     void * pscrC=NULL;
                     menu1(54,6,"",menuC,"n_constr_*",&pscrC, &modeC);
                     switch(modeC)
                     { case 0: put_text(&pscr_mem); break;
                       case 1: show_depend(54,7); break;
                       case 2: show_spectrum(54,9); break;
                     } 
                     if(!modeC) break;
                   } break;
                 }     
        case  5: r=r|qcdmen_();  break;
        case  6: r=r|w_men__();  break;
        case  7: do r=r|(3*edittable(1,4,&compTab,1,"n_comp",0)); while (fillCompositeArray());
	         break;     
        case  8: do r=r|(3*edittable(1,4,&cutTab,1,"n_cut",0)); while (fillCutArray()); 
                 break;             
        case  9: r=r|mappingMenu(); break;                      
        case  10: 
                 if(nout_int==1 && !sf_num[0] && !sf_num[1]  ) 
                 { if(blind)  return 1;                     
                   messanykey(15,15,"Phase space integration for 2->1 processes\n needs distribution functions.");
                   break;
                 }
                                                                        
                 
                 if(checkEnergy())   
                 { if(blind)  return 1;                  
                   messanykey(15,15,"Energy is too small!");                   
                   break;
                 }

                 if(fillCutArray()) 
                 { if(blind) return 2;
                   messanykey(15,15,"Can not evaluate cuts limlts"); 
                   break;
                 }  
        case 11:
                if(mode==11) 
                {  void (*f10_tmp)(int);
                   w_sess__(NULL);
                   f10_tmp=f3_key[7];
                   f3_key[7]=f10_key_prog_for22;
                   if(nin_int==1) decay12(); else
                   { REAL m1,m2, Pcm;
                     pinf_int(Nsub,1,&m1,NULL); 
                     pinf_int(Nsub,2,&m2,NULL);  
                     incomkin(m1,m2,inP1,inP2,NULL,&Pcm,NULL); 
                     if(sf_num[0]||sf_num[1]||nCuts)
                      messanykey(10,10,"Structure functions and cuts are ignored\n");                                       
                     cs_numcalc(Pcm);
                   }
                   f3_key[7]= f10_tmp;
                   r_sess__(NULL); 
                   break;
                } else if(fillRegArray()) 
                {  
                  if(blind) return 3;                
                   messanykey(15,15,
                       "Can not evaluate regularization paremeters");
                   break;    
                }
   
                if(mode==10)  runVegas(); 
                r=0;  
                break;
                 
      }
      if(r) clearEventMax();
      if(r&2) clearGrid();
      if(r&1)newSession();
   }
}
