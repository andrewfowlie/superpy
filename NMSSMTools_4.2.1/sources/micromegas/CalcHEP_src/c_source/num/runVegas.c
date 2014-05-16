/*
 Copyright (C) 2002, Alexander Pukhov
*/

#include<unistd.h>
#include<stdarg.h>
#include<math.h>
#include<stdlib.h>
#include<stdio.h>

#include "interface.h"
#include "cut.h"
#include "4_vector.h"
#include "q_kin.h"
#include "regul.h"
#include "runVegas.h"
#include "rw_sess.h"
#include "subproc.h"
#include "vegas.h"
#include "strfun.h"
#include "alphas2.h"
#include "crt_util.h"
#include "histogram.h"
#include "drandXX.h"
#include "n_calchep_.h"
#include "files.h"
#include "usrfun.h"

//========= old events.c

#include"chep_crt.h"
#include"../../include/version.h"

static FILE * events_;

static void writeEvent(long cCube,  int n, char * rand_state,double*pvect)
{ 
   int i;
   int icc;
   double qF,qR;
#ifdef PARKED 
 fprintf(events_,"%05X|%s|%d\n",cCube,rand_state,n);
#else

   if(*cb_pow_int)
   {  double sum=0;
      for(i=0;i<*cb_pow_int;i++) sum+=fabs((*cb_coeff_int)[i]);
      sum*=drandXX();
      for(i=0;i<*cb_pow_int;i++)
      { sum-=fabs((*cb_coeff_int)[i]);
        if(sum<=0) break;
      } 
      if(i==*cb_pow_int) i--;
      if((*cb_coeff_int)[i]<0) n*=-1;
      icc=i;
   }
   

   fprintf(events_,"%8d ",n);  
   if(nin_int==2) fprintf(events_," %17.10E %17.10E",pvect[3],pvect[7]);

   for(i=0;i<nout_int;i++) fprintf(events_," %17.10E %17.10E %17.10E",
   pvect[4*(i+nin_int)+1],pvect[4*(i+nin_int)+2],pvect[4*(i+nin_int)+3]);

   Scale(pvect,&qF,&qR);
   fprintf(events_,"| %.3E  %.3E ", qF,alpha_2(qR));

   
   if(*cb_pow_int)
   {  int j;
      fprintf(events_,"  ");
      for(j=0;j<*cb_nc_int;j++)
      fprintf(events_,"(%d %d)",(*cb_chains_int)[2*(*cb_nc_int)*icc+2*j], 
                                (*cb_chains_int)[2*(*cb_nc_int)*icc+2*j+1]);
   }


   fprintf(events_,"\n");

#endif   
}

//static long n_cube=1000;
#define n_cube EventGrid
extern long EventGrid;


static long nRandom=100;
static long nSimplex1=50;
static long nEvents=10000;
static double max =1.2;
static double milk=0.1; 

int saveEventSettings(FILE * f)
{
  fprintf(f, "%ld %ld %ld %ld ",n_cube, nRandom, nSimplex1,  nEvents);
  return 0;
}

int  readEventSettings(FILE * f)
{
  fscanf(f, "%ld %ld %ld %ld %ld ", &n_cube, &nRandom, &nSimplex1, &nEvents);
  return 0;
}


static void write_event_cap(void) 
{
  int i,j;

  fprintf(events_,"#%s\n",VERSION_);
  fprintf(events_,"#Type %d -> %d\n",nin_int,nout_int); 
  fprintf(events_,"#Initial_state\n");
  fprintf(events_,"  P1_3=%E" , inP1); 
  if(nin_int>1) fprintf(events_,"  P2_3=%E\n",-inP2);else fprintf(events_,"\n");
  if(nin_int>1) wrt_sf__(events_);

  fprintf(events_,"#PROCESS  ");
  for(i=1;i<=nin_int+nout_int; i++)
  { int pcode;
    char * pname=pinf_int(Nsub,i,NULL,&pcode);
    switch(pcode)
    { case  81: pcode= 1; break;
      case -81: pcode=-1; break;
      case  83: pcode= 3; break;
      case -83: pcode=-3; break;
    }
    fprintf(events_," %d(%s)", pcode, pname);
    if(i==nin_int)  fprintf(events_," ->");
  } 
  fprintf(events_,"\n");    
  fprintf(events_,"#MASSES ");
  for(i=0;i<nin_int+nout_int;i++)
  {  REAL m;
     pinf_int(Nsub,i+1,&m,NULL); 
     fprintf(events_," %.10E", (double)m);
  }   
  fprintf(events_,"\n");
  
  if(integral.n_it) fprintf(events_,"#Cross_section(Width) %E\n",integral.s1/integral.n_it);
  else   fprintf(events_,"#Cross_section(Width) Unknown\n"); 
  fprintf(events_,"#Number_of_events %10d\n",0);

  fprintf(events_,"#Events  "); 
  if(nin_int==2) fprintf(events_,"     P1_3 [Gev]        P2_3 [Gev]   ");
  for(i=1;i<=nout_int; i++) for(j=1;j<=3;j++) 
                          fprintf(events_,"     P%d_%d [Gev]   ",i+nin_int,j);

//  fprintf(events_,"  QCD SCALE    Color  !chains\n");
  fprintf(events_,"  Q_factor   alpha_QCD  Color chains\n");

}

static void improveEvents(vegasGrid * vegPtr,double (*func)(double *,double,double*))
{ int mode;
   void * pscr_=NULL;
   double eff0;

  for(mode=1;mode!=4; )
  {
      char strmen[]="\030"
                    " sub-cubes = N1         "
                    " random search = N2     "
                    " simplex search= N3     "
                    " Start search of maxima "; 
      improveStr(strmen,"N1","%d",n_cube);
      improveStr(strmen,"N2","%d",nRandom);
      improveStr(strmen,"N3","%d",nSimplex1);
      menu1(54,14,"Preparing of generator",strmen,"n_prep_gen_*",&pscr_,&mode);
      switch(mode)
      {
        case 0: return;
        case 1: if(correctLong(50,15,"Number of sub-cubes:",&n_cube,1))
                {
                    free(vegPtr->fMax); 
                    vegPtr->fMax=NULL;
                    vegPtr->nCubes=0;
                }  break;
        case 2: correctLong(50,15,"Random search:",&nRandom,1);break;
        case 3: correctLong(50,15,"Simplex steps :",&nSimplex1,1); break; 

/*                    " milk      = N3         "  
      improveStr(strmen,"N3","%.1f",milk);
        case 4: correctDouble(50,15,"Content of milk:",&milk,1); break;
*/
        case 4: if(n_cube < 1 ) n_cube=1; 
        { 
            int mCheck=vegas_max(vegPtr,n_cube, nRandom,nSimplex1, milk, func,&eff0);
            if(mCheck==0)
            { char mess[50];
              sprintf(mess,"Expected efficiency %f",eff0/max);
              messanykey(25,15,mess);
              put_text(&pscr_);
            } else
            {  mode=1;
               if(mCheck==2) messanykey(25,15,"Not enough memory.\n"
                                    "Decrease the number of sub-cubes");
            }                                                                                                         
        }                                                                             
      }
  }
}
   
static void  generateEvents( vegasGrid * vegPtr,  
                   double (*func)(double *,double,double*),  char *fname,  FILE * iprt)
{                                                   
   int mode=1;
   void * pscr=NULL;
   static int regen=1; 
//   if(!vegPtr->fMax)improveEvents(vegPtr,func);
    
   for(mode=1;;)
   {
     char strmen[]="\032"
                   " Number of events=N1      "
                   " Launch generator         "
		   " Allow weighted events OFF";
		   
     if(!regen) improveStr(strmen,"OFF","%s","ON");		   
     improveStr(strmen,"N1","%d",nEvents);
     menu1(53,10,"",strmen,"n_gen_*",&pscr,&mode);
     switch(mode)
     { case 0: return; 
       case 1: correctLong(50,15,"",&nEvents,1); break;
       case 2: 
       { long  nGenerated=0; 
         double eff;
         int nmax,mult,neg;   
         char mess[200];
         long cEvent;
         long fileEnd;
         
         if(!vegPtr->fMax) { messanykey(15,15,"Generator is not ready."); break;}
         events_= fopen(fname,"a");       
         if(ftell(events_)==0) write_event_cap();
         fileEnd=ftell(events_); 
         build_cb_int(Nsub);
         cEvent= vegas_events(vegPtr,nEvents,max,func,writeEvent,regen,&eff,&nmax,&mult,&neg);
         fclose(events_);

         if(cEvent>0)
         {  int l;
            sprintf(mess,"Statistic\n Events generated: %d\n  efficiency: %.1E\nMax event multiplicity: %d\n"
                      "Multiple events(total): %d \nNegative weight  events: %d \n", cEvent, eff,nmax,mult, neg);

            l=strlen(mess);
            strcat(mess,"---------------\n Accept events? ");    
            if(mess_y_n(25,15,mess)) 
            {  
               long  nEvPos=0;
               integral.old=1;
               mess[l]=0;
               events_=fopen(fname,"r+");
               while(nEvPos==0)
               { char ch;
                 char word[100];
                 do fscanf(events_,"%c",&ch); while(ch !='#');
                 fscanf(events_,"%s",word);
                 if(strcmp(word,"Number_of_events")==0) nEvPos=ftell(events_);
               }
               fscanf(events_,"%ld",&nGenerated);
               nGenerated+=cEvent;
               fseek(events_,nEvPos,SEEK_SET);
               fprintf(events_," %10d",nGenerated);
               fclose(events_);
               fprintf(iprt," %d events are stored in '%s'\n",nGenerated,fname);
               fprintf(iprt,"%s\n",mess);
               fflush(iprt);
            } else  truncate(fname,fileEnd);
         }
         destroy_cb_int();
       } 
       break;

//       case 3: improveEvents(vegPtr,func); put_text(&pscr); break;
       case 3: regen=!regen; 
       
//       goto ret;

     }
   }
}                                                                     


//======== old runVegas.c
int nSess=1;

double inP1=3500, inP2=3500;

static vegasGrid * veg_Ptr=NULL;
static int hFill=0;
vegas_integral integral={{5,5},{10000,10000},0,0.,0.,0.,0.,0.,0.,0,0,0,-1}; 


char * effInfo(void)
{ static char buff[10]; 
  int k;
  double sum;
  if( !integral.freeze || !integral.n_it || !integral.s1 || !veg_Ptr || !veg_Ptr->fMax  ) { buff[0]=0; return buff;}
  for(sum=0,k=0;k<veg_Ptr->nCubes;k++) sum+=veg_Ptr->fMax[k];
  sprintf(buff,"%.1E", integral.s1/integral.n_it*veg_Ptr->nCubes/sum/1.2);
  return buff; 
}

static void clearStatistics(int tp)
{
  integral.In=0;
  integral.dI=0;
  integral.khi2=0;
  integral.s0=0; 
  integral.s1=0; 
  integral.s2=0; 
  integral.n_it=0; 
  integral.nCallTot=0; 
  integral.tp=tp;
  clearHists();  
  { char fname[20];
    sprintf(fname,"distr_%d",nSess);
    unlink(fname);
  }
}                         

void clearGrid(void){ vegas_finish(veg_Ptr); veg_Ptr=NULL;}
void clearEventMax(void)
 { if(veg_Ptr && veg_Ptr->fMax) {free(veg_Ptr->fMax); veg_Ptr->fMax=NULL; veg_Ptr->nCubes=0; }}

void newSession(void)
{
   if(integral.old)
   { char fname[20];

     messanykey(15,15,
     "Some parameters where changed.\nSo integral and statictics for\n"
     "distribushions is forgotten!\nSession number is increased.");
     integral.old=0;
     nSess++;
     clearStatistics(-1);
     sprintf(fname,"prt_%d",nSess);
     unlink(fname);
   }
}



int saveVegasGrid( FILE * f)
{
  if(veg_Ptr)
  {  int i,j;
     double * x=veg_Ptr->x_grid;
     fprintf(f," Vegas_grid: dim=%d  size=%d\n", veg_Ptr->ndim, veg_Ptr->ndmx);
     for(i=0;i<veg_Ptr->ndim;i++)
     { for(j=0;j<=veg_Ptr->ndmx;j++) fprintf(f," %.15E",*(x++));
       fprintf(f,"\n");
     }
     if(veg_Ptr->fMax)
     { long l;
       fprintf(f,"Max(%d):\n",veg_Ptr->nCubes);
       for(l=0;l<veg_Ptr->nCubes;l++) fprintf(f,"%.1E\n",veg_Ptr->fMax[l]);
     } else fprintf(f,"Max(0):\n");
  }else  fprintf(f," Vegas_grid: dim=%d  size=%d\n", 0, 0);
  return 0;
}

int readVegasGrid(FILE * f)
{
  int i,j,ndim,ndmx;
  double * x;
  
  if(veg_Ptr) {vegas_finish(veg_Ptr);veg_Ptr=NULL;}  
  fscanf(f," Vegas_grid: dim=%d  size=%d\n", &ndim, &ndmx);
  if(ndim && ndmx)
  { 
    veg_Ptr=vegas_init(ndim,ndmx);
    x=veg_Ptr->x_grid;
    for(i=0;i<ndim;i++)for(j=0;j<=ndmx;j++) fscanf(f," %lf",(x++));
    fscanf(f," Max(%ld):\n",&(veg_Ptr->nCubes));       
    if(veg_Ptr->nCubes) 
    { long l;
      veg_Ptr->fMax=malloc(sizeof(float)*veg_Ptr->nCubes);
      for(l=0;l<veg_Ptr->nCubes;l++) fscanf(f,"%f",veg_Ptr->fMax+l);
    } else veg_Ptr->fMax=NULL;
  }
  return 0;
}

static int nCall;
static double badPoints;
static double negPoints;

static void printLn(FILE * iprt,int *line,char * format, ...)
{  
   va_list args;
   char dump[STRSIZ];
   va_start(args, format);
   vsprintf(dump,format,args);
   va_end(args);

   goto_xy(1,*line); print("%53s","");
   goto_xy(1,*line);
   print("%s\n",dump); 
   (*line)++;
   if (*line >= maxRow()-2 ) *line=8; else 
   {
     scrcolor(Blue, BGmain);
     print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
   }   
   if(iprt) {fprintf(iprt,"%s\n",dump); fflush(iprt);}
}

static double func_(double *x, double wgt)
{
    double ret_val=0.;
    int err=0,nd,i;
    double factor_0;
    double xx0=x[0],xx1=x[1],x1,x2;    
    double GG,qF,qR;
    nCall++; 
    REAL pvectR[100];
    double pvect[100];
/* get momenta */
    mkmom(x, &factor_0,pvectR);
    if(sf_num[0]){x1=x[0]; x[0]=xx0;}
    if(sf_num[1]){x2=x[1]; x[1]=xx1;}
    if (!factor_0) goto exi;
    nd=4*(nin_int+nout_int);
    for(i=0;i<nd;i++) pvect[i]=pvectR[i];
    factor_0 *= calcCutFactor(pvect)*usrFF(nin_int,nout_int,pvect,p_names,p_codes); 
    if (!factor_0)   goto exi;


    Scale(pvect,&qF,&qR);
/* **  structure function  multiplication */
    if (nin_int == 2) 
    {
	if(sf_num[0]) { factor_0 *= strfun_(1, x1,qF);  if(factor_0==0.) {/*printf("|x1=%.2f|",x1);*/ goto exi;}}
	if(sf_num[1]) { factor_0 *= strfun_(2, x2,qF);  if(factor_0==0.) {/*printf("|x2=%.2f|",x2);*/ goto exi;}} 
    }   
    if (!factor_0)  { goto exi;}
/* ** call for 'running strong coupling constant' */
    GG=sqrt(4*M_PI*alpha_2(qR));    
    ret_val = factor_0 * sqme_int(Nsub,GG,pvectR,&err);

    if(err)       badPoints+=  (ret_val>0 ? ret_val*wgt : - ret_val*wgt); 
    if(ret_val<0) negPoints+=ret_val*wgt;
exi:

    if(hFill) fillHists(ret_val*wgt,pvect);

    return ret_val;
} /* func_ */

static double func_2(double *x, double wgt,double *pout)
{
    double ret_val=0.;
    int err=0,nd,i;
    double factor_0;
    double xx0=x[0],xx1=x[1],x1,x2;    
    double GG,qF,qR;
    nCall++; 
    REAL pvectR[200];
    double pvect[40];
/* get momenta */
    mkmom(x, &factor_0,pvectR);

    if(sf_num[0]){x1=x[0]; x[0]=xx0;}
    if(sf_num[1]){x2=x[1]; x[1]=xx1;}
    
    nd=0;
    if (!factor_0) goto exi;

    nd=4*(nin_int+nout_int);
    for(i=0;i<nd;i++) pvect[i]=pvectR[i];
    factor_0 *= calcCutFactor(pvect);
    if(nin_int>1) factor_0 *= usrFF(nin_int, nout_int,pvect,p_names,p_codes); 
    if (!factor_0)   goto exi;

    Scale(pvect,&qF,&qR);
/* **  structure function  multiplication */
    if (nin_int == 2) 
    {
	if(sf_num[0]) { factor_0 *= strfun_(1, x1,qF);  if(factor_0==0.) {/*printf("|x1=%.2f|",x1);*/ goto exi;}}
	if(sf_num[1]) { factor_0 *= strfun_(2, x2,qF);  if(factor_0==0.) {/*printf("|x2=%.2f|",x2);*/ goto exi;}} 
    }   
    if (!factor_0)  { printf("strf");  goto exi;}
/* ** call for 'running strong coupling constant' */
    GG=sqrt(4*M_PI*alpha_2(qR));
    
    ret_val = factor_0 * sqme_int(Nsub,GG,pvectR,&err);
    
    if(err)       badPoints+=  (ret_val>0 ? ret_val*wgt : - ret_val*wgt); 
    if(ret_val<0) negPoints+=ret_val*wgt;
exi:
    if(pout)for(i=0;i<nd;i++) pout[i]=pvect[i];
    
    return ret_val;
} /* func_ */


extern long EventGrid;

int runVegas(void)
{
    int i;
    double sd;
    double avgi;
    char mess[25];
    FILE * iprt = NULL;
    int mode=1;
    void * pscr=NULL;
    static int n_Line=7;

    i=imkmom(inP1,inP2);
    if(veg_Ptr&&veg_Ptr->ndim!=i)clearGrid();
    if(!veg_Ptr) veg_Ptr=vegas_init(i,50);     

    if(nin_int == 2) strcpy(mess, "Cross section[pb]");
      else           strcpy(mess, "   Width[Gev]    ");
    
/* ** save current session parameters */
     w_sess__(NULL);
/* ** open protocol and resulting files */
       
    {  char fname[50];
       sprintf(fname,"%sprt_%d",outputDir,nSess);
       iprt=fopen(fname,"a");
       if(ftell(iprt)==0) 
       { fprintf(iprt,"    CalcHEP kinematics module \n The session parameters:\n");
         w_sess__(iprt);
         fprintf(iprt,"===================================\n");   
       }
    }

/* **  initkinematics */

    while(correctHistList()) editHist();
    
/* *** Main cycle */
    if(!integral.old || n_Line==7)
    { n_Line=7;
      scrcolor(Blue, BGmain);
//    printLn(iprt,&n_Line," #IT  %20s Error %%    nCall   chi**2",mess); 
      printLn(iprt,&n_Line," #IT %s Error[%%]  nCalls   Eff.  chi^2",mess);
    }

    for(;;)
    {
        static int worn=1;
        char strmen[]="\030"
         " nSess  = N2_1          "
         " nCalls = N1_1          "
         " Set  Distributions     "
         "*Start integration      "
         " Display Distributions  "
         " Clear statistic        "
         " Freeze grid        OFF " 
	 " Clear  grid            "
	 " Event Cubes NCUBE      "
	 " Generate Events        ";

        improveStr(strmen,"N1_1","%d",integral.ncall[0]);
        improveStr(strmen,"N2_1","%d",integral.itmx[0]);
        improveStr(strmen,"NCUBE","%d",EventGrid);

        if(integral.freeze) improveStr(strmen,"OFF","ON");

        menu1(54,7,"",strmen,"n_veg_*",&pscr,&mode);
        switch(mode)
        {     
        case 0:
          if(iprt) fclose(iprt);
          return 0;           
        case 1:  
          correctInt(50,12,"Enter new value ",&integral.itmx[0],1); break;
        case 2: 
          correctLong(50,12,"Enter new value ",&integral.ncall[0],1); break;
        case 3:  editHist(); break;
        case 4:
          if(veg_Ptr->fMax && !integral.freeze)
          {  if(!mess_y_n(15,15,"You have event generator prepared.\n"
             " The  answer 'Y'  will start Vegas session \nwhich destroys it."
             " To save the event generator answer 'N' \nand set "
             " ' Freeze grid' ON")) break;
             else { free(veg_Ptr->fMax); veg_Ptr->fMax=NULL; veg_Ptr->nCubes=0;}  
          }
          for (i = 1; i <= integral.itmx[0]; ++i)                                       
          { double sum;
            char  errtxt[100]="";
            int k;

            if(integral.ncall[0]==0) break;
            nCall=0;                                                                  
            negPoints=0;                                                              
            badPoints=0; 
            hFill=1;  
            if(vegas_int(veg_Ptr, integral.ncall[0],1.5*(!integral.freeze), 
                 func_, &avgi, &sd)        
              ) break;
            integral.old=1;                                              
            negPoints/=nCall;                                                         
            badPoints/=nCall;                                                         
            integral.nCallTot+=nCall;                                                          
            scrcolor(FGmain,BGmain);                                                 
            printLn(iprt,&n_Line,"%4d   %12.4E %10.2E %8d %s",                     
                 ++integral.n_it, avgi,avgi? 100*sd/(double)fabs(avgi):0.,nCall,effInfo());
                                                                   
            if(negPoints<0) sprintf(errtxt+strlen(errtxt)," Negative points %.1G%%;",                
                                      -100*negPoints/(avgi-2*negPoints));             
            if(badPoints)  sprintf(errtxt+strlen(errtxt),                             
                 "Bad Precision %.1G%%;",100*badPoints/(avgi-2*negPoints));           
                                                                                      
            if(errtxt[0])                                                             
            {                                                                         
               scrcolor(Red,BGmain);                                                  
               printLn(iprt,&n_Line,"%s",errtxt);                                     
            }
                                                                                                                                                
            integral.s0+=sd*sd;                                                                  
            integral.s1+=avgi;                                                             
            integral.s2+=avgi*avgi;                                    
          } 
          
          
          integral.In=integral.s1/integral.n_it; 
          integral.dI=sqrt(integral.s0)/integral.n_it;
          if(integral.n_it<=1 || integral.s0==0 ) integral.khi2=0; else 
          integral.khi2=(integral.s2-integral.s1*integral.s1/integral.n_it)*integral.n_it/(integral.n_it-1)/fabs(integral.s0);  
          
          scrcolor(FGmain,BGmain);

          printLn(iprt,&n_Line," < >   %12.4E %10.2E %8d %7.7s %-7.1G" ,
                      integral.In, fabs(integral.In)? 100*integral.dI/(double)fabs(integral.In):0., integral.nCallTot, 
                                                              effInfo(),  integral.khi2);
          if(histTab.strings)
          { char  fname[20];
            FILE * d;
            sprintf(fname,"distr_%d",nSess);
            d=fopen(fname,"w");  
            wrt_hist2(d,Process);
            fclose(d);
          }
                    messanykey(54,11,"Integration is over");
/*          integral.freeze=0; */
          break;

        case 5: showHist(54,10,Process); break;
        case 6: clearStatistics(-1);
                messanykey(54,13,"Old results for integral\n"
                "and distributions\nare deleted.");
                break;
        case 7: integral.freeze=!integral.freeze; break; 
        case 8: if(!integral.freeze || mess_y_n(15,15,"The information for Event Generator will be lost\n OK?"))  
                { int ndim=veg_Ptr->ndim;
                  vegas_finish(veg_Ptr);
                  veg_Ptr=vegas_init(ndim,50);
                  messanykey(57,11,"OK");
                }   
                break;
        case 9: 
           if(correctLong(50,12,"Enter new value ",&EventGrid,1))
           { if(veg_Ptr->fMax) {free(veg_Ptr->fMax); veg_Ptr->fMax=NULL;veg_Ptr->nCubes=0;}} break;
        case 10: 
           if( !veg_Ptr || !veg_Ptr->fMax)
           { char * mess="Before event generation one has to launch  Vegas session with freezed grid\n"
                                           "to prepare generator";
                if(blind) { printf("%s\n",mess); sortie(200);}  else messanykey(4,13,mess);
           }    else    runEvents(); 
       }
    }    
}


int runEvents(void)
{
    FILE * iprt = NULL;
    int i;

    i=imkmom(inP1,inP2);
//    if(veg_Ptr&&veg_Ptr->ndim!=i)clearGrid();
//    if(!veg_Ptr) veg_Ptr=vegas_init(i,50);    

    w_sess__(NULL);
/* ** open protocol and resulting files */
       
    {  char fname[50];
       sprintf(fname,"%sprt_%d",outputDir,nSess);
       iprt=fopen(fname,"a");
       if(ftell(iprt)==0) 
       { fprintf(iprt,"    CalcHEP kinematics module \n The session parameters:\n");
         w_sess__(iprt);
         fprintf(iprt,"===================================\n");   
       }
    }

/* **  initkinematics */


    { char fname[50];

      sprintf(fname,"%sevents_%d.txt",outputDir,nSess);

      hFill=0;
      generateEvents(veg_Ptr,func_2,fname, iprt);
    }
    fclose(iprt);
    return 0;
}
