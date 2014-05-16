#include <stdlib.h>
#include <unistd.h>
#include <dlfcn.h>
#include <math.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <sys/stat.h>
             
#include "VandP.h"
#include "chep_crt.h"
#include "files.h"
#include "plot.h"
#include "dynamic_cs.h"
#include "vp.h"
#include "read_mdl.h"
#include "model.h"
#include "rootDir.h"
#include "simpson.h"
#include "alphas2.h"
#include "screen.h"
#include "viewdir.h"
#include "SLHAplus.h"
#include "physics.h"
#include "sos.h"
#include "numcheck.h"



static int findParam(int X, int Y, int rdf, int shift, int tot, char*title,int*pos)
{
   int  n;
   char*menu;
   char *hlp="";
   
   menu=malloc(2+24*(tot+1)); 
   menu[0]=24; menu[1]=0;
   if(rdf)
   {  sprintf(menu+1,"    READ_FROM_FILE      "); 
      hlp="change_var";
   }   
   for(n=0;n<tot;n++) sprintf(menu+strlen(menu), " %-11.11s%12.4E",
                      varNames[n+shift], (double)(varValues[n+shift]));                     
  
   menu1(X,Y,title,menu,hlp,NULL,pos);   
   free(menu);
   if(pos==0) return -1; else return *pos-1+shift - rdf ;
}

static int changeParam(int X,int Y) 
{ int n,ch;
  char txt[50];
  static char fName[100]="";
  double x;
  int pos;
  
  for(n=0,ch=0,pos=1;n>=-1;)
  { n=findParam(X-1,Y,1,0,nModelVars,"Change Parameter",&pos);
    if(pos==1)
    { FILE *f;
      struct stat buf;
      if(!findCalcHEPfile(fName)) continue;
      if(stat(fName,&buf) ||   !(S_ISREG(buf.st_mode)) )  { messanykey(10,17, "Not a regular file"); continue; }
      f=fopen(fName,"r");
      if(f==NULL) { messanykey(10,17, "Can not open file"); continue; }
      for(;;)
      {   char name[20];
          char txt[40];
          double val; 
          int i;
           
          if(fscanf(f,"%s",name)!=1) break;
          if(name[0]=='#') { fscanf(f,"%*[^\n]"); continue;}
          for(i=0;i<nModelVars;i++) if(strcmp(name,varNames[i])==0) break;
          if(i==nModelVars)
          {
             sprintf(txt,"'%s' - unknown variable",name);  
              messanykey(10,10,txt);
          }
          if(fscanf(f,"%lf",&val)!=1)
          { sprintf(txt," wrong defined number for '%s'",name);
            messanykey(10,10,txt);
          } else  if(i<nModelVars)  { varValues[i]=val; ch=1;}
          
          fscanf(f,"%*[^\n]");
      }  
      fclose(f);
           
    }else  if(n>=0)
    { x=varValues[n];
      sprintf(txt,"%s = ",varNames[n]);
      if(correctDouble(20,20,txt,&x,1)) {  varValues[n]=x; ch=1; }
    }
  }
  return ch;
}

static void show_dependence(int X, int Y)
{ void *pscr1=NULL;
  int i,mPos=1; 
  REAL mem;
  int nc,ni,pos1,pos2;
  char txt[50];
  for(pos1=1;;)
  { nc=findParam(X-1,Y,0,nModelVars,nModelFunc,"Constraint",&pos1);
    if(!pos1) return;
    for(pos2=1;;)
    { double xMin,xMax;
      int nPoints=100;
      sprintf(txt,"check \"%s\" depends on",varNames[nc]); 
      ni=findParam(X-1,Y,0,0,nModelVars,txt,&pos2);
      mem=varValues[ni];
      if(ni<0) break; 
      
      xMin=varValues[ni] - fabs(varValues[ni] )/10;
      xMax=varValues[ni] + fabs(varValues[ni] )/10;
      
      for(;;)
      {  int k3=0; 
         char strmen[]="\026 "
            " x-Min = XXX          "
            " x-Max = YYY          "
            " Npoints = NNN        "
            " Display              ";

         improveStr(strmen,"XXX","%G",xMin);
         improveStr(strmen,"YYY","%G",xMax);
         improveStr(strmen,"NNN","%d",nPoints);
         sprintf(txt,"check %s(%s)",varNames[nc],varNames[ni]);        
         menu1(X,Y+2,txt,strmen,"",NULL,&k3);
         if(!k3) break;
         switch(k3)
         {  case 1: correctDouble(X,Y+12,"xMin = ",&xMin,1); break;
            case 2: correctDouble(X,Y+12,"xMax = ",&xMax,1); break;
            case 3: correctInt(X,Y+12,"nPoints = ",&nPoints,1); break;
            case 4:
            if( xMax>xMin && nPoints>=3 && nPoints<=150)
            {  double dx=(xMax-xMin)/(nPoints-1);
               double f[150];
               int i, NaN=0,Esc=0;
         
               informline(0,nPoints);               
               for(i=0;i<nPoints;i++)
               {  double x=xMin+i*dx;
                  varValues[ni]=x;
                  NaN=calcMainFunc();
                  if(NaN) 
                  {  char mess[100];
                     sprintf(mess,"Can not evaluate constraints for %s=%G",varNames[ni], x);
                     messanykey(16,5,mess);        
                     break;
                  }
                  f[i]=varValues[nc];
                  Esc=informline(i,nPoints);
                  if(Esc) break;  
               }
                  
               varValues[ni]=mem;
               calcMainFunc();

               if(!(NaN||Esc)) plot_1(xMin,xMax,nPoints,f,NULL,"Plot",
                              varNames[ni], varNames[nc]);
                               
            } else messanykey(16,5," Correct input is \n"
                                   "  xMin<xMax,\n"
                                   " 3<=nPoints<=150");
            break;
         }
       }
     }
  }
}


static void  writeSLHA(void)
{  int i;
   FILE *f;
   char fName[100];

   for(i=1;;i++)
   { sprintf(fName,"decaySLHA%d.txt",i);
     if(access(fName,R_OK)) break;     
   }
       
   f=fopen(fName,"w");
   
   fprintf(f,"BLOCK ModelParameters # %s\n",currentModelName());
   for(i=0;i<nModelVars;i++)
   fprintf(f," %3d  %16E # %s\n",i+1, (double)varValues[i], varNames[i]);    
   fprintf(f,"#\n");

   for(i=0;i<nModelParticles;i++)
   {  
    fprintf(f,"BLOCK QNUMBERS %d  # %s\n", ModelPrtcls[i].NPDG, ModelPrtcls[i].name);   
    fprintf(f," 1  %d # 3*el.charge\n 2  %d # 2*spin+1\n 3  %d # color dim\n 4  %d # 0={ self-conjugated}\n#\n",
         ModelPrtcls[i].q3, 
         ModelPrtcls[i].spin2+1, 
         ModelPrtcls[i].cdim, 
         strcmp(ModelPrtcls[i].name,ModelPrtcls[i].aname)? 1:0 );
   }

   fprintf(f,"BLOCK MASS\n");   
   for(i=0;i<nModelParticles;i++) 
   { char *name=ModelPrtcls[i].name;
     fprintf(f," %d  %E # %s\n",ModelPrtcls[i].NPDG,pMass(name),name);
   }
   fprintf(f,"#\n");

   for(i=0;i<nModelParticles;i++)
   {  txtList all=NULL;
      double mass,width;
      char *name;
      
      if( strcmp(ModelPrtcls[i].mass,"0")==0) continue;
      if( strcmp(ModelPrtcls[i].width,"0")==0) continue;
      name=ModelPrtcls[i].name;
      mass=pMass(name);
      if(!mass) continue;
      width=pWidth(name,&all);
      fprintf(f,"DECAY %d  %E # %s\n",ModelPrtcls[i].NPDG,width,name);
      for(;all;all=all->next)
      { int dim; 
        char pn[20], buff[100], *chB,*chE;
        strcpy(buff,all->txt);

        chB=strstr(buff,"->");
        chB+=2;
        for(dim=1;;dim++)
        { 
           chE=strchr(chB,',');
           if(chE)chB=chE+1;else break;
        }

        sscanf(buff,"%s", pn);
        fprintf(f," %s %d  ",pn,dim);
        chB=strstr(buff,"->");
        chB+=2;
        for(;;)
        { 
           chE=strchr(chB,',');
           if(chE)chE[0]=0;
           sscanf(chB,"%s",pn);
           fprintf(f," %d", qNumbers(pn,NULL,NULL,NULL));
           if(chE)chB=chE+1;else break;
        }
        
        
        chB=strstr(all->txt,"->");   
        fprintf(f,"  # %s \n",chB+2);
      }
      fprintf(f,"#\n");          
   }
   fclose(f);
   { char buff[100];
     sprintf(buff,"See results in file '%s'", fName);
     messanykey(16,5,buff); 
   }                  
}


static void show_spectrum(int X, int Y)
{ int i;
  char *menuP=malloc(2+22*(nModelParticles+1));
  int mode=1;  
  menuP[0]=22;
  menuP[1]=0;

  strcpy(menuP+1," All Particles -> SLHA");
       
  for(i=0;i<nModelParticles;i++)
  { char *mass=ModelPrtcls[i].mass;
    char *name=ModelPrtcls[i].name;
    if(!strcmp(mass,"0")) sprintf(menuP+strlen(menuP)," %-6.6s      Zero     ",name);
    else sprintf(menuP+strlen(menuP)," %-6.6s  %12.4E ",name,pMass(name));
  }
  
  while(mode)
  {  menu1(X,Y,"",menuP,"n_qnumbers",NULL, &mode);
     if(mode==1) writeSLHA(); 
     else if(mode>1)  
     { FILE*f=fopen("width.tmp","w");
       int pos=mode-2;
       txtList LL=NULL;
       
       char *mass=ModelPrtcls[pos].mass;
        
       fprintf(f, "Patricle %s(%s),  PDG = %d,  Mass= ", 
       ModelPrtcls[pos].name, ModelPrtcls[pos].aname,
       ModelPrtcls[pos].NPDG);
       if(strcmp(mass,"0")==0)  fprintf(f, "Zero\n"); else
       { 
         double width;
         fprintf(f,"%.3E ", pMass(ModelPrtcls[pos].name));       
         width=pWidth(ModelPrtcls[pos].name,&LL); 
         fprintf(f," Width=%.2E\n",width);
       }
       fprintf(f,"Quantum numbers: ");
       { int spin=ModelPrtcls[pos].spin2;
         int q3=ModelPrtcls[pos].q3; 
         fprintf(f," spin=");
         if(spin&1)  fprintf(f,"%d/2, ",spin); else fprintf(f,"%d, ",spin/2);
         fprintf(f," charge(el.)="); 
         if(q3!=3*(q3/3)) fprintf(f,"%d/3, ",q3);else fprintf(f,"%d ",q3/3); 
         fprintf(f," color=%d\n",ModelPrtcls[pos].cdim);
       }
       if(LL) 
       { txtList ll=LL;
          fprintf(f," Branchings & Decay channels:\n");
          for(;ll;ll=ll->next)
          { char buff[100];
            double br;
            sscanf(ll->txt,"%lf %[^\n]",&br,buff);
            fprintf(f," %.2E     %s\n",br,buff);
          }   
       }  
  
       fclose(f);
       f=fopen("width.tmp","r");
       showtext(2,Y,78,  maxRow()-1,"Particle information",f); 
       fclose(f);                   
       unlink("width.tmp");       
     }
  }   
  free(menuP);
  return;
}


static void localF6(int x){  viewDir("./");}


//#define SO "./results/aux/so_generated/VandP.so"

int numcheck(void)
{  int err,size=100;
     
   for(;;)
   {  compDir=realloc(compDir,size+20);
      if(getcwd(compDir,size)) break; else size*=2;
   }
   modelDir="models";
   modelNum=n_model; 
   strcat(compDir,"/results/aux");
   libDir=malloc(strlen(compDir)+20);
   sprintf(libDir,"%s/so_generated",compDir);
   calchepDir= rootDir;
   if(prepareWorkPlace()) mkdir(libDir,00755); else if(checkWorkPlace())
   { char*command=malloc(strlen(libDir)+20);
      sprintf(command,"rm %s/*.so",libDir);
      system(command);
      free(command);
   }  

   if(getDynamicVP()) return 1;
   cleanDecayTable();
   ForceUG=forceUG;
   {
     char  mmenu[]="\026" 
                   " Parameters           "
                   " All Constraints      "
                   " Masses,Widths,Branch.";
     int m0=1;
     void (*F10)(int);
     void (*F6)(int);
     F10=f3_key[7];    F6=f3_key[3];
     f3_key[3]=localF6;
chdir("results");
outputDir="./"; 
     err=calcMainFunc();
     if(Warnings) messanykey(5,10,Warnings);
     for(;m0;)
     {  menu1(56,7,"",mmenu,"s_num_*",NULL,&m0);
       if(err && (m0==2||m0==3)) 
       { char txt[100];
         sprintf(txt," Can not calculate %s ",varNames[err]);
         messanykey(12,15,txt);
       }   
       switch(m0)
       {                
         case 1: if(changeParam(56,8))
                 {
                   cleanDecayTable();
                   err=calcMainFunc();  
                   if(Warnings) messanykey(5,10,Warnings);        
                 }
                 break;
         case 2: if(nModelFunc) {if(!err ) show_dependence(56,8);} 
                 else messanykey(5,10,"There are no public constraints in this model.");  
                 break;           
         case 3: if(!err) show_spectrum(56,8); break;
       }
     }
chdir("..");
outputDir="results/";
       f3_key[7]=F10;   
       f3_key[3]=F6;       
   }
   return 0;
}


#define  MCHARM   1.3

static double L3=3.125347E-01 , L4=2.763267E-01, L5=1.991586E-01, L6=8.449407E-02;
static int   nf3=3,nf4=4,nf5=5,nf6=6;
static double MbMb,Mtp;


static double alpha(int nf, int odr, double lambda,  double dscale)
{
    double b0 = 11. -  (2./3.)*nf;
    double b1 = 51. - (19./3.)*nf;
    double b2 = 2857. - (5033./9.)*nf + (325./27.)*nf*nf;
    double rl = 2*log(dscale / lambda);
    double alpha0= 4*M_PI/(b0*rl);
    double d__4 = log(rl) - .5;
    double d__2 = 2*b1/(b0*b0*rl);


    if(odr==1) return alpha0;
    else if(odr==2) return  alpha0*(1 - 2*b1*log(rl)/(b0*b0*rl));
    else if(odr==3) return  alpha0*(1 - 2*b1*log(rl)/(b0*b0*rl)
         + d__2*d__2 *(d__4*d__4 + b2*b0 /(8*b1*b1) - 1.25)  );
    else { fprintf(stderr,"Can not evaluate alpha in so large oder (%d).\n",odr);
           exit(1);
         }
}



double alpha_2(double Q)
{

  if(Q<L3*1.5) Q=L3*1.5;
  
        if(Q<MCHARM)  return alpha(nf3, 3, L3, Q);
  else  if(Q<MbMb)    return alpha(nf4, 3, L4, Q);
  else  if(Q<Mtp)     return alpha(nf5, 3, L5, Q);
  else                return alpha(nf6, 3, L6, Q);
}

