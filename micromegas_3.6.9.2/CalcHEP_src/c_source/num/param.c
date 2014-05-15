#include <sys/wait.h>
#include <sys/types.h>
#include <sys/stat.h>
#include"chep_crt.h"
#include"interface.h"
#include"plot.h"
#include"param.h"
#include"rw_sess.h"
#include"subproc.h"
#include"num_in.h"

#include"VandP.h"
#include"dynamic_cs.h"
#include"files.h"
#include"SLHAplus.h"

REAL Pcm22;

/* ************************************************* */
/* Physics model parameters menu                     */
/* ************************************************* */

static void exportParam(void)
{ static REAL ** link=NULL;
  int i,j;

  if(link==NULL)
  { 
    link=malloc(nvar_int*sizeof(REAL *));
    for(i=0;i<nvar_int;i++) link[i]=NULL;
    for(i=0;i<nvar_int;i++)
    { for(j=0;j<nModelVars+nModelFunc;j++)
      if(strcmp(varNames[j],varName_int[i+1])==0) {link[i]=varValues+j;break;}  
    }
  }

  for(i=0;i<nvar_int;i++) if(link[i])va_int[i+1]=*link[i]; 
  else if( strcmp("i", varName_int[i+1]))  printf(" '%s' not found\n",varName_int[i+1]);
  cleanDecayTable();   
}

int checkParam(void)
{ int err;


  err=calcMainFunc();

  if(err>0)
  {  char mess[100];
     sprintf(mess,"Can not evaluate constrained parameter '%s' ",
             varNames[err]);
    
     if(blind) { printf("%s\n",mess); sortie(122);} 
     else     messanykey(10,10, mess);
     return err; 
  }
  exportParam();
  { int err=calcFunc_int();
    if(err>0)
    {  char mess[100];
       sprintf(mess,"Can not evaluate constrained parameter '%s' ",
             varName_int[err]);
                                                                                
       if(blind) { printf("%s\n",mess); sortie(122);}
       else     messanykey(10,10, mess);
      return err;
    }
  }
                        
  return 0;
}

static void param_menu(int sqtrS_on,int rd_on, int vars_on, int func_on,  int * polar, char ** strmen)
{
  int k,pos,npos;
  npos=0;
  if(sqtrS_on)  
  {   npos++; for(k=0;k<2;k++) if(is_polarized(k+1,Nsub))
      {polar[k]=1; npos++;} else polar[k]=0;
  }  else {polar[0]=0;polar[1]=0;}
  
  if(vars_on) npos += nModelVars;
  if(func_on) npos += nModelFunc;
  if(rd_on)   npos++;   
  if(npos==0) {*strmen=NULL; return;}
  
  *strmen=malloc(24*npos+2);
  (*strmen)[0]=24;

  pos=1;    
  for(k=0;k<2;k++) if(polar[k])
  { sprintf((*strmen)+pos," Helicity%d  %-12.4g",k+1,Helicity[k]);
    pos+=24;               
  }
  if(sqtrS_on){ sprintf((*strmen)+pos,"  Pcm       %-12.4g",(double)Pcm22); pos+=24;}
  if(rd_on)
  { sprintf((*strmen)+pos,"%24.24s","     READ_FROM_FILE     ");
    pos+=24;
  }
  if(vars_on) for (k =0; k < nModelVars;k++)
  { 
     sprintf((*strmen)+pos," %-11s%12.4E",varNames[k],(double)varValues[k]);
     pos+=24;
  } 

  if(func_on) for (k =nModelVars; k<nModelVars + nModelFunc;k++)
  {  
     sprintf((*strmen)+pos," %-11s%12.4E",varNames[k],(double)varValues[k]);
     pos+=24;
  } 

  (*strmen)[pos]=0;    
}

int selectParam(int x, int y, char * mess, void ** pscrPtr,  
    int sqtrS_on, int rd_on, int vars_on, int func_on, 
     REAL ** varPos, char * varName,int*mPos)
{
    char* strmen;
    void * pscr;
    int polar[2];
    int position=1;
    
    if(pscrPtr) pscr=*pscrPtr; else pscr=NULL; 

    param_menu( sqtrS_on, rd_on, vars_on,func_on, polar,  &strmen);
    if(strmen)
    { int i, shift=0;
      char *hlp ="";
      if(rd_on) hlp="change_var";
      menu1(x,y,mess,strmen,hlp,&pscr,mPos);
      position=*mPos;
      if(position) sscanf(strmen+1+(position-1)*strmen[0],"%s",varName); 
      free(strmen);
      if (position == 0){ if(pscrPtr) *pscrPtr=NULL;  return 0;}
      if(pscrPtr)   *pscrPtr=pscr; else  put_text(&pscr);
      if(sqtrS_on)
      { shift=1; for(i=0;i<2;i++) if(polar[i]) shift++;   
        if(position<=shift) switch(position)
        { case 1: if(polar[0]) {*varPos=&Helicity[0]; return 1;}
          case 2: if(polar[1]) {*varPos=&Helicity[1]; return 1;}
          case 3:               *varPos=&Pcm22;     return 1;
        }
           
        position -= shift;       
      } 
      if(!vars_on)  position+=nModelVars;
      if(rd_on) position--;
      *varPos= varValues+position-1;
      return 1;
    } else return 0;
} 

int change_parameter(int x,int y, int for22)
{ 
  double val;
  char name[40];
  REAL * vPos;
  int i,err,mPos=1;
  int returnCode=0; 
  void * pscr=NULL;
  REAL *va_mem=(REAL*)malloc(sizeof(REAL)*(nModelVars+1));
  double Pcm_mem=Pcm22;

  for(i=0;i<nModelVars;i++) va_mem[i]=varValues[i];
  
  for(err=1; err;)
  {  
    for(;selectParam(x,y,"Change parameter",&pscr,for22,1,1,0,&vPos,name,&mPos);)
    {  
      if(vPos==varValues-1) 
      { FILE *f;
        static char fName[100];
        int err;
        struct stat buf;
        
        if(!findCalcHEPfile(fName)) continue; 
        if(stat(fName,&buf) ||   !(S_ISREG(buf.st_mode)) )  { messanykey(10,17, "Not a regular file"); continue; }
        f=fopen(fName,"r");
        if(f==NULL) { messanykey(10,17, "Can not open file"); continue; }
        for(;;)
        { char name[20];
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
          } else  { varValues[i]=val; returnCode=1;}
          
          fscanf(f,"%*[^\n]");
        }  
        fclose(f);        
      }
      else
      {
        strcat(name," = ");
        val=*vPos;
        if(correctDouble(x,y+4,name,&val,1)) { *vPos=val; returnCode=1;}
      }    
    }
    if(returnCode)
    { 
      void * pscr;
      get_text(x,y+1,x+25, y+3,&pscr);
      scrcolor(FGmain,BGmain);
      goto_xy(x+5,y+1); print("Be patient:");
      goto_xy(x,y+3); print("Calculation of constraints");
      escpressed();
      err=checkParam();
      if(Warnings) messanykey(5,10,Warnings);
      put_text(&pscr);
    } else  err=0;
    if(err && mess_y_n(10,10, "Restore previous parameter set?"))
    {   
       for(i=0;i<=nModelVars;i++) varValues[i]=va_mem[i];
       Pcm22=Pcm_mem;
    }
  }
  free(va_mem);
  return returnCode;
}

void show_depend(int x, int y)
{ void *pscr1=NULL;
  int i,mPos=1; 

  if(!nModelFunc) { messanykey(5,10,"There are no public constraints in this model"); return;}
  REAL*allfunc=(REAL*) malloc(sizeof(REAL)*nModelFunc);
  for(i=0;i<nModelFunc;i++) allfunc[i]=varValues[1+nModelVars+i];

  for(;;) 
  { void *pscr2=NULL;
    char name1[20];
    REAL * fPos;
    
    if(!selectParam(x,y+1,"Display dependence",&pscr1,0,0,0,1,&fPos,name1,&mPos)) return;
    
    for(;;)
    { char name2[20];
      double val2;
      void *pscr3=NULL;
      double xMin, xMax;
      int  nPoints=100;
      int k3;
      REAL *vPos;
      int mPos_=1;
      
      if(!selectParam(x,y+5,"on parameter",&pscr2,0,0,1,0,&vPos,name2,&mPos_))break;
   
      val2=*vPos; 
      xMin=val2 - fabs(val2)/10;
      xMax=val2 + fabs(val2)/10; 
      
      for(;;)
      {   
         char strmen[]="\030 "
            " x-Min = XXX            "
            " x-Max = YYY            "
            " Npoints = NNN          "
            " Display                ";
     
         improveStr(strmen,"XXX","%G",xMin);
         improveStr(strmen,"YYY","%G",xMax);
         improveStr(strmen,"NNN","%d",nPoints);

         
         menu1(x,y+9,"Plot",strmen,"",&pscr3,&k3);
         if(!k3) break;
         switch(k3)
         {  case 1: correctDouble(x,y+12,"xMin = ",&xMin,1); break;
            case 2: correctDouble(x,y+12,"xMax = ",&xMax,1); break;
            case 3: correctInt(x,y+12,"nPoints = ",&nPoints,1); break;
            case 4:
            if( xMax>xMin && nPoints>=3 && nPoints<=150)
            {  double dx=(xMax-xMin)/(nPoints-1);
               double f[150];
               int i, NaN=0,Esc=0;
         
               informline(0,nPoints);               
               for(i=0;i<nPoints;i++)
               {  double x=xMin+i*dx;
                  *vPos=x;
                  NaN=checkParam();
                  if(NaN) 
                  {  char mess[100];
                     sprintf(mess,"Can not evaluate constraints for %s=%G",name2, x);
                     messanykey(16,5,mess);        
                     break;
                  }
                  f[i]=*fPos;
                  Esc=informline(i,nPoints);
                  if(Esc) break;  
               }
                  
               *vPos=val2; 
               for(i=0;i<nModelFunc;i++) varValues[1+nModelVars+i]=allfunc[i];

               if(!(NaN||Esc)) plot_1(xMin,xMax,nPoints,f,NULL,"Plot",
                               name2,name1);
                               
            } else messanykey(16,5," Correct input is \n"
                                   "  xMin<xMax,\n"
                                   " 3<=nPoints<=150");
            break;
         }
       }
     }
  }
  free(allfunc);
}
