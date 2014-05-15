/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include <unistd.h>
#include"chep_crt.h"
#include <math.h>
#include "syst2.h"
#include "physics.h"
#include "s_files.h"
#include "r_code.h"
#include "drawdiag.h"
#include "draw_ampl.h"
#include "process.h"
#include "showgrph.h"

/*======================================================================*/
  /*
  nn   - first FilePos of subproces;
  nm-1 - last  FilePos of subproces;
  */
 
static   char stc[100];
static   long   nn, nm;
static   int pictureX;
static   int pictureY;  
static   int xn,ynu;

static FILE * diagrFile;
static int diagSize;
static char * mark;
static int squared;

static void  wrttext(int n,int x,int y)
{  
   int ypos=y + ynu -2;
   tg_settextjustify(LeftText,BottomText);
      switch(mark[n-1])
      {  
         case -1:
            scrcolor(Red,White);
            tg_outtextxy(x + 48,ypos,"DEL");
            break;         
         case -2:
            scrcolor(Red,White);
            tg_outtextxy(x+48,ypos,"Out of memory");
         break;
         case  1:
            scrcolor(LightBlue,White);
            tg_outtextxy(x+48,ypos,"CALC");
         break;
         case  2:
            scrcolor(LightBlue,White);
            tg_outtextxy(x+48,ypos,"ZERO ");
      }
      scrcolor(Black,White);
}

static void  setPictureScale_(int squared, int*pictureX,int*pictureY)
{  
  if(squared && nin+nout< 7)setPictureScale(squared, pictureX,pictureY);
  else {*pictureX=amplitudeFrameX(nin,nout); *pictureY=amplitudeFrameY(nin,nout);
  if(squared) *pictureX =2*(*pictureX);
  }
}

static void picture_(int squared, char* buff, int x, int y)
{
    if(squared &&  nin+nout< 7) picture(squared,buff,x,y);
    else if(!squared) 
   {  vampl V;
      mkverts((particleNumType *)buff,&V);
      drawAmpitudeDiagram(&V,0,x,y);
   }  
   else 
/*
ypedef struct csdiagram
   {
      decayDiagram  dgrm1,dgrm2;
      int        lnk[MAXINOUT];
      int           mult;
      unsigned          del;
      char          status;   // -2-outOfMemory,-1-deleted, 0-Rest, 1-calculated,2-Zero  
      int nsub;
   }  csdiagram;
*/   
   {    csdiagram*sqD=(csdiagram*)buff;
        vampl V;
        int pictureX;
         
        mkverts(sqD->dgrm1,&V);
        drawAmpitudeDiagram(&V,0,x,y);
        mkverts(sqD->dgrm2,&V);
        pictureX=amplitudeFrameX(nin,nout);
        drawAmpitudeDiagram(&V,1,x+2*pictureX,y); 
        return;
   }  
}
 

static void new_picture(int n, int x, int y)
{  char buff[MAX(sizeof(csdiagram),sizeof(adiagram))];
   wrttext(n,x,y);   
   n+=(nn-1);
   fseek(diagrFile,n*diagSize,SEEK_SET);
   fread(buff,diagSize,1,diagrFile);
   picture_(squared,buff,x,y);
}  

/* tex tex tex tex tex tex tex tex tex tex */

static void texinvoke(char* pname)
{  
   char buff[MAX(sizeof(csdiagram),sizeof(adiagram))];   
   char dd[30],txt[200];
   int  n_,dn_;
   int flag;

   char    f_name[STRSIZ];

   if (squared)  strcpy(dd,"csd_"); else  strcpy(dd,"fd_");
   nextFileName(f_name,dd,".tex"); 
   sprintf(txt,"diagrams for process %s",pname);
   texStart(f_name,txt,"scriptsize" );    

   texxscale=1;
   texyscale=1;

   setPictureScale_(squared, &pictureX,&pictureY);

   texxshift=0;
   texymax1=pictureY;

   pictureX*=texxscale;
   pictureY*=texyscale;

   dn_=0;

   for (n_ = nn; n_ < nm; n_++)
   {
     flag= ( mark[n_-nn] !=-1);
     if (flag)
     {  fseek(diagrFile,n_*diagSize,SEEK_SET) ;
        fread(buff,diagSize,1,diagrFile);       
        dn_++;
        f_printf(out_tex,"{} \\qquad\\allowbreak\n");
        f_printf(out_tex,"%%  diagram # %d\n",(n_ -(int)nn)+1);
        f_printf(out_tex,"\\begin{picture}(%d,%d)(0,0)\n",pictureX,pictureY);             
        picture_(squared,buff,0,0); 
        f_printf(out_tex,"\\end{picture} \\ \n");
     }
   }

   texFinish();
   setPictureScale_(squared, &pictureX,&pictureY);
/*pictureY/=2; */
   sprintf(txt,"LaTeX output is written in file\n %s",f_name);
   messanykey(20,14,txt);

}


/* end tex */

static int comList( int n, char key)
{
 int l;

switch( key)
{

  case 'L':  /* LaTex output for all undel. diagrs. */
    texinvoke(stc);   
    return 0 /*nothing to redraw */;

  case 'D':        /* Del -   Delete all subprocess  */
    for(l=0;l<nm-nn;l++)  {if(mark[l]==0|| mark[l]==-2) mark[l]=-1;}    
    return 2 /* redrawScreen */;
  
  case 'R': /*  Ins -  Restore all deleted  */
    for(l=0;l<nm-nn;l++) { if(mark[l]==-1) mark[l]=0;}   
    return 2 /* redrawScreen */;

  case 'O':   /*  Space  Bar   del/ins proceess  */
    if(mark[n-1]==0||mark[n-1]==-2 ) mark[n-1]=-1;
    else if(mark[n-1]==-1) mark[n-1]=0; 
    return 1;
    
  case 'G':
    {        
       FILE * txt;
       makeghostdiagr(nn+n-1,"view.tmp"); 
       txt=fopen("view.tmp","r");
       showtext(1, 1 , 80,1,"Ghosts",txt);
       fclose(txt);
       unlink("view.tmp");
       return 0;    
    }    
  }
}


void  showgraphs(char upravl)
{ 
   int   ndel, nrest, ncalc;
   int onlyview;
   int i;
   int delMarkPos;
   int ntot;
   char help[80]; 
   char comLine[80];
   
   onlyview = (upravl < 0);
   if (onlyview) upravl = -upravl;

   squared=(upravl>1);
   
   setPictureScale_(squared, &xn,&ynu);

   {  adiagram  varc;
      csdiagram ars;
     
      if (squared)
      {
        menuq=fopen(MENUQ_NAME,"r+b");
        diagrFile=fopen(DIAGRQ_NAME,"r+b");
        diagSize=sizeof(csdiagram);
        delMarkPos = (char*)&(ars.status) - (char*)&ars;
      }
      else
      { 
        menup=fopen(MENUP_NAME,"r+b");    
        diagrFile=fopen(DIAGRP_NAME,"r+b");
        diagSize=sizeof(adiagram);
        delMarkPos=(char*)&(varc.delMark)-(char*)&varc;
      }
   } 
          
/*  nn   -  Position of first diagrams in subproces  */

   rd_menu(upravl,nsub,stc,&ndel,&ncalc,&nrest,&nn);
   ntot=ndel+ ncalc+nrest; 
   nm=nn+ntot;



   if(squared)  
   {  strcpy(help,"s_diag_s");
      strcpy(comLine,"Delete,On/off,Restore,Latex,Ghosts");
   }
   else if (onlyview) {strcpy(help,"s_diag_v");strcpy(comLine,"Latex");}
   else {strcpy(help,"s_diag_e");strcpy(comLine,"Delete,On/off,Restore,Latex");}
      
   mark=malloc(ntot); 
   for(i=0;i<ntot;i++) mark[i]=0;
 
   
   for (i=nn;i<nm;i++)
   {  
     fseek(diagrFile,i*diagSize+delMarkPos,SEEK_SET);
     fread(mark+i-nn,1,1,diagrFile);
   }
   
   pictures(ntot,xn,ynu, new_picture,comLine,comList,help);

   ndel=0;
   ncalc=0;
   for (i=nn;i<nm;i++)
   {
     fseek(diagrFile,i*diagSize+delMarkPos,SEEK_SET);
     fwrite(mark+i-nn,1,1,diagrFile);
     if(mark[i-nn]>=1) ncalc++;
     else if(mark[i-nn]==-1) ndel++;
   }

   nrest=ntot-ndel-ncalc;
   wrt_menu(upravl,nsub,stc,ndel,ncalc,nrest,nn);

   if (squared) fclose(menuq); else  fclose(menup);   
   free(mark);
   fclose(diagrFile);
}
