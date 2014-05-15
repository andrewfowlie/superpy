/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include "help.h"
#include "crt_util.h"
#include <ctype.h>
#include "pictures.h" 

void  pictures(int  n_diag_tot, int sizeX,int sizeY, 
void (*pict)(int,int,int),char* commandStr, int (*command)(int, char),
char * help)                    
{
/* to  sreen restore */ 
 void *   prtscr;

/* up command line */
#define DNMAX 8 
int  dnKey[DNMAX] = {KB_F1,KB_F2,KB_PAGEU,KB_PAGED,KB_HOME,KB_END,'#',KB_ESC};
int  dnLen[DNMAX], imax_dn;
char dnStr[] = "F1-Help,F2-Man,PgUp,PgDn,Home,End,# ,Esc";

/* down command line */
#define UPMAX 10
int  upKey[DNMAX], upLen[DNMAX],imax_up;
char * upStr=commandStr;

 int n=1;   /* current position */
 int nFirst=1;  /* first diagram on screen */
 
 int cHeight=tg_textheight("H"); /* Letter's height */

 int maxx,maxy; /* screen size */
 int  ndi_x,ndi_y; /* number of rows and columns */
 int lx,ly; /* cell size */
 int x0=1,y0=cHeight+2; /* left-top corner of the first cell */
 int dx,dy; /* chink */
 int  key; /* for inkey() */

 int redraw;
 int  i,j,i0,j0,l;
 char st[100];
 char *p;

/*----------------------------------------*/
/* save screen */
   get_text(1,1,maxCol(),maxRow(),&prtscr);
   
/* down command line preparing */
   imax_dn=0;
   p=dnStr-1;
      
   do
   { p++;
     if (*p==','||*p==';'|| *p==0)
     { char cc=*p;
       *p=0;
       dnLen[imax_dn++]=tg_textwidth(dnStr);
       *p=cc;
     }         
   } while(*p && imax_dn<DNMAX);

/* upper command line preparing */
   imax_up=0;
   p=upStr-1;
      
   do
   { p++;
     if (isupper(*p)) upKey[imax_up]=*p ;
     if (*p==','||*p==';'|| *p=='\0')
     { char cc=*p;
       *p=0;
       upLen[imax_up++]=tg_textwidth(upStr);
       *p=cc;
     }         
   } while(*p && imax_up<UPMAX);
   
   x0=1;
   y0=cHeight+2;

redrawScreen:
   redraw=1;
   dy=4;
   dx=4; 

   maxx=tg_getmaxx();
   maxy=tg_getmaxy();
   
/* ndi_x,ndi_y - number of cell */
   dx=3;  ndi_x=(maxx-2*x0-dx)/(dx+ sizeX );
   dy=3;  ndi_y=(maxy-2*y0 -dy)/(dy+sizeY);

/* lx,ly - cell size */   
   lx=(maxx-2*x0)/ndi_x;
   ly=(maxy-2*y0)/ndi_y;

   dx=(lx - sizeX)/2;
   dy=(ly - sizeY)/2;
  
/* begin to write */
   clr_scr(Black,White);
   tg_setviewport(0,0,maxx,maxy);
   
/* up command string */
   scrcolor(White,Black);
   tg_bar(0,0,maxx,cHeight);
   tg_settextjustify(LeftText,TopText);
   tg_outtextxy(0,0,upStr);

/*down command string */   
   tg_bar(0,maxy-cHeight ,maxx,maxy);
   tg_settextjustify(LeftText,BottomText);
   tg_outtextxy(0,maxy,dnStr);
   
/* cells drawing */
   scrcolor(Blue,White);
   tg_setlinestyle(SolidLn,NormWidth);
   for(i=0;i<=ndi_x;i++) { tg_line(x0+i*lx,y0,x0+i*lx,y0+ndi_y*ly);}
   for(i=0;i<=ndi_y;i++) { tg_line(x0 ,y0+i*ly, x0+ndi_x*lx,y0+i*ly);}

   st[0]=0;
   
   for(;;)
   {
/* new n correction */
      if(n<1) n=1;
      if(n>n_diag_tot) n=n_diag_tot;
      
/* find new page beginning */      
      l= 1+(ndi_x*ndi_y)*((n-1)/(ndi_x*ndi_y));
      if (l!=nFirst) {redraw=1; nFirst=l;}

/* writing the number of the current diagram */
      tg_settextjustify(RightText,TopText);
      tg_setviewport(0,0,maxx,maxy);
/* remove old number */      
      scrcolor(Black,White); tg_outtextxy(maxx,0,st);
/* new number */
      sprintf(st,"%d/%d ",n,n_diag_tot);
      scrcolor(White,Black); tg_outtextxy(maxx,0,st);       
/*diagrams drawing*/
      l=nFirst;
      for(j=0;j<ndi_y; j++) for(i=0;i<ndi_x;i++)
      { 
        if( redraw)
        {
           scrcolor(Black,White);
           tg_bar(x0+i*lx+1,      y0+j*ly+1,
                  x0+i*lx+lx-1,   y0+j*ly+ly -1);
           if(l<=n_diag_tot && redraw) (*pict)(l,x0+i*lx+dx, y0+j*ly+dy);
          
        }
        scrcolor(Black,White);
        tg_setlinestyle(SolidLn,ThickWidth);
        if(l==n) 
        {  tg_rectangle(x0+i*lx+2,    y0+j*ly+2,
                        x0+i*lx+lx-2, y0+j*ly+ly-2); 
            i0=i;j0=j;
        }
        l++;
      }  
      redraw=0; 
      
      key = inkey();

/* mouse event*/
      if (key == KB_MOUSE && mouse_info.but1 == 2)
      {  key=0;
/* up line */       
         if (mouse_info.y <= cHeight) 
         { 
            i = 0;
            while (i < imax_up && upLen[i] < mouse_info.x) i++;
            if(i!=imax_up)  key =upKey[i];               
         }
/* diagrams area */         
         else  if (mouse_info.y <= maxy-y0)
         {  if (mouse_info.x > x0 && mouse_info.x < maxx-x0)
            {  
               i = (mouse_info.x - x0 )/lx ;
               j = (mouse_info.y -  y0)/ly ;
               n=nFirst+i+j*ndi_x;               
            }
         }
/* down line */         
         else
         {  i = 0;
            while (i < imax_dn && dnLen[i] < mouse_info.x) i++;
            if(i!=imax_dn) key =dnKey[i];               
         }
      }
/* mane switch */     
      switch (key)
      { 
         case KB_F1:   show_help(help);     break;
         case KB_F2:   show_help("diagram");break;
         case KB_ESC:  clr_scr(FGmain,BGmain); put_text(&prtscr);
                                            return;
         case KB_HOME: n=1;                 break;
         case KB_END:  n=n_diag_tot;        break; 
         case KB_PAGED:n += (ndi_x*ndi_y);  break;
         case KB_PAGEU:n -= (ndi_x*ndi_y);  break;
         case KB_LEFT: n--;                 break;
         case KB_RIGHT:n++;                 break; 
         case KB_UP:   n -= ndi_x;          break;
         case KB_DOWN: n += ndi_x;          break;
         case '#': 
               correctInt( 10,10," New position = ",&n,1);         
                                            break;
         default:
         if(isalpha(key))
         { key=toupper(key);
           if(command)
           for(i=0;i<imax_up;i++) if (upKey[i]==key) {(*command)(n, key); break;}
           redraw=1;
         }                                    
      } /* end of switch */

      if(maxx != tg_getmaxx() || maxy != tg_getmaxy() ) goto redrawScreen;      

/* frame box clearing*/
      scrcolor(White,Black); 
      tg_setlinestyle(SolidLn,ThickWidth);
      i=i0;j=j0;      
      tg_rectangle(x0+i*lx+2,    y0+j*ly+2,
                   x0+i*lx+lx-2, y0+j*ly+ly-2);       
   }  /* end of while */
}
