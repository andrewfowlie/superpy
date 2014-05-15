/*
 Copyright (C) 1997, Victor Edneral, e-mail: edneral@theory.npi.msu.su
 Corrected by A.Pukhov 
*/

#include "syst.h"
#include "crt.h"
#include <math.h>
#include"crt0.h"
int texflag=0;
int blind = 1;
char    boxFrame[9] = {1,2,3,4,5,6,7,8,9};
int nColors=15;
int fColor=White;
int bColor=Black;
char * inkeyString;

 double  texxscale=0.5;
 double  texyscale=0.5;
 int     texxshift,texymax1;
 FILE *  out_tex=NULL;


static int  xc = 0, yc = 0;
static double texCharScale;
static int bgi_x=0,bgi_y=0, bgi_maxx=639,bgi_maxy=479,CdX=8, CdY=19;

static int bgi_hj=LeftText,bgi_vj=TopText;

static struct tg_viewporttype bgi_vport={0,0,639,479};
static struct tg_linesettingstype  bgi_line={SolidLn,NormWidth};

enum graphObj { line,bar,string,t_xt};


typedef struct graphCommand
{ struct graphCommand * next;int type;
  int x1, y1, x2, y2, color,size,aux;
  char * text;} graphCommand;

typedef  graphCommand * graphComPtr;

static graphComPtr graphList=NULL;

static void print0(int xc,int yc, int color,int bkcolor, char * s);

static void texline(int x1,int y1,int x2,int y2, int arrowflag) /* by A.Pukhov */
{ struct tg_linesettingstype ls;
  double x,y,fx,fy;
  double dashsiz=3, dotsiz=1;

   tg_getlinesettings(&ls);
   x = (x1+texxshift)*texxscale;
   y = (texymax1-y1)*texyscale;
   fx = (x2+texxshift)*texxscale;
   fy = (texymax1-y2)*texyscale;

   switch(ls.linestyle)
   {  case SolidLn:
      switch(arrowflag)
      {
      case 1  :
         f_printf(out_tex,
            "\\ArrowLine(%.1f,%.1f)(%.1f,%.1f) \n" ,x,y,fx,fy);
      break;
      case -1 :
         f_printf(out_tex,
            "\\ArrowLine(%.1f,%.1f)(%.1f,%.1f) \n" ,fx,fy,x,y);
      break;
      case 0  :
         f_printf(out_tex,
            "\\Line(%.1f,%.1f)(%.1f,%.1f) \n" ,x,y,fx,fy);
      }
      break;
      case DottedLn:
      switch(arrowflag)
      {
      case 1  :
         f_printf(out_tex,
          "\\DashArrowLine(%.1f,%.1f)(%.1f,%.1f){%.1f} \n" ,x,y,fx,fy,dotsiz);
      break;
      case -1 :
         f_printf(out_tex,
          "\\DashArrowLine(%.1f,%.1f)(%.1f,%.1f){%.1f} \n" ,fx,fy,x,y,dotsiz);
      break;
      case 0  :
         f_printf(out_tex,
            "\\DashLine(%.1f,%.1f)(%.1f,%.1f){%.1f}\n" ,x,y,fx,fy,dotsiz);
      }

      break;
      case DashedLn:
      switch(arrowflag)
      {
      case 1  :
         f_printf(out_tex,
          "\\DashArrowLine(%.1f,%.1f)(%.1f,%.1f){%.1f} \n" ,x,y,fx,fy,dashsiz);
      break;
      case -1 :
          f_printf(out_tex,
          "\\DashArrowLine(%.1f,%.1f)(%.1f,%.1f){%.1f} \n" ,fx,fy,x,y,dashsiz);
      break;
      case 0  :
         f_printf(out_tex,
          "\\DashLine(%.1f,%.1f)(%.1f,%.1f){%.1f} \n" ,x,y,fx,fy,dashsiz);
      }
   }
}

static void texBar(int x1,int y1,int x2,int y2) /* by A.Pukhov */
{ struct tg_linesettingstype ls;
  double x,y,fx,fy;

   tg_getlinesettings(&ls);
   x = (x1+texxshift)*texxscale;
   y = (texymax1-y1)*texyscale;
   fx = (x2+texxshift)*texxscale;
   fy = (texymax1-y2)*texyscale;
   f_printf(out_tex,"\\BBox(%.1f,%.1f)(%.1f,%.1f) \n" ,x,y,fx,fy);
}



static void texouttext(int x,int y,char* txt)   /* by A.Pukhov */
{  struct tg_textsettingstype ts;
   f_printf(out_tex,"\\Text(%.1f,%.1f)[",(x+texxshift)*texxscale,
        (texymax1-y)*texyscale);

   tg_gettextsettings(&ts);

   switch(ts.horiz)
   {  case LeftText:  f_printf(out_tex,"l"); break;
      case RightText: f_printf(out_tex,"r"); break;
   }

   switch(ts.vert)
   {  case BottomText: f_printf(out_tex,"b"); break;
      case TopText:  f_printf(out_tex,"t");   break;
   }
   f_printf(out_tex,"]{$%s$}\n",txt);
}


int setTexCharSize(char * font)   /* by A.Pukhov */
{
        if(strcmp(font,"tiny") ==0 ) texCharScale=0.672;
   else if(strcmp(font,"scriptsize") ==0 ) texCharScale=0.782;
   else if(strcmp(font,"footnotesize") ==0 ) texCharScale=0.837;
   else if(strcmp(font,"small") ==0 ) texCharScale=0.927;
   else if(strcmp(font,"normalsize") ==0 ) texCharScale=1.0;
   else if(strcmp(font,"large") ==0 ) texCharScale=1.18;
   else if(strcmp(font,"Large") ==0 ) texCharScale=1.44;
   else if(strcmp(font,"LARGE") ==0 ) texCharScale=1.74;
   else if(strcmp(font,"huge") ==0 ) texCharScale=2.07;
   else if(strcmp(font,"Huge") ==0 ) texCharScale=2.07;
   else { texCharScale=0.8; return 0;}
   return 1;  

}





int   tg_getmaxx(void)    {	return bgi_maxx;   }
int   tg_getmaxy(void)    {	return bgi_maxy;   }
int   tg_getx(void)       {	return bgi_x;      }
int   tg_gety(void)       {	return bgi_y;      }
int   maxCol(void)        {     return MIN( (bgi_maxx+1)/CdX,STRSIZ-1);}
int   maxRow(void)        {     return (bgi_maxy+1)/CdY;}   

void	tg_gettextsettings(struct tg_textsettingstype * txi)
{
   txi->horiz=bgi_hj;
   txi->vert=bgi_vj;
}

void  tg_getviewsettings(struct tg_viewporttype *   vp)
{
   vp->left=bgi_vport.left;
   vp->right=bgi_vport.right;
   vp->top=bgi_vport.top;
   vp->bottom=bgi_vport.bottom;
}

void	 tg_setviewport( int left,  int top, int right, int bottom)
{
		bgi_vport.left=left;
		bgi_vport.top=top;
		bgi_vport.right=right;
		bgi_vport.bottom=bottom;
		bgi_x=0;
		bgi_y=0;
}

void	tg_setlinestyle(int linestyle, int thickness)
{
	bgi_line.linestyle=linestyle;
	bgi_line.thickness=thickness;
}

static void   addToGraphList(int type,int x1,int y1,int x2,int y2,
int color,int aux,char*text)             /* by A.Pukhov */
{ graphComPtr p=graphList ,pred=NULL;

  while(p)
  {  int del=0;

     switch(type)
     { case bar:
	del=( MIN(p->x1,p->x2) >=MIN(x1,x2) && MAX(p->x1,p->x2) <=MAX(x1,x2)&&
	   MIN(p->y1,p->y2) >=MIN(y1,y2) && MAX(p->y1,p->y2) <=MAX(y1,y2)  );
	    break;
       case t_xt:
        if(p->type==t_xt && p->y1==y1)
        {
           del=( MIN(p->x1,p->x2) >=MIN(x1,x2) && MAX(p->x1,p->x2) <=MAX(x1,x2)&&
	   MIN(p->y1,p->y2) >=MIN(y1,y2) && MAX(p->y1,p->y2) <=MAX(y1,y2)  );
	}
            break;
       case line: del=p->type==line && p->x1==x1 && p->y1==y1
				     && p->x2==x2 && p->y2==y2;
	    break;
       case string:  del=p->type==string && p->x1==x1 && p->y1==y1 &&
		   strcmp(p->text,text)==0;
	    break;
     }
     if(del)
     {
	 if(pred==NULL)
	 {  graphList=p->next;
	    if (p->text) free(p->text);
            free(p);
            p=graphList;
         }else
	 {  pred->next=p->next;
	    if (p->text) free(p->text);
            free(p);
            p=pred->next;
         }
     }else
     {  pred=p;
	p=p->next;
     }
  }
  
  p=(graphComPtr)m_alloc(sizeof(graphCommand));
  p->type=type;
  p->x1=x1;
  p->y1=y1;
  p->x2=x2;
  p->y2=y2;
  p->color=color;
  p->aux=aux;
  if(text==NULL) p->text=NULL; else
  { p->text=(char*)m_alloc(strlen(text)+1);
    strcpy(p->text,text);
  }
  p->next=NULL;
  if(pred==NULL) graphList=p; else pred->next=p;
}

static void repeatGrCom(int x1,int y1,int x2,int y2,int all) /* by A.Pukhov */
{  graphComPtr p = graphList;
   graphComPtr pred = NULL;
   graphComPtr p1,p2;
   int dn,xx1,xx2,yy1,yy2;
   int count =0;

   if(blind==1) return;

/*   printf(" repeatGrCom   "); */
   sg_screenSize(&bgi_maxx,&bgi_maxy);
   if( graphList !=NULL && graphList->type == bar)
   { 
      graphList->x2=bgi_maxx;
      graphList->y2=bgi_maxy;
   }     
   sg_setClip(x1,y1,x2,y2);
   while (p!=NULL)
   { count ++;
     xx1=MIN(p->x1,p->x2);
     xx2=MAX(p->x1,p->x2);
     yy1=MIN(p->y1,p->y2);
     yy2=MAX(p->y1,p->y2);
     if(yy2 >= y1 && yy1 <= y2 && xx1<=x2 && xx2>=x1)
     switch(p->type)
     {
      case bar:  sg_drawBox(p->x1,p->y1,p->x2,p->y2,p->color);
	 break;
      case line: sg_drawLine(p->x1,p->y1,p->x2,p->y2,p->color,
	p->aux&15,p->aux/16);
	 break;
      case string:  sg_outText(p->x1,p->y1,p->color,p->text);
         break;
      case t_xt:
        if (all) print0(p->x1,p->y1,p->color,p->aux, p->text);else
        {
           if(p->x1<x1)
           { p1=(graphComPtr)m_alloc(sizeof(graphCommand));
             (*p1)=(*p);
             p1->x2=x1-1;
             dn=(x1 - p->x1)/CdX;
             p1->text=(char*)m_alloc(dn+1);
             memcpy(p1->text,p->text,dn+1);
             p1->text[dn]=0;
           } else p1=NULL;

           if(p->x2>x2)
           { p2=(graphComPtr)m_alloc(sizeof(graphCommand));
             (*p2)=(*p);
             p2->x1=x2+1;
             dn=(p->x2-x2)/CdX;
             p2->text=(char*)m_alloc(dn+1);
             strcpy(p2->text,p->text+ strlen(p->text)-dn );
           } else p2=NULL;
           free(p->text);

           if(p1!=NULL)
           { (*p)=(*p1);
             free(p1);
             if(p2!=NULL) p->next=p2;
           } else
           { if(p2!=NULL)
             { (*p)=(*p2);
               free(p2);
             } else
             {
               if(pred==NULL)
               {  graphList=p->next;
                  free(p);
                  p=graphList;
               }else
               {  pred->next=p->next;
                  free(p);
                  p=pred->next;
               }
               goto contin;
             }
           }
        }
     }
     pred=p;
     p=p->next;
contin:;
   }
   sg_setClip(0,0,bgi_maxx,bgi_maxy);
/*   printf( "   ( %d commands )\n",count); */
}



static void crt_expose(int x1,int y1,int w,int h)  /* by A.Pukhov */
{   
    repeatGrCom(x1,y1,x1+w,y1+h,1);
    refresh_scr();     
}

void tg_bar(int left,  int top, int right, int bottom)
{ int x1,x2,y1,y2;
  if(texflag) 
  {
    texBar(left,top,right,bottom);
  }
  else if(blind!=1)
  {
    x1=left +bgi_vport.left;
    y1=top+bgi_vport.top;
    x2=right+bgi_vport.left;
    y2=bottom+bgi_vport.top;
    sg_drawBox(x1,y1,x2,y2,bColor);
    addToGraphList(bar,x1,y1,x2,y2,bColor,0,NULL);
  }
}

void  tg_clearviewport(void)
{
   tg_bar(0,0,bgi_vport.right-bgi_vport.left,bgi_vport.bottom-bgi_vport.top);
   bgi_x=0;
   bgi_y=0;
}

void  tg_line(int x1, int y1, int x2, int y2)
{
  if(texflag) texline(x1,y1,x2,y2,0); 
  else
  {
     x1+=bgi_vport.left;
     y1+=bgi_vport.top;
     x2+=bgi_vport.left;
     y2+=bgi_vport.top;
     if(blind!=1)
     { sg_drawLine(x1,y1,x2,y2,fColor,bgi_line.thickness,bgi_line.linestyle);
       addToGraphList(line,x1,y1,x2,y2,fColor,
		 bgi_line.thickness+16*bgi_line.linestyle,NULL);
     }
  }
}


void  tg_arrowline(int x1, int y1, int x2, int y2)
{
   if(texflag) texline(x1,y1,x2,y2,1);
   else
   {  double dx,dy,dl,xc,yc;
      int bgi_line_linestyle;
      tg_line(x1,y1,x2,y2);
      dx=x2-x1;
      dy=y2-y1;
      dl=sqrt(dx*dx+dy*dy);
      dx=2*dx/dl;
      dy=2*dy/dl;
      xc=(x1+x2+2*dx)/2;
      yc=(y1+y2+2*dy)/2;
      bgi_line_linestyle=bgi_line.linestyle;
      bgi_line.linestyle=SolidLn;
      tg_line((int)xc,(int)yc,(int)(xc-2*dx-2*dy),(int)(yc-2*dy+2*dx));
      tg_line((int)xc,(int)yc,(int)(xc-2*dx+2*dy),(int)(yc-2*dy-2*dx));
      bgi_line.linestyle=bgi_line_linestyle;
   }
}

void	tg_linerel( int dx, int dy)
{
	tg_line(bgi_x,  bgi_y,
		  bgi_x+dx,	bgi_y+dy);
	bgi_x+=dx;
	bgi_y+=dy;
}

void	tg_lineto( int x, int y)
{
	tg_line( bgi_x,	bgi_y,
				 x,       y);
	bgi_x=x;
	bgi_y=y;
}

void	 tg_moverel( int dx, int dy)
{
	bgi_x+=dx;
	bgi_y+=dy;
}

void tg_moveto(int x, int y)
{
	bgi_x=x;
	bgi_y=y;
}

void tg_settextjustify( int horiz, int vert)
{
	bgi_hj=horiz;
	bgi_vj=vert;
}

int tg_textwidth(char* textstring)
{
  int dx,dy;
  if(texflag)
  { int sum=0; int i=0;
    while(textstring[i]!=0)
    {      if (textstring[i]<'A') sum+=5;
      else if (textstring[i]<'a') sum+=8;
      else  sum+=5;
      i++;
    }
    return (int)(sum*texCharScale/texxscale);
  }
  if(blind!=1) sg_textSize(textstring , &dx, &dy);
  else       dx=strlen(textstring)*10;
  return  dx;
}

int tg_textheight(char* textstring)
{
  int dx,dy;
  if(texflag) return (int)(10*texCharScale/texyscale);
  if(blind!=1) sg_textSize(textstring , &dx, &dy); else dy=8;
  return  dy;
}

void  tg_outtextxy( int x, int y, char* textstring)
{ int dx,dy;

  if(texflag) texouttext(x,y,textstring); else
  {  x+=bgi_vport.left;
     y+=bgi_vport.top;
     if(blind!=1) 
     {  sg_textSize(textstring,&dx,&dy);
        x-=(dx*bgi_hj)/2 ;
        y+=(dy* bgi_vj)/2 ;
        sg_outText(x,y,fColor,textstring);
	addToGraphList(string,x,y,x+dx-1,y-dy+1,fColor,0,textstring);
     }
  }
}


void  tg_outtext(char* textstring)
{
	tg_outtextxy(bgi_x,bgi_y,textstring);
	if(bgi_hj==LeftText)
			bgi_x+=tg_textwidth(textstring);
}

void	 tg_rectangle( int left, int top, int right, int bottom )
{
	int xx,yy;
	xx=bgi_x;
	yy=bgi_y;

	tg_moveto(left, top);
	tg_lineto(right,top);
	tg_lineto(right,bottom);
	tg_lineto(left, bottom);
	tg_lineto(left,top);
	tg_moveto(xx,yy);
}

void tg_getlinesettings(struct tg_linesettingstype *  sls)
{
   sls->linestyle=bgi_line.linestyle;
   sls->thickness=bgi_line.thickness;
}

void refresh_scr(void) { if(blind!=1) crt0_keypressed(); }

int  start1(char * window_name,char * icon_file,char * ini_file, 
   void(*xw_errorExit)(void))   /* by A.Pukhov */
{  int err=0;
   if(blind!=1)
   { err=crt0_start(window_name,icon_file,ini_file,&nColors,xw_errorExit);
      xw_expose=(*crt_expose);
      crt0_charSize(&CdX,&CdY);
      clr_scr(FGmain,BGmain);
      if(err==1)
      {  printf("You have launched CalcHEP compiled without X11\n"
                "Current version can work only in the 'blind' mode\n"
                "To use CalcHEP in graphic interface mode one has to\n"
                "update your operation system and recompile CalcHEP\n"
                "PS: more likely your directory /usr/include/X11/ is empty\n");
         exit(1);       
      }    
   }
   if(blind==2) printf("\n\"");
   return err;
}

void finish(void){ if(blind!=1) crt0_finish(); if(blind==2) printf("\"\n");}

void goto_xy(int x,int y) {  xc = x-1; yc = y-1; }

void clr_scr(int color, int bkcolor)
{
   graphComPtr p;

   if(texflag) return;
   while (graphList !=NULL)
   { p=graphList;
     graphList=graphList->next;
     if(p->text != NULL) free(p->text);
     free(p);
   }
   scrcolor(color,bkcolor);
   tg_setlinestyle(SolidLn,NormWidth);
   tg_settextjustify(LeftText,TopText);
   if(blind!=1)
   {
     sg_screenSize(&bgi_maxx,&bgi_maxy);
     tg_setviewport(0,0,bgi_maxx,bgi_maxy);
     tg_bar(0,0,bgi_maxx,bgi_maxy);
   }
   xc =0;
   yc = 0;
   bgi_x=0;
   bgi_y=0;
}

int where_x(void) { return (int)(xc + 1); }

int where_y(void) {  return (yc + 1); }

void clr_eol(void) { clrbox(xc+1,yc+1,maxCol(),yc+1); }

static void print0(int xc,int yc, int color,int bkcolor, char * s) /* by A.Pukhov*/
{char buff[STRSIZ];
 int i,x,y,h1,w1;
 if (xc/CdX>maxCol())return;
 
 strcpy(buff,s);
 if(strlen(s)+xc/CdX>maxCol()) buff[maxCol()-xc/CdX]=0;

 for(i=0;i<strlen(buff);i++) {if (buff[i]<9) buff[i]=' ';}

 crt0_puts(xc,yc,color,bkcolor,buff);

   h1=CdY/2;
   w1=CdX/2;
   x=xc;
   y=yc;

 for(i=0;i<strlen(buff);i++)
 {  x=xc+i*CdX;

    switch(s[i])
    {
      case 1:
	sg_drawLine(x+w1,y+CdY,x+w1,y+h1,color,1,SolidLn);
	sg_drawLine(x+w1,y+h1,x+CdX,y+h1,color,1,SolidLn);
	break;
      case 2:
      case 6:
	sg_drawLine(x   ,y+h1,x+CdX,y+h1,color,1,SolidLn);
	break;
      case 3:
	sg_drawLine(x   ,y+h1,x+w1,y+h1,color,1,SolidLn);
	sg_drawLine(x+w1,y+h1,x+w1,y+CdY,color,1,SolidLn);
	break;
      case 4:
      case 8:
	sg_drawLine(x+w1,y   ,x+w1,y+CdY,color,1,SolidLn);
	break;
      case 5:
	sg_drawLine(x+w1,y   ,x+w1,y+h1,color,1,SolidLn);
	sg_drawLine(x+w1,y+h1,x   ,y+h1,color,1,SolidLn);
	break;
      case 7:
	sg_drawLine(x+w1,y   ,x+w1,y+h1,color,1,SolidLn);
	sg_drawLine(x+w1,y+h1,x+CdX,y+h1,color,1,SolidLn);
    }
  }
}

void get_text(int x1,int y1,int x2,int y2,void ** dump) /* by A.pukhov */
{
 int dn;
 graphComPtr p,q;

 x1=(x1-1)*CdX;
 y1=(y1-1)*CdY;
#ifdef _WIN32
  x2=x2*CdX+1;
#else
  x2=x2*CdX-1;
#endif
 y2=y2*CdY-1;


 p=(graphComPtr)m_alloc(sizeof(graphCommand));
 p->x1=x1;
 p->y1=y1;
 p->x2=x2;
 p->y2=y2;
 p->text=NULL;
 *dump=(void *)p;

 q=graphList;

 while(q!=NULL)
 {
   if(q->type==t_xt && q->y2 <= y2 && q->y1 >= y1 )
   {
     if (q->x1<x2 && q->x2>x1)
     {
        p->next=(graphComPtr)m_alloc(sizeof(graphCommand));
        p=p->next;
        (*p)=(*q);
        p->text=(char*)m_alloc(strlen(q->text)+1);
        dn=(x1-p->x1)/CdX;
        if(dn<0) dn=0;
        strcpy(p->text,q->text +dn);
        p->x1=q->x1+dn*CdX;
        dn=(q->x2-x2)/CdX;
        if (dn>0)
        { p->text[strlen(p->text)-dn]=0;
          p->x2 -= dn*CdX;
        }
     }
   }
   q=q->next;
 }

 p->next=NULL;
/* printf(" OK\n"); */
}

void put_text(void ** dump)  /* by A.Pukhov */
{
   graphComPtr p,q;

   p= (graphComPtr)(*dump);
   if(p==NULL) return;
   repeatGrCom(p->x1,p->y1,p->x2,p->y2,0);

   q=p;
   p=p->next;
   free(q);

   while(p!=NULL)
   {
     print0(p->x1,p->y1,p->color,p->aux, p->text);
     addToGraphList(t_xt,p->x1,p->y1,p->x2,p->y2,p->color,p->aux, p->text);
     q=p;
     p=p->next;
     free(q->text);
     free(q);
   }
   *dump=NULL;
}

void del_text(void ** dump)
{
   graphComPtr p,q;
   p=( graphComPtr)(*dump);  
   while (p!=NULL)
   { q=p;
     p=p->next;
     if(q->text != NULL) free(q->text);
     free(q);
   }
}

int print(char * format, ...)  /* by A.Pukhov */
{  va_list args;
   char dump[STRSIZ], *d;
   va_start(args, format);
   vsprintf(dump,format,args);
   va_end(args);
   d = dump;

   for(;;)
   {  int q;
      char * d1 = strchr(d,'\n');
      if (d1) d1[0] =0;

      q = (int) strlen(d);
      if(q)
      {
         if(blind!=1)
         {
           print0(xc*CdX,yc*CdY,fColor,bColor,d);
           addToGraphList(t_xt,xc*CdX,yc*CdY,(xc+q)*CdX-1,(yc+1)*CdY-1,fColor,bColor,d);
         }
         xc += q;
      }

      if (d1 == NULL)  return 0;
      yc++;
      xc = 0;
      d = d1 + 1;
      if (!d[0])  return 0;
   }
}

void scrcolor(int f_col,int b_col)
{
   if (!nColors)
   {      if(f_col ==b_col)  { if(f_col >7) {f_col=White;b_col=White;}
                               else         {f_col=Black;b_col=Black;}
                             }
     else if (f_col >b_col)  { f_col=White; b_col=Black;}
     else                    { f_col=Black; b_col=White;}
   }
   fColor=f_col;
   bColor=b_col;
}

void clrbox(int x1,int y1,int x2,int y2)
{
  repeatGrCom((x1-1)*CdX,(y1-1)*CdY,x2*CdX-1,y2*CdY-1,0);
  xc = x1-1;
  yc = y1-1;
}

/* LaTeX   by A.Pukhov  */
double texX(double x)
{
  return (x+texxshift)*texxscale;
}

double texY(double y)
{
  return (texymax1-y)*texyscale;
}


/* Sound support */
void be_be(void){ if(blind!=1) crt0_beep();}

void clearTypeAhead(void)
{  if(blind!=1) while (crt0_keypressed()) crt0_inkey(); }


int  inkey(void)   /* by A.Pukhov */
{ int  key; 

  if(!inkeyString) for(;;) 
  {    
     key = crt0_inkey(); 
     if(blind==2) 
     { 
       if(strchr("{}[]\\", key) || key==KB_MOUSE ) continue; 
       if(key==KB_ESC)        { printf("}"); return key;}
       else if(key==KB_DOWN)  { printf("["); return key;}
       else if(key==KB_UP)    { printf("]"); return key;}
       else if(key==KB_ENTER) { printf("{"); return key;}
       else if(isprint(key))  { printf("%c",key); return key;}
       else { printf("\\%.2x",key); return key;}
       continue;
     }
      
     if((key == KB_MOUSE)&&mouse_info.but1!=2) continue;  
     break;
  } 
  else      
  { 
    key= inkeyString[0];  
    if(key==0) { if(blind==1) sortie(101); else inkeyString=NULL; }
    else   inkeyString++;   
  } 

  if(inkeyString)
  { 
    switch (key)
    { 
      case '}':   key=KB_ESC ;  break;
      case '[':   key=KB_DOWN;  break;
      case ']':   key=KB_UP;    break;
      case '{':   key=KB_ENTER; break;
      case '\\':   sscanf(inkeyString,"%2x",&key); inkeyString+=2; break;
    }
  }
  
/*
if(blind==1)
{
if (key>30)  f_printf(stderr,"inkey='%c'\n",key);
 else         f_printf(stderr,"inkey=%d\n",key);
}
*/
  return key;
}

int  escpressed(void)  /* by A.Pukhov */
{  int key;
   if(blind) return 0;
   if(crt0_keypressed()) key=crt0_inkey();else return 0;
   for(;;)
   {
     if ((key == KB_ESC)||(key == KB_BACKSP)||(key == 5) ) return 1;
     if(crt0_keypressed())  key=crt0_inkey(); else return 0;
   }
}
