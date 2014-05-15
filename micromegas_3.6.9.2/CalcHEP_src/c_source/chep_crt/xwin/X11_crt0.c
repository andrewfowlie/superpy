/*
 Copyright (C) 1997, Andrei Semenov 
*/
/*-------------------------------------------------*
 *               XWindows -version                 *
 *      by Andrew Semenov,  Moscow University      *
 *        The latest revision of 21.08.1994        *
 *	Restricted version for Small Graphics	   *
 *-------------------------------------------------*/

#include <math.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "crt0.h" 


#define CRT_X   80
#define CRT_Y   25
#define EXITCODE 80

struct mouse_info_struct mouse_info;

static int isProportional;

static Display *display;
static int	screen;
static Window  win;
static GC	gc;
static Colormap cmap;
static XFontStruct   *bgi_font;
static XFontStruct   *text_font;

static char  BGI_DEFAULT_FONT [200]="-adobe-courier-bold-r-normal--14-*"; 
static char   BGI_TEXT_FONT   [200]="-adobe-courier-bold-r-normal--14-*";
static int BeepOn=1;

void (*xw_expose)(int x,int y, int width,int height) =NULL; 

static int bgi_maxcolor=15;
static int bgi_maxx, bgi_maxy;
static  int     bgi_status = 0;
static int xw_graphic;
static unsigned long bgi_colors[16];
static int xw_fn_w, xw_fn_h;

static void(*xw_error_stat)(void)=NULL;

static int  read_ini(char * fname)
{
  char str[200],name[200],parm1[200],parm2[200];
  FILE * ini;
 
  ini=fopen(fname,"r");
  if (ini==NULL) return 0;

  while( fgets(str,195,ini) !=NULL )
  { name[0]=0; parm1[0]=0; parm2[0]=0;
    sscanf(str,"%s%s%s",name,parm1,parm2); 
         if(!strcmp(name,"font")) { strcpy(BGI_TEXT_FONT,parm1);strcpy(BGI_DEFAULT_FONT,parm1);}
    else if(!strcmp(name,"colors")){if (strcmp(parm1,"off")==0) bgi_maxcolor=0;} 
    else if(!strcmp(name,"sound")) {if(strcmp(parm1,"off")==0) BeepOn=0;}
  }
  fclose(ini);
  return 1;
}

static int  X_IO_Error_H(Display * display) 
{ if(xw_error_stat)(*xw_error_stat)(); else exit(EXITCODE); return 0;}
  
static int  X_Error_H( Display* display ,XErrorEvent* error_event) 
{return X_IO_Error_H(display) ;}

static void setpalette(int   colornum, int  color)
{
	unsigned int r,g,b;
	XColor xc;
	if( colornum > bgi_maxcolor)  return;
	xc.flags=DoRed|DoGreen|DoBlue;
	r=(color & 04 )/2 + (color & 040)/040;
	g=(color & 02 )   + (color & 020)/020;
	b=(color & 01 )*2 + (color & 010)/010;
	xc.red=   r << 14;
	if( r != 0 ) xc.red |= 0x3fff;
	xc.green= g << 14;
	if( g != 0 ) xc.green |= 0x3fff;
	xc.blue=  b << 14;
	if( b != 0 ) xc.blue |= 0x3fff;

        XAllocColor(display,cmap,&xc);
	bgi_colors[colornum]=xc.pixel;
}

static void  setallpalette(void)
{
   int i;
   int palette_colors[16]={0,1,2,3,4,5,036,7,070,071,072,073,074,075,076,077}; 
   if(!bgi_maxcolor)
   {
      for(i=0;i<8;i++)  bgi_colors[i]=BlackPixel(display,screen);
      for(i=8;i<16;i++) bgi_colors[i]=WhitePixel(display,screen);
   }
   else
   {
      for(i=1;i<=bgi_maxcolor;i++)setpalette(i,palette_colors[i]);
      for(i=bgi_maxcolor+1;i<15;i++) palette_colors[i]=palette_colors[i-bgi_maxcolor];
      bgi_colors[0]=BlackPixel(display,screen);
      bgi_colors[15]=WhitePixel(display,screen);
   }
}


static int bgi_0to2(char * window_name,char * icon_file)
{
	unsigned int    depth;
	XSizeHints	size_hints;
	char		*dsp_name = NULL;
	Pixmap icon_pixmap;
	XSetWindowAttributes setwinattr;   
	if((display = XOpenDisplay (dsp_name))==NULL) return -3;

	screen=DefaultScreen(display);
	cmap=DefaultColormap(display,screen);
	depth=DefaultDepth(display,screen);
	if((bgi_font=XLoadQueryFont(display,BGI_DEFAULT_FONT))==NULL)
	{
	   fprintf(stderr,"Warning Graphical font not found : %s\n",BGI_DEFAULT_FONT);
	   fprintf(stderr,"Font  \"fixed\" is used instead\n");
	   if((bgi_font=XLoadQueryFont(display,"fixed"))==NULL)
	   { fprintf(stderr,"Font \"fixed\" not found too. Program terminated\n");  
	     X_IO_Error_H(NULL);
	   }
	}

	if((text_font=XLoadQueryFont(display, BGI_TEXT_FONT))==NULL)
	{ 
	   fprintf(stderr,"Text font not found : %s\n",BGI_TEXT_FONT);
	   fprintf(stderr,"Font  \"fixed\" is used instead\n");
	   text_font=XLoadQueryFont(display,"fixed");
	}
	xw_fn_w=text_font->max_bounds.width;
	xw_fn_h=text_font->ascent+text_font->descent;
	isProportional=(xw_fn_w!=text_font->min_bounds.width);
	if (isProportional)
	{ 
	   fprintf(stderr, "\n Warning: Text font is proportional\n");
	}
	bgi_maxx=text_font->max_bounds.width*CRT_X-1;
	bgi_maxy=(text_font->ascent+text_font->descent)*CRT_Y-1;
	if(bgi_maxcolor && depth >= 3 )
	{
            bgi_maxcolor= (1 << depth)-1;
            if(bgi_maxcolor>15) bgi_maxcolor=15;
        }
	else  bgi_maxcolor=0;
	win=XCreateWindow(	display,
	    RootWindow(display,screen),
	    0,
	    0,
	    bgi_maxx+1,
	    bgi_maxy+1,
	    3,
	    depth,
	    InputOutput,
	    CopyFromParent,
	    0,
	    &setwinattr);

        {  int hx,hy;
           unsigned w,h;     
           if(XReadBitmapFile(display,win,icon_file,&w,&h,&icon_pixmap,&hx,&hy) 
              !=BitmapSuccess)icon_pixmap=None;
        } 

	size_hints.flags=PSize|PMinSize;
	size_hints.width=bgi_maxx+1;
	size_hints.height=bgi_maxy+1;
	size_hints.min_width=bgi_maxx+1;
	size_hints.min_height=bgi_maxy+1;
	XSetStandardProperties(	display,
	    win,
	    window_name,
	    window_name,
	    icon_pixmap,
	    NULL,
	    0,
	    &size_hints);
	XSelectInput( display,
	    win,
	    ExposureMask|
	    KeyPressMask|
	    StructureNotifyMask|
	    ButtonPressMask|
	    ButtonReleaseMask);

	XSetErrorHandler(X_Error_H);
	XSetIOErrorHandler(X_IO_Error_H);


	{  unsigned long valuemask = 0;
	   XGCValues values;
	   XEvent report;
	   gc=XCreateGC(display,win,valuemask,&values);
           XMapWindow(display,win);
           XMaskEvent(display,ExposureMask,&report) ;
        }

	XSetFont(display,gc,text_font->fid);
	XSetForeground(display,gc,BlackPixel(display,screen));  
        xw_graphic=0;          
        setallpalette();	
	XFillRectangle(display,win,gc,0,0,bgi_maxx+1,bgi_maxy+1);
	return 0;	
}

int  crt0_start(char * window_name,char * icon_file,char * ini_file,int*colors,
                void (*xw_error)(void))
{ char inif[200];
  xw_error_stat=xw_error;
  if(bgi_status==0)
  {  char * f;
     sprintf(inif,"%.190s",ini_file);
     f=strtok(inif,";");
     do{ if(read_ini(f)) f=NULL;else f=strtok(NULL,";");} while (f);
     if(bgi_0to2(window_name,icon_file))
     { fprintf(stderr,"X Windows not detected.\n"); X_IO_Error_H(NULL);}     
  }
  bgi_status=1;
  *colors=bgi_maxcolor;
  return 0;
}

void  crt0_finish(void)
{
  XUnmapWindow(display,win);
  XFreeGC(display,gc);
  XFreeFont(display,bgi_font);
  XCloseDisplay(display);
  bgi_status=0;
  xw_error_stat=NULL;
}

void  crt0_puts(int xc,int yc, int color,int bkcolor,char* s)
{ 
  if(xw_graphic){XSetFont(display,gc,text_font->fid);xw_graphic=0;}
  if(isProportional)
  { int i;
     XSetForeground(display,gc,bgi_colors[bkcolor]);
     XFillRectangle(display,win,gc,xc,yc,
            xw_fn_w*strlen(s),xw_fn_h);

    XSetForeground(display,gc,bgi_colors[color]);
    for(i=0;i<strlen(s);i++) 
    XDrawString(display,win,gc,xc+i*xw_fn_w,yc+text_font->ascent,s+i,1);  
  } else
  {
     XSetForeground(display,gc,bgi_colors[color]);
     XSetBackground(display,gc,bgi_colors[bkcolor]);
     XDrawImageString(display,win,gc,xc,yc+text_font->ascent,
     s,strlen(s));
  }  
}


static int sizeChanged=0;

int crt0_keypressed(void)
{

  XEvent report;

  while(XCheckMaskEvent(display,ExposureMask|StructureNotifyMask,&report))
  { 
     switch(report.type)
     {
	case Expose:

           if (xw_expose !=NULL) (*xw_expose)(report.xexpose.x,report.xexpose.y,
	       report.xexpose.width,report.xexpose.height);
	   break;
	case ConfigureNotify:
	   bgi_maxx=report.xconfigure.width-1;
	   bgi_maxy=report.xconfigure.height-1;		
           if (xw_expose !=NULL) (*xw_expose)(0,0,bgi_maxx,bgi_maxy);
           sizeChanged=1;		
     }

  }


	
  if( XCheckMaskEvent(display,ButtonPressMask|ButtonReleaseMask|KeyPressMask|
  VisibilityChangeMask| ResizeRedirectMask|PropertyChangeMask|ColormapChangeMask,
	    &report))
  {	   
	switch(report.type)
	{
	case KeyPress:
	case ButtonPress:
	case ButtonRelease:
		XPutBackEvent(display,&report);
		return 1;
	default:
	 return 0;	
	}

  } return sizeChanged;
}

int  crt0_inkey(void)
{
        if (sizeChanged) 
        { sizeChanged=0; return KB_SIZE;} 
	for(;;)
	{
		XEvent report;
		XNextEvent(display,&report);

		switch(report.type)
		{
		case KeyPress:
			{
				char c; 
				KeySym keysym;
				c=0;
      				XLookupString(&report.xkey,&c,1,&keysym,NULL);
				if(c!=0 && c<127 )  return c;
				switch(keysym)
				{
				case XK_Left: return KB_LEFT;
				case XK_Down: return  KB_DOWN;
				case XK_Up:   return   KB_UP ;
				case XK_Right:return KB_RIGHT;
				case XK_Home: return KB_HOME;
				case XK_Tab:  return KB_TAB; 

				case XK_F1:        return  KB_F1;
				case XK_F2:    	   return  KB_F2;
				case XK_F3:        return  KB_F3;
				case XK_F4:        return  KB_F4;
				case XK_F5:        return  KB_F5;
				case XK_F6:        return  KB_F6;
				case XK_F7:        return  KB_F7;
				case XK_F8:        return  KB_F8;
				case XK_F9:        return  KB_F9;
				case XK_F10:       return  KB_F10;
                                case XK_F11:       return  KB_ESC;
                                case XK_Delete: /*   return  KB_DC; */
                                case XK_F12: 
				case XK_BackSpace: return  KB_BACKSP;
				case XK_Insert:    return  KB_IC;
				case XK_Next:      return  KB_PAGED;
				case XK_Prior:     return  KB_PAGEU;
				case XK_Begin:     return  KB_HOME;
				case XK_End:       return  KB_END;

				default:  break;				}
			}
			break;

		case Expose:
                  if (xw_expose !=NULL) (*xw_expose)(report.xexpose.x,
                  report.xexpose.y,report.xexpose.width,report.xexpose.height);
		  break;

		case ButtonPress:
			mouse_info.x=report.xbutton.x;
			mouse_info.y=report.xbutton.y;
			mouse_info.col=mouse_info.x/xw_fn_w+1;
			mouse_info.row=mouse_info.y/xw_fn_h+1;
			mouse_info.but1=mouse_info.but2=mouse_info.but3=0;
			mouse_info.but4=mouse_info.but5=0;
			if(report.xbutton.button==Button1) mouse_info.but1= 1;
			if(report.xbutton.button==Button2) mouse_info.but2= 1;
			if(report.xbutton.button==Button3) mouse_info.but3= 1;
			if(report.xbutton.button==Button4) mouse_info.but4= 1;
			if(report.xbutton.button==Button5) mouse_info.but5= 1;
			return KB_MOUSE;
			
		case ButtonRelease:
			mouse_info.x=report.xbutton.x;
			mouse_info.y=report.xbutton.y;
			mouse_info.col=mouse_info.x/xw_fn_w+1;
			mouse_info.row=mouse_info.y/xw_fn_h+1;
			mouse_info.but1=mouse_info.but2=mouse_info.but3=0;
			mouse_info.but4=mouse_info.but5=0;
			if(report.xbutton.button==Button1) mouse_info.but1= 2;
			if(report.xbutton.button==Button2) mouse_info.but2= 2;
			if(report.xbutton.button==Button3) mouse_info.but3= 2;
			if(report.xbutton.button==Button4) mouse_info.but4= 2;
			if(report.xbutton.button==Button5) mouse_info.but5= 2;
			return KB_MOUSE;
			
		case ConfigureNotify:
			bgi_maxx=report.xconfigure.width-1;
			bgi_maxy=report.xconfigure.height-1;
			if (xw_expose !=NULL) (*xw_expose)(0,0,bgi_maxx,bgi_maxy);
			return KB_SIZE;
		}
	}

}

void  crt0_beep(void) {if(BeepOn)XBell(display,0);}

void sg_drawBox(int x1, int y1, int x2, int y2, int color)
{  
   XSetForeground(display,gc,bgi_colors[color]);
   XFillRectangle(display,win,gc,x1,y1,x2-x1+1,y2-y1+1);
}

static void bgi_set_dash (unsigned short rp)
{
        int i,j,k,offset,dn;
        unsigned short p;
        char bgi_dashes[16];
        
        p=0;
        for(i=0;i<16;i++)
        {
                p = p << 1;
                if( ( rp & 1 ) == 1)
                        p|=1;
                rp = rp >> 1;
        }
        i=0;
        while((p&1)==1 && i<16)
        {
                i++;
                p = p >> 1;
        }
        offset=i;

        dn=0;

        while(i<16)
        {
                j=0;
                k=0;
                while((p&1) == 0 && (i<16))
                {
                        i++;
                        j++;
                        p = p >> 1;
                };
                while((p&1) == 1 && (i<16))
                {
                        i++;
                        k++;
                        p = p >> 1;
                };
                if(i==16)
                        k+=offset;
                bgi_dashes[dn++]=j;
                bgi_dashes[dn++]=k;
        }
        XSetDashes(display,gc, offset==0 ? offset : 16-offset ,bgi_dashes,dn);
}



void sg_drawLine(int x1, int y1, int x2, int y2, int color,
					  int thickness, int style)
{
	if(thickness==NormWidth) thickness=0;else thickness=ThickWidth;
	if(style==SolidLn)
	XSetLineAttributes(display,gc,thickness,LineSolid,CapRound,JoinRound);
	else
	{
	  XSetLineAttributes(display,gc,thickness,LineOnOffDash,CapRound,JoinRound);
	  switch(style)
	  {
	     case DottedLn: bgi_set_dash((unsigned short)0x1111); break;
	     case DashedLn: bgi_set_dash((unsigned short)0x0F0F); break;
	  }
	}
	XSetForeground(display,gc,bgi_colors[color]);
        XDrawLine(display,win,gc,x1,y1,x2,y2);
}

void sg_outText(int x, int y, int color, char*  s)
{
  if (!xw_graphic){XSetFont(display,gc,bgi_font->fid);xw_graphic=1;} 
  x-=bgi_font->min_bounds.lbearing;
  y-=bgi_font->descent;
  XSetForeground(display,gc,bgi_colors[color]);
  XDrawString(display,win,gc,x,y,s,strlen(s));
}

void sg_textSize(char* s, int* dx, int* dy)
{
    *dy=bgi_font->ascent+bgi_font->descent;
    *dx=XTextWidth(bgi_font,s,strlen(s));	
}

void sg_screenSize(int*  x, int* y)
{ *x=bgi_maxx;
  *y=bgi_maxy;
}

void crt0_charSize(int * dx,int *dy)
{ *dx=xw_fn_w;
  *dy=xw_fn_h;
}


void  sg_setClip(int x1,int y1,int x2,int y2)
{
  XRectangle rects;
  rects.x=x1;
  rects.y=y1;
  rects.width=x2-x1+1;
  rects.height=y2-y1+1;
  XSetClipRectangles(display,gc,0,0,&rects,1,Unsorted);
}
