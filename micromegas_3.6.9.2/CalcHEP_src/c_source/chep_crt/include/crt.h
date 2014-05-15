#ifndef __CRT_
#define __CRT_

#include "crt0.h"
#include <stdio.h>
/* interface variables   */

extern int  blind;
extern char * inkeyString;
extern int     windmin;
extern int     windmax;

extern char boxFrame[9];
extern int  nColors;
extern int  fColor,bColor;

#define FGmain  Black
#define BGmain  LightGray


extern   int print(char * format, ... );
extern int  start1(char * window_name,char * icon_file, char * ini_file, 
                    void(*error)(void));
extern void  finish(void);
extern void  goto_xy(int x,int y);
extern int   where_x(void);
extern int   where_y(void);
extern int   maxCol(void);   
extern int   maxRow(void); 
extern void  clr_scr(int color, int backgraund);
extern void  clr_eol(void);
extern void  chepbox(int x1,int y1,int x2,int y2);
extern void  scrcolor(int t_col,int b_col);
extern void  clrbox(int x1,int y1,int x2,int y2);
extern void  get_text(int x1,int y1,int x2,int y2,void ** dump);
extern void  put_text(void ** dump);
extern void del_text(void ** dump);
extern void  refresh_scr(void);
extern int setTexCharSize(char * font);


/* Horizontal and vertical justification for settextjustify */
enum pas_text_just
{
    LeftText   = 0,
    CenterText = 1,
    RightText  = 2,

    BottomText = 0,
    TopText    = 2
};

struct tg_linesettingstype {
	int linestyle;
	int thickness;
};


struct tg_textsettingstype {
	int horiz;
	int vert;
};


struct tg_viewporttype {
	int left, top, right, bottom;
};



extern  int     texflag;
extern  double  texxscale,texyscale;
extern  int     texxshift,texymax1;
extern  FILE *  out_tex;
extern  double texX(double);
extern  double texY(double);


extern  void tg_arrowline(int x1,int y1,int x2,int y2);

extern void tg_getviewsettings(struct tg_viewporttype * viewport);
extern int tg_textheight(char *textstring);
extern int tg_textwidth(char *textstring);
extern void  tg_bar( int left, int top, int right, int bottom );
extern void  tg_clearviewport( void );

extern int   tg_getmaxx( void );
extern int   tg_getmaxy( void );
extern void  tg_gettextsettings(struct tg_textsettingstype * text);
extern void   tg_line( int x1, int y1, int x2, int y2 );
extern void   tg_lineto( int x, int y );
extern void   tg_moveto( int x, int y );
extern void   tg_outtext( char *textstring );
extern void   tg_outtextxy( int x, int y, char *textstring );
extern void   putpixel( int x, int y, int color );
extern void   tg_rectangle( int left, int top, int right, int bottom );
extern void   tg_setlinestyle( int linestyle, int thickness );
extern void   tg_settextjustify( int horiz, int vert );
extern void   tg_setviewport( int left, int top, int right, int bottom);

extern void tg_getlinesettings (struct tg_linesettingstype * sls);
extern int   tg_getx(void);
extern int   tg_gety(void);
extern void  tg_linerel( int dx, int dy);
extern void  tg_moverel( int dx, int dy);


/* Sound support    */
extern void  be_be(void);


  extern int inkey(void);
  extern int escpressed(void);
  extern void clearTypeAhead(void);



#endif
