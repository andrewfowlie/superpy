#ifndef __CRT0__
#define __CRT0__

extern int  crt0_start(char * window_name, char * icon_file,char * ini_file,
                        int * colors, void(*error)(void));
extern void  crt0_finish(void);
extern void  crt0_puts(int xc,int yc,int tColor,int bColor,  char * s);
extern void  crt0_beep(void);
extern int   crt0_inkey(void);
extern int   crt0_keypressed(void);
extern void  crt0_refresh(void);
extern void  crt0_charSize(int * dx,int *dy);

extern  void sg_drawBox(int x1,int y1,int x2,int y2, int color);
extern  void sg_drawLine(int x1,int y1,int x2,int y2, int color,int thickness,int style);

extern  void sg_outText(int x,int y,int color, char * s);
extern  void sg_textSize(char * s, int * dx, int * dy);
extern  void sg_screenSize(int * x,int * y);
extern  void sg_setClip(int x1,int y1,int x2,int y2);


struct mouse_info_struct
	{
	int but1,but2,but3,but4,but5;
	int x,y;
	int col,row;
	};
	
extern struct mouse_info_struct mouse_info;

/*

#define KB_DOWN      (-80)
#define KB_UP        (-72)
#define KB_LEFT      (-75)
#define KB_RIGHT     (-77)
#define KB_HOME      (-71)
#define KB_END       (-79)
#define KB_BACKSP    (  8)
#define KB_TAB       (  9)
#define KB_F1        (-59)
#define KB_F2        (-60)
#define KB_F3        (-61)
#define KB_F4        (-62)
#define KB_F5        (-63)
#define KB_F6        (-64)
#define KB_F7        (-65)
#define KB_F8        (-66)
#define KB_F9        (-67)
#define KB_F10       (-68)
#define KB_DC        (-83)
#define KB_IC        (-82)
#define KB_PAGED     (-81)
#define KB_PAGEU     (-73)
#define KB_ENTER     ( 13)
#define KB_ESC       ( 27)
#define KB_SIZE      (-2)
#define KB_MOUSE     (-3)

*/

#define KB_BACKSP    (  8)
#define KB_TAB       (  9)
#define KB_ENTER     ( 13)
#define KB_ESC       ( 27)


#define KB_F1        (129)
#define KB_F2        (130)
#define KB_F3        (131)
#define KB_F4        (132)
#define KB_F5        (133)
#define KB_F6        (134)
#define KB_F7        (135)
#define KB_F8        (136)
#define KB_F9        (137)
#define KB_F10       (138)

#define KB_DOWN      (139)
#define KB_UP        (140)
#define KB_LEFT      (141)
#define KB_RIGHT     (142)
#define KB_HOME      (143)
#define KB_END       (144)
#define KB_DC        (145)
#define KB_IC        (146)
#define KB_PAGED     (147)
#define KB_PAGEU     (148)
#define KB_SIZE      (149)
#define KB_MOUSE     (150)



/*  Foreground and background CRT color constants  */

enum CRT_COLORS
{
    Black=0,
    Blue,
    Green,
    Cyan,
    Red,
    Magenta,
    Brown,
    LightGray,

    DarkGray,
    LightBlue,
    LightGreen,
    LightCyan,
    LightRed,
    LightMagenta,
    Yellow,
    White
};



enum line_styles
{
  SolidLn   = 0,
  DottedLn  = 1,  
  DashedLn  = 3
};
                
enum line_widths
{
   NormWidth  = 1,
   ThickWidth = 3
};
                                                               
extern void (*xw_expose)(int x,int y, int width,int height);

                                                               
#endif
