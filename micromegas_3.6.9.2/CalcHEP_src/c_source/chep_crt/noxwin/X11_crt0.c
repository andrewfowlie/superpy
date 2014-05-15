#include<stdlib.h>
#include "crt0.h" 

struct mouse_info_struct mouse_info;
void (*xw_expose)(int x,int y, int width,int height) =NULL; 
int  crt0_start(char * window_name,char * icon_file,char * ini_file,
                         int*colors, void(*error)(void)){return 1;}     
void  crt0_finish(void){}
void  crt0_puts(int xc,int yc, int color,int bkcolor,char* s){ }
int crt0_keypressed(void) { return 0; }
int  crt0_inkey(void){return 0;}
void  crt0_beep(void) {}
void sg_drawBox(int x1, int y1, int x2, int y2, int color){  }
static void bgi_set_dash (unsigned short rp){}
void sg_drawLine(int x1, int y1, int x2, int y2, int color,int thickness, int style) {}
void sg_outText(int x, int y, int color, char*  s){ }
void sg_textSize(char* s, int* dx, int* dy){ }
void sg_screenSize(int*  x, int* y){}
void crt0_charSize(int * dx,int *dy) { }
void  sg_setClip(int x1,int y1,int x2,int y2){ }
