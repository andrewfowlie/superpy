#ifndef __CRT_UTIL_
#define __CRT_UTIL_

#include "crt.h"

extern void  menu0(int col,int row,char* label, char* filename,
        void (**funcKey)(int i) , char** funcKeyMess, void **  hscr, int* kk);

extern int   mess_y_n(int    x1,int    y1,	char  * txtstr);
extern void  messanykey(int    x1,	int    y1, char  * txtstr);
extern int   yesnokey(void);
extern int   informline(long curent,long total);
extern int   str_redact(char * txt, int  npos, int maxLen);

extern int   correctDouble(int x,int y,char* txt,double * var,int clear);
extern int   correctLong(int x,int y,char* txt,long   * var,int clear);
extern int   correctInt(int x,int y,char* txt,int* var,int clear);
extern int   correctStr(int x,int y,char* txt,char* var, int maxLen, int clear);

extern void  menu1(int col,int row,char * name,  char* menustr,
                                    char * help, void*hscr,int * kk);
extern int   improveStr(char * str, char * mark, char * format, ...);

extern char *f3_mess[8] ;
extern void (*f3_key[8])(int);


#endif
