#ifndef __FILE_SCR_
#define __FILE_SCR_

#include "syst.h"

typedef struct linerec
   {
		struct linerec * next;
		struct linerec * pred;
		char         line[STRSIZ];
	}  linerec;
typedef struct linerec * linelist;


typedef struct table
	{  char   mdlName[30];
           char   headln[80];
           char   format[STRSIZ];
           linelist   strings;
           int pos;
	}  table;


extern char errorText[256];
	extern void  showtext (int x1, int y1, int x2, int y2,
        char* headline,	FILE * fileptr);
        extern int   readtable0(table * tab, FILE* f);
        extern void  writetable0(table * tab,FILE* f);

	extern int  readtable(table * tab, char* fname);
	extern void  writetable(table * tab,char* fname);
	extern void  cleartab(table * tab);

#endif
