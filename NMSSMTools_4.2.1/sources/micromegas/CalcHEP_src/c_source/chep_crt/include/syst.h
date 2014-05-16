#ifndef __SYST_
#define __SYST_

#ifndef STRSIZ
#define STRSIZ 4096
#endif



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include "getmem.h"

#ifndef MAX
#define MAX(a,b)        (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b)        (((a) < (b)) ? (a) : (b))
#endif

#ifndef ABS
#define ABS(x)          (((x) < 0 ? -(x) : (x)))
#endif


#ifndef SEEK_SET 
#include<unistd.h>
#endif 

extern   void (*diskerror) (void);
extern   int f_printf(FILE *fp,char * format, ... );
extern   size_t f_write(void *ptr,size_t size,size_t n,FILE *fp);

extern  char * trim(char *);

#define strlen(x) (long)strlen(x)

    extern int  setLockFile(char * fname);
    extern void unLockFile(int id); 
    extern int writeLockFile(char * fname);
    extern void sortie(int code);  


#endif
