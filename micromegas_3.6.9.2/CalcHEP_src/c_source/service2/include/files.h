#ifndef __FILES_
#define __FILES_
#include"syst.h"
    extern char * outputDir;
    extern char pathtocalchep[STRSIZ];
    extern void  copyfile(char *   namefrom, char *  nameto);
    extern void nextFileName(char* f_name,char * firstname,char * ext);
    extern char * unixPath(char *name);
    extern int findCalcHEPfile(char * name);
#define f_slash '/'
#endif
