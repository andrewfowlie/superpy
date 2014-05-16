#ifndef __OUT_SERV_
#define __OUT_SERV_

#include "s_files.h"
#include "physics.h"
#include "polynom.h"

extern int   outputLanguage;

extern void momentToString(char * moment, char * outstr);

extern void readDenominators(void);

extern void  rewritepolynom(void);
extern void  findPrtclNum (char * procName,int * prtclNum);
extern void  emitconvlow(int * prtclNum);
extern void  writeLabel(char  comment);
extern void DiagramToOutFile(vcsect * vcs, int label,char comment);  
extern void  makeOutput(  void (*startOutput)(int,int*,int),
                          void (*diagramOutput)(vcsect*, catrec* ),         
                          void (*endOutput)(int*)              
                       );
                       
extern void seekArchiv(long n);

extern void readvardef(FILE*f);
extern void clearvardef(void);

extern poly readBuff;
extern int readSize;

#define wrtBuffSize 2000                                                 
#endif
