#ifndef __WRITEF__ 
#define __WRITEF__

extern int xpos;
extern void  wrt_0(char *s);
extern void  writeF(char * format,...);
extern FILE* outFile;
extern void  outFileOpen(char * format, ...);
extern void  outFileClose(void);

#endif
