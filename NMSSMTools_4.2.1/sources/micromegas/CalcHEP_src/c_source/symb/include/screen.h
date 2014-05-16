#ifndef __SCREEN_
#define __SCREEN_

extern int  menulevel;

extern void menuhelp(void);
extern void modelinfo(void);
extern void processinfo(void);
extern void diagramsinfo(void);
extern void sq_diagramsinfo(void);
extern void viewsqdiagr(void);
extern void sqdiagrmenu(void);
extern void viewfeyndiag(int  del_mode);
extern int  viewresults(void);
extern void showheap(void);
extern char* currentModelName(void);

extern void f3_key_prog(int x);
extern void f4_key_prog(int x);
extern void f5_key_prog(int x);
extern void f6_key_prog(int x);
extern void f9_key_prog(int x);
extern void f10_key_prog(int x);

           						
extern void editModel(int edit);
extern void writeModelFiles( int l);

#endif
