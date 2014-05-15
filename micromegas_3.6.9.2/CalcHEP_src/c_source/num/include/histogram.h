#ifndef __HIST__
#define __HIST__

#include"file_scr.h"
#include"phys_val.h"

extern table histTab;
extern void editHist(void);
extern void showHist(int X, int Y,char *title);
extern int wrt_hist(FILE *nchan);
extern int rdr_hist(FILE *nchan);
extern int wrt_hist2(FILE *nchan, char * comment);
extern int rdr_hist2(FILE *nchan,char ** comment);
extern int clearHists(void);
extern void fillHists(double w,double*V);
extern int correctHistList(void);
extern int add_hist(FILE *f, char **procname);
extern void xUnit(char key, char * units);
/*
typedef  struct  histRec 
{ struct histRec * next;
  linelist  mother;
  long  nPoints;
  char key[2][4];               
  physValRec* pList[2];
  char title[2][50];
  double hMin[2],hMax[2];
  double f[900];
  double ff[900];   
} histRec;
*/

#endif
