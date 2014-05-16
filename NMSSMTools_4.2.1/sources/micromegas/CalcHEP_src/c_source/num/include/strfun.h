#ifndef __STRFUN__
#define __STRFUN__

extern int sf_menu(int i);
extern double strfun_(int i, double x, double q);

extern int rd_sf__(FILE *mode);
extern int wrt_sf__(FILE *mode);
extern int loadStrFun(char *  name1, char*name2);

extern int initStrFun(int mode);

extern void strFunName(int i, char * mess);

extern double sf_mass[2];
extern double sf_be[2];
extern int sf_num[2];

#endif
