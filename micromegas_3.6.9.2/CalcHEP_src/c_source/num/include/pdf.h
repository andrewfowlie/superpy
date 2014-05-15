
#ifndef __PDF__ 
#define __PDF__

#define NTPYMX 4
#define NGRMAX 9

typedef struct { char SFNAME[100][NGRMAX][NTPYMX][8]; } w505110_str;
typedef struct { int NPGSMX[NGRMAX][NTPYMX];} w505120_str;
extern w505110_str w505110_;
extern w505120_str w505120_;

extern void pdfset_( char *, double *, int);
extern void structm_(double *x, double * scale, double * upv, double *dnv,
double * usea, double *dsea, double *str, double *chm, double *bot, double *top,
double *gl); 

extern double  alphas2_(double *);
#endif
