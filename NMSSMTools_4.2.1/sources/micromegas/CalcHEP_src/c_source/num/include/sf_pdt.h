#ifndef __SF_PDT__
#define __SF_PDT__

extern int     p_pdt(int* pNum);
extern void    n_pdt(int i, char *name);
extern int     r_pdt(int i, char *name);
extern int     m_pdt(int i,int*);
extern int     mc_pdt(int i);
extern int     init_pdt(int i,double * be, double * mass);
extern double  c_pdt(int i, double x,double q);
#endif
