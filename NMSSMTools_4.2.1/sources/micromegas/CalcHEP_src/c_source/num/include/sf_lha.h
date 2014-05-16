#ifndef __SF_LHA__
#define __SF_LHA__

extern int     p_lha(int* pNum);
extern void    n_lha(int i, char *name);
extern int     r_lha(int i, char *name);
extern int     m_lha(int i,int *);
extern int     mc_lha(int i);
extern int     init_lha(int i,double * be, double * mass);
extern double  c_lha(int i, double x,double q);

#endif
