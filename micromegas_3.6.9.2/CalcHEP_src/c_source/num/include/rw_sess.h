#ifndef __RW_SESS__
#define __RW_SESS__ 
#include<stdio.h>

extern int w_sess__(FILE * mode);
extern int r_sess__(FILE * mode);

/*extern void inipar_(void); */
/*extern void clearSession(void); for interpreter */

extern int wrtprc_(void); 
extern int is_polarized(int k,int  nsub);
extern char* p_names[20];
extern int   p_codes[20];
#endif
