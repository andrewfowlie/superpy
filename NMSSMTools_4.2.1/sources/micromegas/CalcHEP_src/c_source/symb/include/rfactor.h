#ifndef __RFACTOR_
#define __RFACTOR_
#include"lnum.h"
#include"procvar.h"

#define vmptr struct vmrec * 
typedef struct vmrec 
   { 
      vmptr        next; 
      char         name[VAR_NAME_SIZE_EXT];
      unsigned         deg;
   }  vmrec;
#undef vmptr
typedef struct vmrec *vmptr; 

typedef struct s_monom 
   { 
      NUM_TYPE      c; 
      vmptr        v; 
   }  s_monom; 
typedef s_monom *smptr; 

typedef struct r_monom 
   { 
      s_monom      n, d; 
   }  r_monom; 
typedef r_monom *rmptr; 

#define s_listptr struct s_listrec * 
typedef struct s_listrec 
   { 
      s_listptr    next; 
      s_monom      monom; 
   }  s_listrec;
#undef s_listptr
typedef struct s_listrec *s_listptr; 


extern void  eraseslist(s_listptr  s);

extern void  clrvm(vmptr  c);

extern void  sew_vm(vmptr *  p1,
                    vmptr    p2,
                    int  mlt);

extern void  reduce_s(s_monom *  s1,
                      s_monom *  s2);

extern void  mult_s(s_monom *  s1,
                    s_monom *  s2);

extern void  mult_rptr(rmptr *  m1,
                       rmptr *  m2);

extern char * smonomtxt(s_monom  s);

extern char * rmonomtxt(r_monom  r);

extern void  copysmonom(s_monom   src,
                        s_monom * dst);

extern void  diagramsrfactors(hlpcsptr     gst,
                              s_listptr *  s,
                              rmptr     *  totf);

extern void*  read_rmonom(char *  txt);
#endif
