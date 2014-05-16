#ifndef __PHYS_VAL__
#define __PHYS_VAL__

extern  double calcPhysVal(char key,char * lv,double*V);
extern  int  checkPhysVal(char * name, char * key, char *plist);

typedef   struct physValRec
{ struct  physValRec * next;
  char    pstr[30];
}physValRec;

extern void cleanPVlist(physValRec * p);
extern int  checkPhysValN(char * name, char * key, physValRec **pLists);  

#endif
