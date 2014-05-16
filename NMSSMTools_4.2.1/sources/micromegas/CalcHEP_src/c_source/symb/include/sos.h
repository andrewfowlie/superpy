#ifndef __SOS_
#define __SOS_

extern void  save_sos(int  ercode);
extern void  restoreent(int * exitlevel);
extern void  saveent(int exitlevel);

extern int forceUG;
extern int newCodes;
#endif
