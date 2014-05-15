
#ifndef __SETS__
#define __SETS__
     
#define SETLEN 32

#define _E  (-1)           /* end of set marker */
#define UpTo  (-2)         /* UpTo is analog of .. in Pascal */

typedef struct  set { unsigned char field[SETLEN];} set;

extern int  set_first(set* a,int i);
extern set  set_constr( int i, ...);
extern int  set_in(int a,set sp);
extern set  set_or( set a, set b);
extern set  set_aun( set a, set b);
extern set  set_and( set a, set b);
extern int  set_eq0(set a);
extern int  set_eq(set a,set b);

extern void set_add1(set *a, int k);
extern void set_del1(set *a, int k);
extern void set_print(set s);

#endif
