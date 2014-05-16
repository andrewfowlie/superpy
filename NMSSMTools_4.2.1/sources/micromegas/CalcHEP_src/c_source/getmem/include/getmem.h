#ifndef __GETMEM_
#define __GETMEM_

#include<stdlib.h>
#include<stdio.h>
#define blocksize  8192


typedef struct  marktp { void * blk_; unsigned pos_;} marktp;

extern void *  getmem_(unsigned  size);
extern void    release_(marktp * mrk);
extern void    mark_(marktp * mrk);
extern int     blockrest(int size); 


extern void    (*memerror) (void);
extern long    usedmemory;


extern   void (*memerror) (void);
extern   void * m_alloc(size_t size);
extern   void * re_alloc(void * prt,size_t size);

#endif
