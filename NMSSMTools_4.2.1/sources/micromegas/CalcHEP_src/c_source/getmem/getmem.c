/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include "getmem.h"
#include "syst.h"
#define ALIGNMENT 

long  usedmemory = 0;



#define markreleaseptr struct markreleaserec *
typedef struct markreleaserec
   {
		markreleaseptr    next;
		void *            dummy;
      unsigned char              mfield[blocksize];
   }  markreleaserec;

#undef markreleaseptr
typedef struct markreleaserec *markreleaseptr;

static markreleaseptr   blockptr = NULL;
static unsigned             currentpos=0;

int blockrest(int  size)
{
#ifdef ALIGNMENT
	size= ((size+7)>>3)<<3;
#endif
	return (blocksize - currentpos)/size ;
}



void  mark_(marktp*  mrk)
{
	mrk->blk_  = (void *)blockptr;
	mrk->pos_  = currentpos;
}

void release_(marktp* mrk)
{
  markreleaseptr   q;
	while (mrk->blk_ != (void *)blockptr)
	{
		q=blockptr;
		blockptr=q->next;
		free(q);
		usedmemory=usedmemory-sizeof(markreleaserec);
	}
	currentpos=mrk->pos_;
}


void *  getmem_(unsigned size)
{  markreleaseptr   q;
	void *          p;
#ifdef ALIGNMENT
	size= ((size+7)>>3)<<3;
#endif
	if( (blockptr == NULL) || (currentpos + size >= blocksize))
   {
      q = blockptr;
      blockptr = (markreleaseptr)m_alloc(sizeof(markreleaserec));
      if(blockptr == NULL)
      { blockptr=q;
        fprintf(stderr,"NOT ENOUGH MEMORY!\n");
        return NULL;
      }
      blockptr->next = q;
      currentpos = 0;
      usedmemory=usedmemory+sizeof(markreleaserec);
   }
   p = (void *)&(blockptr->mfield[currentpos]);
   currentpos += size;
   return p;
}


void  (*memerror)(void)  = NULL;

void *  m_alloc(size_t size)
{  void * p;
   p=malloc(size);
   if(!p) { if(memerror) (*memerror)(); else sortie(70);}
   return p;
}

void *  re_alloc(void * prt, size_t size)
{ 
    if(!prt)  return m_alloc(size);else 
    {   void * p= realloc( prt,size); 
         
	if(!p) { if(memerror) (*memerror)(); else sortie(70);}
	return p;
    }
}
