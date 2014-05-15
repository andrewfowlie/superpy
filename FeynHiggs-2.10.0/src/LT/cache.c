/*
	cache.c
		caching of tensor coefficients in
		dynamically allocated memory
		this file is part of FeynHiggs
		last modified 21 Sep 12 th
*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "externals.h"


#if NOUNDERSCORE
#define cachelookup_ cachelookup
#define cachecopy_ cachecopy
#endif

#ifndef REALSIZE
#define REALSIZE 8
#endif

#ifndef BIGENDIAN
#define BIGENDIAN 0
#endif

#if REALSIZE > 8
#define CMPBITS 64
#define MSB (1-BIGENDIAN)
#else
#define CMPBITS 62
#define MSB 0
#endif


typedef long long dblint;

typedef unsigned long long udblint;

typedef struct { dblint part[REALSIZE/8]; } RealType;

typedef const RealType cRealType;

typedef struct { RealType re, im; } ComplexType;

typedef long long memindex;


static inline int SignBit(const dblint i)
{
  return (udblint)i >> (8*sizeof(i) - 1);
}


static inline memindex PtrDiff(const void *a, const void *b)
{
  return (char *)a - (char *)b;
}


#if CMPBITS <= 64

static dblint CmpPara(cRealType *para1, cRealType *para2, int n)
{
  const dblint mask = -(1ULL << (64 - CMPBITS));
  while( n-- ) {
    const dblint c = (mask & para1->part[MSB]) -
                     (mask & para2->part[MSB]);
    if( c ) return c;
    ++para1;
    ++para2;
  }
  return 0;
}

#else

static dblint CmpPara(cRealType *para1, cRealType *para2, int n)
{
  const dblint mask = -(1ULL << (128 - CMPBITS));
  while( n-- ) {
    dblint c = para1->part[MSB] - para2->part[MSB];
    if( c ) return c;
    c = (mask & para1->part[1-MSB]) - (mask & para2->part[1-MSB]);
    if( c ) return c;
    ++para1;
    ++para2;
  }
  return 0;
}

#endif


static void *Lookup(cRealType *para, double *base,
  void (*calc)(RealType *, cRealType *),
  const int npara, const int nval)
{
  typedef struct node {
    struct node *next[2], *succ;
    int serial;
    RealType para[2];
  } Node;

#define base_valid (int *)&base[0]
#define base_last (Node ***)&base[1]
#define base_first (Node **)&base[2]

  const int valid = *base_valid;
  Node **last = *base_last;
  Node **next = base_first;
  Node *node;

  if( last == NULL ) last = next;

  while( (node = *next) && node->serial < valid ) {
    const dblint i = CmpPara(para, node->para, npara);
    if( i == 0 ) return &node->para[npara];
    next = &node->next[SignBit(i)];
  }

  node = *last;

  if( node == NULL ) {
	/* The "RealType para[2]" bit in Node is effectively an extra
	   Complex for alignment so that node can be reached with
	   an integer index into base */
    assert( node = malloc(sizeof(Node) +
      npara*sizeof(RealType) + nval*sizeof(ComplexType)) );
    node = (Node *)((char *)node +
      (PtrDiff(base, &node->para[npara]) & (sizeof(ComplexType) - 1)));
    node->succ = NULL;
    node->serial = valid;
    *last = node;
  }

  *next = node;
  *base_last = &node->succ;
  *base_valid = valid + 1;

  node->next[0] = NULL;
  node->next[1] = NULL;

  memcpy(node->para, para, npara*sizeof(RealType));
  calc(&node->para[npara], para);

  return &node->para[npara];
}


memindex cacheindex_(cRealType *para, double *base,
  void (*calc)(RealType *, cRealType *),
  const int *pnpara, const int *pnval)
{
  ComplexType *val = Lookup(para, base, calc, *pnpara, *pnval);
  return PtrDiff(val, base)/(long)sizeof(ComplexType);
}


void cachecopy_(ComplexType *dest, cRealType *para, double *base,
  void (*calc)(RealType *, cRealType *),
  const int *pnpara, const int *pnval)
{
  ComplexType *val = Lookup(para, base, calc, *pnpara, *pnval);
  memcpy(dest, val, *pnval*sizeof *dest);
}

