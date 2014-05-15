#ifndef  __MATRIX__
#define  __MATRIX__

#include "my_complex.h"

//typedef struct{ double r;  double i; } my_complex;

typedef struct
	{ int size; 
           my_complex e[4][4];
        }  matrix;

typedef struct
	{int size;
	my_complex e[4];}
	vecteur;

extern matrix identite(int size);

extern matrix initm(int size);

extern vecteur initv(int size);

extern void affichem(matrix mat);

extern matrix adj(matrix A);

extern matrix transp(matrix A);

extern matrix prodm(matrix A,matrix B);

extern vecteur multg(vecteur A,matrix B);

extern my_complex scal(vecteur A,vecteur B);   

extern vecteur adjv(vecteur A);
#endif
