#include "matrix.h"
#include<stdlib.h>
#include<stdio.h>

matrix initm(int size)
{matrix result;
int i,j;
result.size=size;
for (i=0;i<size;i++)
	{for (j=0;j<size;j++)
		{result.e[i][j].r=0; result.e[i][j].i=0;}
		}
return result;
}

vecteur initv(int size)
{vecteur result;
int i;
result.size=size;
for (i=0;i<size;i++)
	{result.e[i].r=0; result.e[i].i=0;}
	
return result;
}

matrix identite(int size)
{matrix result;
int i,j;
result.size=size;
for (i=0;i<size;i++)
	{for (j=0;j<size;j++)
		{result.e[i][j].r=0; result.e[i][j].i=0;}
	result.e[i][i].r=1;
	}
return result;
}



void affichem(matrix mat)
{int i,j;
for (i=0;i<mat.size;i++)
	{for (j=0;j<mat.size;j++)
		printf("%e ",mat.e[i][j].r);
	printf("\n");}
}

/* because not used
void affichev(vecteur vec)
{int i;
printf("(");
for (i=0;i<vec.size;i++)
	printf("%e+i%e, ",vec.e[i].r,vec.e[i].i);
	
	printf(")\n");
}
*/

matrix prodm(matrix A,matrix B)
{matrix result;
int i,j,k;
result.size=A.size;
for (i=0;i<A.size;i++)
	for (j=0;j<A.size;j++)
		{
		result.e[i][j].r=0;
		result.e[i][j].i=0;
		for (k=0;k<2;k++)	
			result.e[i][j]=somme(result.e[i][j],prod(A.e[i][k],B.e[k][j]));
		}
return(result);
}

matrix adj(matrix A)
{
int i,j;
matrix result;
result.size=A.size;
for (i=0;i<A.size;i++)
	for (j=0;j<A.size;j++)
		result.e[i][j]=conjug(A.e[j][i]);
return result;
}

matrix transp(matrix A)
{
int i,j;
matrix result;
result.size=A.size;
for (i=0;i<A.size;i++)
	for (j=0;j<A.size;j++)
		result.e[i][j]=A.e[j][i];
return result;
}

vecteur multg(vecteur A,matrix B)
{vecteur result;
int i,j;
result.size=A.size;
for (i=0;i<A.size;i++)
	{result.e[i].r=0; 
	result.e[i].i=0;
	for (j=0;j<A.size;j++)
		result.e[i]=somme(result.e[i],prod(A.e[j],B.e[j][i]));
		}
return result;
}

my_complex scal(vecteur A,vecteur B)
{my_complex result;
int i;
result.r=0; result.i=0;
for (i=0;i<A.size;i++)
	result=somme(result,prod(A.e[i],B.e[i]));
return result;
}

vecteur adjv(vecteur A)
{vecteur result;
int i;
result.size=A.size;
for(i=0;i<A.size;i++)
	result.e[i]=conjug(A.e[i]);
return result;
}
