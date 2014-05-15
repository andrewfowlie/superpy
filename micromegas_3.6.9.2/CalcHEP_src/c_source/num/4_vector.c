/*
Copyright (C) 1997, Slava Ilyin
*/

#include<math.h>
#include"4_vector.h"
#include"nType.h"

//REAL pvect[400];  /*was [4][100] */

void lvtonv(char *lv, int nin, int nv, REAL * V)
{
  int i,n;
  REAL*q=V+4*(nv-1);
  for(i=0;i<4;i++) q[i]=0;
  for(i=0;n=lv[i] ;i++) if(n>nin) vsum4(nv,n,nv,1,V); else vsum4(nv,n,nv,-1,V);
} 

/* ****************************************** */
/*    Scalar product of two 4-vectors:     * */
/* ****************************************** */
REAL vdot4(int i,int j,REAL*V)
{
  i=4*i-4;
  j=4*j-4;
  return  V[i]*V[j]-V[i+1]*V[j+1]-V[i+2]*V[j+2]-V[i+3]*V[j+3];             
} 

/* ******************************************************* */
/*       SUM or Difference of two 4-vectors:            *  */
/* ISG=1  sum   and  ISG=-1   difference                *  */
/*             P(I) + ISG*P(J)=>P(K)                    *  */
/* ******************************************************* */

void vsum4(int i, int j, int k, int isg,REAL*V)
{   int l;
    i = 4*i - 4;
    j = 4*j - 4;
    k = 4*k - 4;
    if (isg == 1) {for(l=0;l<4;l++) V[k+l] = V[i+l] + V[j+l];}
    else          {for(l=0;l<4;l++) V[k+l] = V[i+l] - V[j+l];}
} /* vsum4_ */

/* ****************************************** */
/* NULLification of a 4-vector:  P(I) => 0 * */
/* ****************************************** */



void eps4(int n1, int  n2, int n3, int n4,REAL*V)
{
    REAL a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, 
	    a33, d1021, d1022, d1023, d1122, d1123, d1223;

n1*=4; n2*=4; n3*=4; n4*=4;

    a10 = V[n1 - 4];
    a11 = V[n1 - 3];
    a12 = V[n1 - 2];
    a13 = V[n1 - 1];
    a20 = V[n2 - 4];
    a21 = V[n2 - 3];
    a22 = V[n2 - 2];
    a23 = V[n2 - 1];
    a30 = V[n3 - 4];
    a31 = V[n3 - 3];
    a32 = V[n3 - 2];
    a33 = V[n3 - 1];
/*                               A10  A20  A30  X0 */
/*                               A11  A21  A31  X1 */
/*                               A12  A22  A32  X2 */
/*                               A13  A23  A33  X3 */
    d1021 = a10 * a21 - a20 * a11;
    d1022 = a10 * a22 - a20 * a12;
    d1023 = a10 * a23 - a20 * a13;
    d1122 = a11 * a22 - a21 * a12;
    d1123 = a11 * a23 - a21 * a13;
    d1223 = a12 * a23 - a22 * a13;
    V[n4 - 4] =   a31 * d1223 - a32 * d1123 + a33 * d1122;
    V[n4 - 3] =   a30 * d1223 - a32 * d1023 + a33 * d1022;
    V[n4 - 2] = -(a30 * d1123 - a31 * d1023 + a33 * d1021);
    V[n4 - 1] =   a30 * d1122 - a31 * d1022 + a32 * d1021;
} /* eps4_ */

void pvFill(REAL mass, REAL * mom4, int pos,REAL*V)
{
  int i,i0=4*(pos-1);
  V[i0]=mass;
  V[i0]*=V[i0];
  mass=mass*mass;

  for(i=1;i<4;i++) {V[i0+i]=mom4[i]; V[i0]+=V[i0+i]*V[i0+i];}
  V[i0]=sqrt(V[i0]);  
}



void incomkin(REAL m1, REAL m2, REAL p1, REAL p2, 
           REAL *sqrt_S_, REAL *Pcm_, REAL * rapidity_)
{
  REAL sqrt_S, Pcm,rapidity;
  REAL e1=sqrt(m1*m1+p1*p1);
  REAL e2=sqrt(m2*m2+p2*p2);
   
  sqrt_S=(e1+e2)*(e1+e2)-(p1-p2)*(p1-p2);
  
  rapidity= atanh((p1-p2)/(e1+e2));

  Pcm=p1*cosh(rapidity)-e1*sinh(rapidity);

  if(sqrt_S_) *sqrt_S_=sqrt(sqrt_S);
  if(Pcm_) *Pcm_=Pcm; 
  if(rapidity_) *rapidity_=rapidity;
}


