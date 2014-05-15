!--------------------------------------------------------------------
! This file is part of HiggsSignals (TS 04/03/2013)
!--------------------------------------------------------------------
module numerics
!--------------------------------------------------------------------
!-- This module contains useful numerical functions inspired by the book
!-- NUMERICAL RECIPES in FORTRAN (http://www.haoli.org/nr/bookfpdf.html)

!--cof, stp are needed for gammln(xx)
double precision, dimension(6) :: cof = (/76.18009172947146d0, &
	-86.50532032941677d0,24.01409824083091d0,-1.231739572450155d0,&
	.1208650973866179d-2,-.5395239384953d-5/)
double precision, parameter :: stp = 2.5066282746310005d0

contains
!--------------------------------------------------------------------
function factrl(n)
! USES gammln
! Returns the value n! as a floating-point number.
!--------------------------------------------------------------------
integer, intent(in) :: n
double precision :: factrl
integer :: j,ntop 
REAL a(33)!	Table to be filled in only as required.
SAVE ntop,a 
DATA ntop,a(1)/0,1./!	Table initialized with 0! only.

if (n.lt.0) then
 write(*,*) "WARNING : negative factorial in factrl" 
else if (n.le.ntop) then	!Already in table.
 factrl=a(n+1) 
else if (n.le.32) then	!Fill in table up to desired value.
 do j=ntop+1,n
  a(j+1)=j*a(j)
 enddo
 ntop=n
 factrl=a(n+1)
else
! Larger value than size of table is required. 
! Actually, this big a value is going to overflow on many computers,
! but no harm in trying.
 factrl=exp(gammln(dble(n+1.)))
endif
end function factrl
!--------------------------------------------------------------------
function gammp(a,x)
! USES gcf,gser
! Returns the incomplete gamma function P(a,x). 
!--------------------------------------------------------------------
implicit none
double precision, intent(in) :: a,x
double precision :: gammp,gammcf,gamser,gln
if(x.lt.0..or.a.le.0.) write(*,*) "WARNING: bad arguments in gammp"
if(x.lt.a+1.)then		!	Use the series representation.
 call gser(gamser,a,x,gln)
 gammp=gamser 
else					!   Use the continued fraction representation
 call gcf(gammcf,a,x,gln)
 gammp=1.-gammcf		!	and take its complement. 
endif
end function gammp
!--------------------------------------------------------------------
function gammq(a,x)
! USES gcf,gser
! Returns the incomplete gamma function Q(a, x) ≡ 1 − P (a, x). 
!--------------------------------------------------------------------
implicit none
double precision, intent(in) :: a,x
double precision :: gammq,gammcf,gamser,gln
if(x.lt.0..or.a.le.0.) write(*,*) "WARNING: bad arguments in gammq"
if(x.lt.a+1.)then		!	Use the series representation
 call gser(gamser,a,x,gln)
 gammq=1.-gamser		!	and take its complement. 
else					!	Use the continued fraction representation.
 call gcf(gammcf,a,x,gln)
 gammq=gammcf 
endif
end function gammq
!--------------------------------------------------------------------
function gammln(xx)
!--Returns the value ln[Γ(xx)] for xx > 0.
!--------------------------------------------------------------------
implicit none
double precision, intent(in) :: xx
integer :: j
double precision :: gammln,ser,tmp,x,y
x=xx
y=x 
tmp=x+5.5d0 
tmp=(x+0.5d0)*log(tmp)-tmp 
ser=1.000000000190015d0 
do j=1,6
 y=y+1.d0
 ser=ser+cof(j)/y
enddo
gammln=tmp+log(stp*ser/x)
end function gammln
!--------------------------------------------------------------------
subroutine gser(gamser,a,x,gln)
! USES gammln
! Returns the incomplete gamma function P(a,x) evaluated by its series 
! representation as gamser. Also returns lnΓ(a) as gln. 
!--------------------------------------------------------------------
implicit none
double precision, intent(in) :: a,x
integer, parameter :: itmax=100
double precision, parameter :: eps=3.e-7
integer :: n
double precision :: ap,del,sum,gamser,gln 
gln=gammln(a) 
if(x.le.0.)then
 if(x.lt.0.) write(*,*) "WARNING: x < 0 in gser"
 gamser=0.
else
 ap=a
 sum=1./a
 del=sum
 do n=1,itmax
  ap=ap+1.
  del=del*x/ap 
  sum=sum+del 
  if(abs(del).lt.abs(sum)*eps) exit
  if(n.eq.itmax) write(*,*) "WARNING: a too large, ITMAX too small in gser"
 enddo
 gamser=sum*exp(-x+a*log(x)-gln)
endif
end subroutine gser
!--------------------------------------------------------------------
subroutine gcf(gammcf,a,x,gln)
! USES gammln
! Returns the incomplete gamma function Q(a, x) evaluated by its continued
! fraction representation as gammcf. Also returns ln Γ(a) as gln. 
!Parameters: 
!--ITMAX is the maximum allowed number of iterations;
!--EPS is the relative accuracy;
!--FPMIN is a number near the smallest representable floating-point number.
!--------------------------------------------------------------------
implicit none
double precision, intent(in) :: a,x
integer, parameter :: itmax=100
double precision, parameter :: eps=3.e-7
double precision, parameter :: fpmin=1.e-30
integer :: i
double precision :: an,b,c,d,del,h,gammcf,gln
gln=gammln(a)
b=x+1.-a
c=1./fpmin
d=1./b
h=d
do i=1,itmax
  an=-i*(i-a)
  b=b+2.
  d=an*d+b
  if(abs(d).lt.fpmin)d=fpmin
  c=b+an/c
  if(abs(c).lt.fpmin)c=fpmin
  d=1./d
  del=d*c
  h=h*del
  if(abs(del-1.).lt.eps) exit
  if(i.eq.itmax) write(*,*) "WARNING: a too large, ITMAX too small in gcf"
enddo
gammcf=exp(-x+a*log(x)-gln)*h
end subroutine gcf
!--------------------------------------------------------------------
function my_erf(x)
implicit none
double precision, intent(in) :: x
double precision :: my_erf
if(x.lt.0.) then
 my_erf=-gammp(0.5D0,x**2)
else 
 my_erf=gammp(0.5D0,x**2)
endif 
end function my_erf

!--------------------------------------------------------------------
subroutine invmatrix(A,Y)
! Calculates the inverse of matrix A and outputs it as Y.
!--------------------------------------------------------------------
 double precision, dimension(:,:), intent(inout) :: A
 double precision, allocatable, intent(out) :: Y(:,:)

 double precision, allocatable :: A1(:,:)
 double precision, allocatable :: temp(:) 
 integer, allocatable :: INDX(:) 
 integer d, rc, n, i, j
 character*12 input, output
 character*8 s
 double precision, parameter :: small = 1.0D-6
  
 n = size(A,dim=1) 
  
  allocate(A1(n,n),stat=ialloc)
  allocate(Y(n,n),stat=ialloc)
  allocate(temp(n),stat=ialloc)
  allocate(INDX(n),stat=ialloc)

  do i=1, n 
    do j=1, n      
	  A1(i,j) = A(i,j)
	  Y(i,j) = 0.d0
    end do
    Y(i,i) = 1.d0
  end do

!call LU decomposition routine (only once)
  call LUDCMP(A1,n,INDX,D,rc)

!call solver if previous return code is ok
!to obtain inverse of A1 one column at a time
  if (rc.eq.0) then
    do j=1, n
      call LUBKSB(A1,n,INDX,Y(:,J))
    end do
  endif
!the inverse matrix is now in matrix Y
!the original matrix A1 is destroyed

!print results or error message
  if (rc.eq.1) then
    write(*,*) ' The matrix is singular, no solution !'
!  else
!    write(*,*) ' '
!    write(*,*) '  Inverted matrix Y:'
!    write(*,*)
!    do i=1, n
!     write(*,'(10E14.6)') (Y(i,j),j=1,n)
!    enddo
  end if

!verify A x Y = I (result put in A) 
  call MATMULT(A,Y,A1,n,n)
!A should now contain identity matrix
!  write(*,*) ' '
!  write(*,*) '  Verification A*Y = I:'
!  write(*,*)

!  do i=1, n
!   write(*,*) (A1(i,j),j=1,n)
!  enddo
  do i=1,n
   if(abs(A1(i,i)-1).ge.small) then
!    write(*,*) 'WARNING in subroutine invmatrix: Difficulties in inversion!'
!    write(*,*) 'Deviation from 1 is:',(A1(i,i)-1)
   endif 
  enddo

  end subroutine invmatrix
!--------------------------------------------------------------------
!*******************************************************
!*    LU decomposition routines						   *
!*                                                     *
!*                 F90 version by J-P Moreau, Paris    *
!* --------------------------------------------------- *
!* Reference:                                          *
!*                                                     *
!* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
!*  University Press, 1986" [BIBLI 08].                *
!*                                                     * 
!*******************************************************
 Subroutine LUDCMP(A,N,INDX,D,CODE)
!--------------------------------------------------------------------
!  * Given an N x N matrix A, this routine replaces it by the LU *
!  * decomposition of a rowwise permutation of itself. A and N   *
!  * are input. INDX is an output vector which records the row   *
!  * permutation effected by the partial pivoting; D is output   *
!  * as -1 or 1, depending on whether the number of row inter-   *
!  * changes was even or odd, respectively. This routine is used *
!  * in combination with LUBKSB to solve linear equations or to  *
!  * invert a matrix. Return code is 1, if matrix is singular.   *
!--------------------------------------------------------------------
 PARAMETER(NMAX=100,TINY=1.5D-16)
 REAL*8  AMAX,DUM, SUM, A(N,N),VV(NMAX)
 INTEGER CODE, D, INDX(N)

 D=1; CODE=0

 DO I=1,N
   AMAX=0.d0
   DO J=1,N
     IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
   END DO ! j loop
   IF(AMAX.LT.TINY) THEN
     CODE = 1
     RETURN
   END IF
   VV(I) = 1.d0 / AMAX
 END DO ! i loop

 DO J=1,N
   DO I=1,J-1
     SUM = A(I,J)
     DO K=1,I-1
       SUM = SUM - A(I,K)*A(K,J) 
     END DO ! k loop
     A(I,J) = SUM
   END DO ! i loop
   AMAX = 0.d0
   DO I=J,N
     SUM = A(I,J)
     DO K=1,J-1
       SUM = SUM - A(I,K)*A(K,J) 
     END DO ! k loop
     A(I,J) = SUM
     DUM = VV(I)*DABS(SUM)
     IF(DUM.GE.AMAX) THEN
       IMAX = I
       AMAX = DUM
     END IF
   END DO ! i loop  
   
   IF(J.NE.IMAX) THEN
     DO K=1,N
       DUM = A(IMAX,K)
       A(IMAX,K) = A(J,K)
       A(J,K) = DUM
     END DO ! k loop
     D = -D
     VV(IMAX) = VV(J)
   END IF

   INDX(J) = IMAX
   IF(DABS(A(J,J)) < TINY) A(J,J) = TINY

   IF(J.NE.N) THEN
     DUM = 1.d0 / A(J,J)
     DO I=J+1,N
       A(I,J) = A(I,J)*DUM
     END DO ! i loop
   END IF 
 END DO ! j loop

 RETURN
 END subroutine LUDCMP
!--------------------------------------------------------------------
 Subroutine LUBKSB(A,N,INDX,B)
!--------------------------------------------------------------------
!  * Solves the set of N linear equations A . X = B.  Here A is     *
!  * input, not as the matrix A but rather as its LU decomposition, *
!  * determined by the routine LUDCMP. INDX is input as the permuta-*
!  * tion vector returned by LUDCMP. B is input as the right-hand   *
!  * side vector B, and returns with the solution vector X. A, N and*
!  * INDX are not modified by this routine and can be used for suc- *
!  * cessive calls with different right-hand sides. This routine is *
!  * also efficient for plain matrix inversion.                     *
!--------------------------------------------------------------------
 REAL*8  SUM, A(N,N),B(N)
 INTEGER INDX(N)

 II = 0

 DO I=1,N
   LL = INDX(I)
   SUM = B(LL)
   B(LL) = B(I)
   IF(II.NE.0) THEN
     DO J=II,I-1
       SUM = SUM - A(I,J)*B(J)
     END DO ! j loop
   ELSE IF(SUM.NE.0.d0) THEN
     II = I
   END IF
   B(I) = SUM
 END DO ! i loop

 DO I=N,1,-1
   SUM = B(I)
   IF(I < N) THEN
     DO J=I+1,N
       SUM = SUM - A(I,J)*B(J)
     END DO ! j loop
   END IF
   B(I) = SUM / A(I,I)
 END DO ! i loop

 RETURN
 END subroutine LUBKSB
!--------------------------------------------------------------------
  SUBROUTINE MATMULT(A,B,C,N,M)                                              
!*******************************************                                     
!*       MULTIPLY TWO REAL MATRICES        *
!* --------------------------------------- *                                     
!* INPUTS:    A  MATRIX N*N                *                                     
!*            B  MATRIX N*M                *                                     
!*            N  INTEGER                   *                                     
!*            M  INTEGER                   *                                     
!* --------------------------------------- *                                     
!* OUTPUT:    C  MATRIX N*M, PRODUCT A*B   *                                     
!*                                         *                                     
!******************************************* 
  REAL*8 A(N,N),B(N,M),C(N,M),SUM                                           
  DO I=1,N                                                                  
    DO J=1,M                                                                
      SUM=0.                                                                
      DO K=1,N                                                              
        SUM=SUM+A(I,K)*B(K,J)                                               
      ENDDO                                                                 
      C(I,J)=SUM                                                            
    ENDDO                                                                   
  ENDDO                                                                     
  RETURN                                                                    
  END subroutine matmult 
!--------------------------------------------------------------------
end module numerics
!--------------------------------------------------------------------
