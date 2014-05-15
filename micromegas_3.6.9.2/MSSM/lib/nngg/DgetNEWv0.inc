* Dget.F
* the four-point tensor coefficients
* this file is part of LoopTools
* last modified 16 Jun 04 th


#define QEXT DBLE


#include "lt.h"


	integer*8 function Dget(p1, p2, p3, p4, p1p2, p2p3,
     &    m1, m2, m3, m4)
	implicit none
	double precision p1, p2, p3, p4, p1p2, p2p3
	double precision m1, m2, m3, m4
	double precision p12, p13, p14, p23, p24, p34

	integer*8 cachelookup
	external cachelookup, Dijkl

	double precision para(10)

	double precision Dcache(2)
	common /Dbase/ Dcache

* The following data statement initializes two *pointers* to NULL
* (see cache.c).  This is sneaky but ok because IEEE says 0D0 is
* represented as 8 0-bytes.
	data Dcache /0, 0/

	para(1) = p1
	para(2) = p2
	para(3) = p3
	para(4) = p4
	para(5) = p1p2
	para(6) = p2p3
	para(7) = m1
	para(8) = m2
	para(9) = m3
	para(10) = m4

c       call wrtdarg(para)
c	   write(*,*) '   '
	p12=(p1p2-p1-p2)/2.0D0
	p23=(p2p3-p2-p3)/2.0D0
	p13=(p4-p1-p2-p3-p12*2D0-p23*2D0)/2D0
	
	p14=-(p1+p12+p13)
	p24=-(p12+p2+p23)
	p34=-(p13+p23+p3)
	
	
c	write(*,*) p1,p2,p3,p4
c	write(*,*) p12,p13,p14,p23,p24,p34
c	write(*,*) '12: ',p1+p2-2.0D0*p12, p1*p2-p12**2
c	write(*,*) '13: ',p1+p3-2.0D0*p13, p1*p3-p13**2
c	write(*,*) '14: ',p1+p4-2.0D0*p14, p1*p4-p14**2
c	write(*,*) '23: ',p2+p3-2.0D0*p23, p2*p3-p23**2
c	write(*,*) '24: ',p2+p4-2.0D0*p24, p2*p4-p24**2
c	write(*,*) '34: ',p3+p4-2.0D0*p34, p3*p4-p34**2
c	write(*,*) 'k12', p1*p2-p12**2
c	write(*,*) 'k13', p1*p2-p12**2+p1*p3-p13**2+2*(p1*p23-p12*p13)
c	write(*,*) 'k23', p1*p3-p13**2+p2*p3-p23**2+2*(p3*p12-p13*p23)
	Dget = cachelookup(para, Dcache, Dijkl, 10, 46)
	end

************************************************************************

	subroutine Dijkl(P, D)
	implicit none
	double precision P(10)
	double complex D(46)

	integer*8 Cget
	double complex D0
	external Cget, D0

	double complex Ccache(1)
	common /Cbase/ Ccache

	XREAL M11, M12, M13, M22, M23, M33
	common /matrix3/ M11, M12, M13, M22, M23, M33

	XREAL f1, f2, f3, det3, A1, AC
	XREAL p12,p13,p23
	double complex d1, d2, d3, B1, B2
	double complex DD1, DD2, DD3
	double complex DD00, DD11, DD12, DD13, DD22, DD23, DD33, DD111b
	double complex DD001, DD002, DD003, DD111, DD112, DD113, DD122
	double complex DD123, DD133, DD222, DD223, DD233, DD333, DD111a
	double complex c1, c2, c3, c4, c5
        double precision e1, e2, e3, e4 
	integer*8 C234, C134, C124, C123, indx

	f1 = .5D0*(QEXT(P(1)) + QEXT(P(5)) - P(2))
	f2 = .5D0*(QEXT(P(1)) + QEXT(P(4)) - P(6))
	f3 = .5D0*(QEXT(P(5)) + QEXT(P(4)) - P(3))
	M22 = QEXT(P(1))*P(4) - f2*f2
	M23 = f1*f2 - P(1)*f3
	det3 = 2*(P(5)*M22 - P(4)*f1*f1 + f3*(f1*f2 + M23))
	M11 = (QEXT(P(4))*P(5) - f3*f3)/det3
	M12 = (f2*f3 - P(4)*f1)/det3
	M13 = (f1*f3 - P(5)*f2)/det3
	M22 = M22/det3
	M23 = M23/det3
	M33 = (QEXT(P(1))*P(5) - f1*f1)/det3

	f1 = QEXT(P(1)) + QEXT(P(7)) - P(8)
	f2 = QEXT(P(5)) + QEXT(P(7)) - P(9)
	f3 = QEXT(P(4)) + QEXT(P(7)) - P(10)

	C234 = Cget(P(2), P(3), P(6), P(8), P(9), P(10))
	C134 = Cget(P(5), P(3), P(4), P(7), P(9), P(10))
	C124 = Cget(P(1), P(6), P(4), P(7), P(8), P(10))
	C123 = Cget(P(1), P(2), P(5), P(7), P(8), P(9))

	c1 = Ccache(C234 + cc0) + Ccache(C234 + cc1) +
     &    Ccache(C234 + cc2)
	c2 = Ccache(C234 + cc12) + Ccache(C234 + cc22) +
     &    Ccache(C234 + cc2)
	c3 = 2*(Ccache(C234 + cc1) + c2) + Ccache(C234 + cc11) -
     &    Ccache(C234 + cc22) + Ccache(C234 + cc0)
	c4 = Ccache(C234 + cc1) + Ccache(C234 + cc11) +
     &    Ccache(C234 + cc12)
	c5 = Ccache(C234 + cc111) + Ccache(C234 + cc112) +
     &    Ccache(C234 + cc11)

	D(dd0) = D0(P(1), P(2), P(3), P(4), P(5), P(6),
     &    P(7), P(8), P(9), P(10))



	p12 = 0.5*(P(5)-P(1)-P(2))
	p23 = 0.5*(P(6)-P(2)-P(3))
	p13 = 0.5*(P(4)+P(2)-P(5)-P(6))


*       Case 1: p1, p2 independent


c	if(dabs(p(1)*p(2)-p12**2) .gt. dabs(p(1)*p(3)-p13**2+2*(p(1)*p23-p12*p13)) ) then

	if(dabs(p(1)*p(2)-p12**2) .gt. dabs(p(1)*p(2)-p12**2+p(1)*p(3)-p13**2+2*(p(1)*p23-p12*p13)) ) then


	e1 = 1D0
	e2 = (p13*(P(2)+p12)-(P(1)+p12)*p23)/(P(1)*P(2)-p12*p12)
	e3 = (P(1)*(P(2)+p23)-p12*(p12+p13))/(P(1)*P(2)-p12*p12)
	e4 = e2+e3-1D0
        B1 = 1D0
        B2 = 0D0
	AC = P(1)*e2*(e2-1)+(P(1)+P(2)+2*p12)*e3*(e3-1D0)+(P(1)+p12)*2*e2*e3-e4*P(7)+e2*P(8)+e3*P(9)-P(10)
c	write(*,*) 'p1p2 indep'

	else
*       Case 2: p1, p3 independent


	
	e1 = - (P(1)*(P(2)+p23)-p12*(p12+p13))/(P(1)*(P(2)+P(3)+2*p23)-p12*p12-p13*p13-2*p12*p13)
	e2 = ((P(1)+p12)*(P(3)+p23)-p13*(P(2)+p12+p13+p23))/(P(1)*(P(2)+P(3)+2*p23)-p12*p12-p13*p13-2*p12*p13)
	e3 = -1D0
	e4 = e2-e1-1D0
	B1 = 0D0
	B2 = 1D0
	AC = P(1)*(e2*e2+e1*e1-2*e1*e2-e2+e1)+e1*(e1+1D0)*(P(2)+P(3)+2*p12+2*p13+2*p23)-2*e1*e2*(p12+p13)-P(9)+e2*P(8)-e1*P(10)-e4*P(7)
c	write(*,*) 'p1p3 indep'

	endif
	
	if(dabs(e1).lt.1.D-14) then
	do indx=cc2,cc222
	Ccache(C123+indx)=0D0
	enddo
	endif
	

	if(dabs(e2).lt.1.D-14) then
	do indx=cc2,cc222
	Ccache(C134+indx)=0D0
	enddo
	endif


	if(dabs(e3).lt.1.D-14) then
	do indx=cc2,cc222
	Ccache(C124+indx)=0D0
	enddo
	endif


	if(dabs(e4).lt.1.D-14) then
        do indx=cc2,cc222
        Ccache(C234+indx)=0D0
	enddo
	endif





	DD1 = (1/AC)*(e1*Ccache(C123 + cc1)-e3*Ccache(C124 + cc1)-e4*(Ccache(C234 + cc0)+Ccache(C234 + cc1)+Ccache(C234 + cc2)))

	DD2 = (1/AC)*(e1*Ccache(C123 + cc2)-e2*Ccache(C134 + cc1)+e4*Ccache(C234 + cc1))

	DD3 = (1/AC)*(-e2*Ccache(C134 + cc2)-e3*Ccache(C124 + cc2)+e4*Ccache(C234 + cc2))


	D(dd1) = DD1 + B1*e2*DD3 + B2*e2*DD2


	D(dd2) = B1*(DD2+e3*DD3)


	D(dd3)= B2*(DD3-e1*DD2)


c	write(*,*) 'dd1=', Ccache(C234 + cc2)
c	write(*,*) 'e1=', e1 
c	write(*,*) 'e2=', e2 
c	write(*,*) 'e3=', e3 
c	write(*,*) 'e4=', e4 
c	write(*,*) 'p12=', p(1)*p(2)-p12**2 
c	write(*,*) 'p13=', p(1)*p(3)-p13**2 
c	write(*,*) 'p23=', p(2)*p(3)-p23**2 
c	write(*,*) 'p13=', p(1)*p(3)-p13**2+2*(p(1)*p23-p12*p13)
c	write(*,*) 'seg=', P(1)*(P(2)+P(3)+2*p23)-p12*p12-p13*p13-2*p12*p13

  	call rotate3(D(dd1), D(dd2), D(dd3),
     &    Ccache(C134 + cc0) - Ccache(C234 + cc0) - f1*D(dd0),
     &    Ccache(C124 + cc0) - Ccache(C234 + cc0) - f2*D(dd0),
     &    Ccache(C123 + cc0) - Ccache(C234 + cc0) - f3*D(dd0))



	DD00 = (1/AC)*(e1*Ccache(C123 + cc00)-e2*Ccache(C134 + cc00)-e3*Ccache(C124 + cc00)+e4*Ccache(C234 + cc00))

	DD11 = (1/AC)*(e1*Ccache(C123 + cc11)-e3*Ccache(C124 + cc11)+e4*(Ccache(C234 + cc11)+2D0*Ccache(C234 + cc12)+Ccache(C234 + cc22)+2D0*Ccache(C234 + cc1)+2D0*Ccache(C234 + cc2)+Ccache(C234 + cc0)))

	DD12 = (1/AC)*(e1*Ccache(C123 + cc12)-e4*(Ccache(C234 + cc11)+Ccache(C234 + cc12)+Ccache(C234 + cc1)))

	DD13 = (1/AC)*(-e3*Ccache(C124 + cc12)-e4*(Ccache(C234 + cc12)+Ccache(C234 + cc22)+Ccache(C234 + cc2)))

	DD22 = (1/AC)*(e1*Ccache(C123 + cc22)-e2*Ccache(C134 + cc11)+e4*Ccache(C234 + cc11))
	
	DD23 = (1/AC)*(-e2*Ccache(C134 + cc12)+e4*Ccache(C234 + cc12))

	DD33 = (1/AC)*(-e2*Ccache(C134 + cc22)-e3*Ccache(C124 + cc22)+e4*Ccache(C234 + cc22))


	D(dd00) = DD00

	D(dd11) = DD11 + B1*(2D0*e2*DD13+e2*e2*DD33) + B2*(2D0*e2*DD12+e2*e2*DD22)

	D(dd12) = B1*(DD12+e3*DD13+e2*DD23+e2*e3*DD33)

	D(dd13) = B2*(DD13-e1*DD12+e2*DD23-e1*e2*DD22)

	D(dd22) = B1*(DD22+2D0*e3*DD23+e3*e3*DD33)
	
	D(dd23) = 0D0

	D(dd33) = B2*(DD33-2D0*e1*DD23+e1*e1*DD22)






c	D(dd00) = P(7)*D(dd0) + .5D0*(Ccache(C234 + cc0) +
c     &    D(dd1)*f1 + D(dd2)*f2 + D(dd3)*f3)
	call rotate3(D(dd11), d1, d2,
     &    c1 - f1*D(dd1) - 2*D(dd00),
     &    Ccache(C124 + cc1) + c1 - f2*D(dd1),
     &    Ccache(C123 + cc1) + c1 - f3*D(dd1))
	call rotate3(D(dd12), D(dd22), d3,
     &    Ccache(C134 + cc1) - Ccache(C234 + cc1) - f1*D(dd2),
     &    -(Ccache(C234 + cc1) + f2*D(dd2)) - 2*D(dd00),
     &    Ccache(C123 + cc2) - Ccache(C234 + cc1) - f3*D(dd2))
	call rotate3(D(dd13), D(dd23), D(dd33),
     &    Ccache(C134 + cc2) - Ccache(C234 + cc2) - f1*D(dd3),
     &    Ccache(C124 + cc2) - Ccache(C234 + cc2) - f2*D(dd3),
     &    -(Ccache(C234 + cc2) + f3*D(dd3)) - 2*D(dd00))
c	D(dd12) = .5D0*(D(dd12) + d1)
c	D(dd13) = .5D0*(D(dd13) + d2)
c	D(dd23) = .5D0*(D(dd23) + d3)





	DD001 = (1/AC)*(e1*Ccache(C123 + cc001)-e3*Ccache(C124 + cc001)-e4*(Ccache(C234 + cc001)+Ccache(C234 + cc002)+Ccache(C234 + cc00)))

	
	DD002 = (1/AC)*(e1*Ccache(C123 + cc002)-e2*Ccache(C134 + cc001)+e4*Ccache(C234 + cc001))

	DD003 = (1/AC)*(-e2*Ccache(C134 + cc002)-e3*Ccache(C124 + cc002)+e4*Ccache(C234 + cc002))

	DD111a = (1/AC)*(e1*Ccache(C123 + cc111)-e3*Ccache(C124 + cc111)-e4*(Ccache(C234 + cc111)+3D0*Ccache(C234 + cc112)+3D0*Ccache(C234 + cc122)))


	DD111b = -(1/AC)*(e4*(Ccache(C234 + cc222)+3D0*Ccache(C234 + cc11)+6D0*Ccache(C234 + cc12)+3D0*Ccache(C234 + cc22)+3D0*Ccache(C234 + cc1)+3D0*Ccache(C234 + cc2)+Ccache(C234 + cc0)))


	DD111 = DD111a + DD111b


	DD112 = (1/AC)*(e1*Ccache(C123 + cc112)+e4*(Ccache(C234 + cc111)+2D0*Ccache(C234 + cc112)+2D0*Ccache(C234 + cc11)+2D0*Ccache(C234 + cc12)+Ccache(C234 + cc1)+Ccache(C234 + cc122)))

	DD113 = (1/AC)*(-e3*Ccache(C124 + cc112)+e4*(2D0*Ccache(C234 + cc122)+Ccache(C234 + cc222)+2D0*Ccache(C234 + cc12)+2D0*Ccache(C234 + cc22)+Ccache(C234 + cc2)+Ccache(C234 + cc112)))


	DD122 = (1/AC)*(e1*Ccache(C123 + cc122)-e4*(Ccache(C234 + cc111)+Ccache(C234 + cc112)+Ccache(C234 + cc11)))
	
	DD123 = -(1/AC)*(e4*(Ccache(C234 + cc112)+Ccache(C234 + cc122)+Ccache(C234 + cc12)))

	DD133 = (1/AC)*(-e3*Ccache(C124 + cc122)-e4*(Ccache(C234 + cc122)+Ccache(C234 + cc222)+Ccache(C234 + cc22)))

	DD222 = (1/AC)*(e1*Ccache(C123 + cc222)-e2*Ccache(C134 + cc111)+e4*Ccache(C234 + cc111))

	DD223 = (1/AC)*(-e2*Ccache(C134 + cc112)+e4*Ccache(C234 + cc112))

	DD233 = (1/AC)*(-e2*Ccache(C134 + cc122)+e4*Ccache(C234 + cc122))

	DD333 = (1/AC)*(-e2*Ccache(C134 + cc222)-e3*Ccache(C124 + cc222)+e4*Ccache(C234 + cc222))


	D(dd001) = DD001 + B1*e2*DD003 + B2*e2*DD002

	D(dd002) = B1*(DD002 + e3*DD003)

	D(dd003) = B2*(DD003 - e1*DD002)

	D(dd111) = DD111 + B1*(3D0*e2*DD113 + 3D0*e2*e2*DD133 + e2*e2*e2*DD333) + B2*(3D0*e2*DD112+3D0*e2*e2*DD122+e2*e2*e2*DD222)

	D(dd112) = B1*(DD112+e3*DD113+2D0*e2*e3*DD133+e2*e2*DD233+e2*e2*e3*DD333+2D0*e2*DD123)

	D(dd113) = B2*(DD113-e1*DD112-2D0*e2*e1*DD122+e2*e2*DD223-e2*e2*e1*DD222+2D0*e2*DD123)

	D(dd122) = B1*(DD122+2D0*e3*DD123+e3*e3*DD133+e2*DD223+2D0*e2*e3*DD233+e2*e3*e3*DD333)

        D(dd123) = 0D0

	D(dd133) = B2*(DD133-2D0*e1*DD123+e1*e1*DD122+e2*DD233-2D0*e2*e1*DD223+e2*e1*e1*DD222)

	D(dd222) = B1*(DD222+3D0*e3*DD223+3D0*e3*e3*DD233+e3*e3*e3*DD333)

	D(dd223) = 0D0

	D(dd233) = 0D0

	D(dd333) = B2*(DD333-3D0*e1*DD233+3D0*e1*e1*DD223-e1*e1*e1*DD222)


c	write(*,*) 'd1=', D(dd1)
c	write(*,*) 'd2=', D(dd2)
c	write(*,*) 'd3=', D(dd3)
c	write(*,*) 'C234=', Ccache(C234+cc1)
c	write(*,*) 'd00=', D(dd00)
c	write(*,*) 'd11=', D(dd11)
c	write(*,*) 'd12=', D(dd12)
c	write(*,*) 'd13=', D(dd13)
c	write(*,*) 'd22=', D(dd22) 
c	write(*,*) 'd23=', D(dd23)
c	write(*,*) 'd33=', D(dd33)




	call rotate3(D(dd001), D(dd002), D(dd003),
     &    Ccache(C134 + cc00) - Ccache(C234 + cc00) - f1*D(dd00),
     &    Ccache(C124 + cc00) - Ccache(C234 + cc00) - f2*D(dd00),
     &    Ccache(C123 + cc00) - Ccache(C234 + cc00) - f3*D(dd00))
	call rotate3(D(dd111), D(dd112), D(dd113),
     &    -f1*D(dd11) - c3 - 4*D(dd001),
     &    Ccache(C124 + cc11) - f2*D(dd11) - c3,
     &    Ccache(C123 + cc11) - f3*D(dd11) - c3)
	call rotate3(D(dd122), D(dd222), D(dd223),
     &    Ccache(C134 + cc11) - Ccache(C234 + cc11) - f1*D(dd22),
     &    -(Ccache(C234 + cc11) + f2*D(dd22)) - 4*D(dd002),
     &    Ccache(C123 + cc22) - Ccache(C234 + cc11) - f3*D(dd22))
	call rotate3(D(dd133), D(dd233), D(dd333),
     &    Ccache(C134 + cc22) - Ccache(C234 + cc22) - f1*D(dd33),
     &    Ccache(C124 + cc22) - Ccache(C234 + cc22) - f2*D(dd33),
     &    -(Ccache(C234 + cc22) + f3*D(dd33)) - 4*D(dd003))
	call rotate3(d1, D(dd123), d2,
     &    c2 - f1*D(dd13) - 2*D(dd003),
     &    Ccache(C124 + cc12) + c2 - f2*D(dd13),
     &    c2 - f3*D(dd13) - 2*D(dd001))
c	D(dd113) = .5D0*(D(dd113) + d1)
c	D(dd133) = .5D0*(D(dd133) + d2)

	D(dd0000) = (P(7)*D(dd00) + .5D0*(Ccache(C234 + cc00) +
     &    f1*D(dd001) + f2*D(dd002) + f3*D(dd003)))/3D0 + 1/36D0
	D(dd0011) = (P(7)*D(dd11) + .5D0*(c3 + f1*D(dd111) +
     &    f2*D(dd112) + f3*D(dd113)))/3D0
	D(dd0012) = (P(7)*D(dd12) - .5D0*(c4 - f1*D(dd112) -
     &    f2*D(dd122) - f3*D(dd123)))/3D0
	D(dd0013) = (P(7)*D(dd13) - .5D0*(c2 - f1*D(dd113) -
     &    f2*D(dd123) - f3*D(dd133)))/3D0
	D(dd0022) = (P(7)*D(dd22) + .5D0*(Ccache(C234 + cc11)+
     &    f1*D(dd122) + f2*D(dd222) + f3*D(dd223)))/3D0
	D(dd0023) = (P(7)*D(dd23) + .5D0*(Ccache(C234 + cc12)+
     &    f1*D(dd123) + f2*D(dd223) + f3*D(dd233)))/3D0
	D(dd0033) = (P(7)*D(dd33) + .5D0*(Ccache(C234 + cc22)+
     &    f1*D(dd133) + f2*D(dd233) + f3*D(dd333)))/3D0

	c3 = c2 + Ccache(C234 + cc112) + Ccache(C234 + cc122)
	c2 = c2 + c3 - Ccache(C234 + cc2) + Ccache(C234 + cc122) +
     &    Ccache(C234 + cc222)
	c4 = c4 + c1 + c2 + c3 + c4 + c5 + Ccache(C234 + cc12)

	call rotate3(D(dd1111), D(dd1112), D(dd1113),
     &    c4 - f1*D(dd111) - 6*D(dd0011),
     &    Ccache(C124 + cc111) + c4 - f2*D(dd111),
     &    Ccache(C123 + cc111) + c4 - f3*D(dd111))
	call rotate3(d1, D(dd1123), D(dd1133),
     &    -(c2 + f1*D(dd113)) - 4*D(dd0013),
     &    Ccache(C124 + cc112) - c2 - f2*D(dd113),
     &    -(c2 + f3*D(dd113)) - 2*D(dd0011))
	call rotate3(D(dd1122), d2, D(dd1223),
     &    c5 - f1*D(dd122) - 2*D(dd0022),
     &    c5 - f2*D(dd122) - 4*D(dd0012),
     &    Ccache(C123 + cc122) + c5 - f3*D(dd122))
	call rotate3(D(dd1222), D(dd2222), D(dd2223),
     &    Ccache(C134 + cc111) - Ccache(C234 + cc111) - f1*D(dd222),
     &    -(Ccache(C234 + cc111) + f2*D(dd222)) - 6*D(dd0022),
     &    Ccache(C123 + cc222) - Ccache(C234 + cc111) - f3*D(dd222))
	call rotate3(D(dd1233), D(dd2233), d3,
     &    Ccache(C134 + cc122) - Ccache(C234 + cc122) - f1*D(dd233),
     &    -(Ccache(C234 + cc122) + f2*D(dd233)) - 2*D(dd0033),
     &    -(Ccache(C234 + cc122) + f3*D(dd233)) - 4*D(dd0023))
	call rotate3(D(dd1333), D(dd2333), D(dd3333),
     &    Ccache(C134 + cc222) - Ccache(C234 + cc222) - f1*D(dd333),
     &    Ccache(C124 + cc222) - Ccache(C234 + cc222) - f2*D(dd333),
     &    -(Ccache(C234 + cc222) + f3*D(dd333)) - 6*D(dd0033))
	D(dd1113) = .5D0*(D(dd1113) + d1)
	D(dd1222) = .5D0*(D(dd1222) + d2)
	D(dd2333) = .5D0*(D(dd2333) + d3)
	
c	D(dd11)=0D0
c	D(dd33)=0D0
c	D(dd22)=0D0

c	do indx=dd00,dd33
c	D(indx)=0D0
c	enddo
	
	end

************************************************************************

	subroutine rotate3(out1, out2, out3, in1, in2, in3)
	implicit none
	double complex out1, out2, out3, in1, in2, in3

	XREAL M11, M12, M13, M22, M23, M33
	common /matrix3/ M11, M12, M13, M22, M23, M33

	return
	
	out1 = dcmplx(
     &    M11*dble(in1) + M12*dble(in2) + M13*dble(in3),
     &    M11*dimag(in1) + M12*dimag(in2) + M13*dimag(in3))
	out2 = dcmplx(
     &    M12*dble(in1) + M22*dble(in2) + M23*dble(in3),
     &    M12*dimag(in1) + M22*dimag(in2) + M23*dimag(in3))
	out3 = dcmplx(
     &    M13*dble(in1) + M23*dble(in2) + M33*dble(in3),
     &    M13*dimag(in1) + M23*dimag(in2) + M33*dimag(in3))
	end

