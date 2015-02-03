c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: B_FUN.F
c     Revised: 28: 3:1994(J.R.)
c     Next order of F expansion around s=0 added
c     Revised:  8: 6:2013(J.R.)
c     New common /cmem/ added to comply with changes in c_fun.f
c     Revised:  25: 2:2014(J.R.)
c     Additional order of s/m^2 added in f expansion

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains two-point standard loop integrals         c
c     a(m1), b0(s,m1,m2), b1(s,m1,m2), b21(s,m1,m2), b22(s,m1,m2)  c
c     These functions are defined as in: A.Axelrod,                c
c     Nucl.Phys.B209(1982)p.349 except that b functions have       c
c     opposite signs (compare: Chankowski,Pokorski,Rosiek,         c
c     Nucl.Phys.B423(1994)p.437 available as hep-ph/9303309.       c           
c     The parameter del in the common/renorm/ helps to check the   c
c     cancelation of the divergences by counterterms. The  amiu2   c
c     parameter in the common/renorm/  is the 't Hooft mass        c
c     unit which also should cancel in the final expressions       c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function a(x)
c     A0 Veltman function
      implicit double precision (a-h,o-z)
      logical infstat
      common/renorm/del,amiu2,infstat
      if (infstat) then
        a = - x*x*del
        return
      end if
      if (x.eq.0.d0) then
        a = 0
      else
        a = - x*x*(del + 1 - log(x*x/amiu2))
      end if
      return
      end

      double complex function f(s,a1,a2)
c     F function - general case
      implicit double precision (a-h,o-z)
      double complex feq,f0,fneq
      if ((a1.eq.0.d0).and.(a2.eq.0.d0)) stop '2 x ZERO MASSES IN F'
      if (a1.eq.a2) then
        f = feq(s,a1)
        return
      else if ((a1*a2).eq.0.d0) then
        f = f0(s,a1 + a2)
        return
      else
        f = fneq(s,a1,a2)
      end if
      return
      end

      double complex function feq(s,a1)
c     F function - equal masses
      implicit double precision (a-h,o-z)
      common/spence/ber(9),pi,pi6,eps
      x = 4*a1*a1/s
      if (abs(1/x).le.eps) then
        feq = 2/x/3.d0*(1 + 2/5.d0/x + 8/35.d0/x/x)
        return
      end if
      if (abs(x).le.eps) then
        feq = 2 + (1 - x/2 + x*x/16)*dcmplx(log(abs(x/4)),pi*theta(s))
        return
      end if
      if (x.lt.0.d0) then
        sq = sqrt(1 - x)
        feq = 2 - sq*log((sq + 1)/(sq - 1))
        return
      else if (x.gt.1.d0) then
        sq = sqrt(x - 1)
        feq = 2 - 2*sq*atan(1/sq)
        return
      else if (x.eq.1.d0) then
        feq = 2
        return
      else
        sq = sqrt(1 - x)
        feq = 2 - sq*dcmplx(log((1 + sq)/(1 - sq)), -pi)
      end if
      return
      end

      double complex function f0(s,a1)
c     F function - one mass zero
      implicit double precision (a-h,o-z)
      common/spence/ber(9),pi,pi6,eps
      x = a1*a1/s
      if (abs(1/x).le.eps) then
        f0 = (1 + 1/3.d0/x + 1/6.d0/x/x)/2/x
        return
      end if
      if (abs(x).le.eps) then
        f0 = 1 + (1 - x + x*x/4)*dcmplx(log(abs(x)),pi*theta(s))
        return
      end if
      if ((x.lt.0.d0).or.(x.gt.1.d0)) then
        f0 = 1 + (x - 1)*log(1 - 1/x)
        return
      else if (x.eq.1.d0) then
        f0 = 1
        return
      else
        f0 = 1 + (x - 1)*dcmplx(log(1/x - 1),-pi)
      end if
      return
      end

      double complex function fneq(s,a1,a2)
c     F function - non equal masses
      implicit double precision (a-h,o-z)
      common/spence/ber(9),pi,pi6,eps
      external init_spence
      s1 = a1*a1
      s2 = a2*a2
      sum = s1 + s2
      dif = s1 - s2
      if (abs(s/dif).le.eps) then
         x = s/dif/dif
         f1 = sum - 4*s1*s2/dif*log(a1/a2)
         f2 = 3*sum*sum - 2*dif*dif - 12*s1*s2*sum/dif*log(a1/a2)
         f3 = sum*(5*sum*sum - 13/3.d0*dif*dif) 
     $        - (dif**4 - 6*(dif*sum)**2 + 5*sum**4)/dif*log(a1/a2)
         fneq = x/2*(f1 + x*f2/3 + x*x*f3/4)
         return
      end if
      fneq = 1 + (dif/s - sum/dif)*log(a2/a1)
      if (abs(dif/s).le.eps) then
        fneq = fneq - (1 - sum/s + sum*sum/s/s)
     $        *dcmplx(log(abs(s)/a1/a2),-pi*theta(s))
        return
      end if
      sum = (a1 + a2)**2
      dif = (a1 - a2)**2
      if (s.lt.dif) then
        sq1 = sqrt(sum - s)
        sq2 = sqrt(dif - s)
        fneq = fneq + sq1*sq2/s*log((sq1 + sq2)/(sq1 - sq2))
        return
      else if (s.lt.sum) then
        sq1 = sqrt(sum - s)
        sq2 = sqrt(s - dif)
        fneq = fneq - 2*sq1*sq2/s*atan(sq2/sq1)
        return
      else
        sq1 = sqrt(s - sum)
        sq2 = sqrt(s - dif)
        fneq = fneq
     1       - sq1*sq2/s*dcmplx(log((sq2 + sq1)/(sq2 - sq1)),-pi)
      end if
      return
      end

      double complex function b0(s,a1,a2)
c     B0 Veltman function
      implicit double precision (a-h,o-z)
      double complex f,db0
      logical dbstat,infstat
      common/spence/ber(9),pi,pi6,eps
      common/renorm/del,amiu2,infstat
      common/bdif/dbstat
      if (dbstat) then
        b0 = db0(s,a1,a2)
        return
      end if
      if (infstat) then
        b0 = del
        return
      end if
      if (s.eq.0.d0) then
        if (abs(a1) + abs(a2).eq.0.d0) then
          b0 = del
          return
        end if
        if (a1*a2.eq.0.d0) then
          b0 = del - log((a1 + a2)*(a1 + a2)/amiu2) + 1
        else
          b0 = del - log(a1*a2/amiu2)
          if (a1.ne.a2) then
            sum = a1*a1 + a2*a2
            dif = a1*a1 - a2*a2
            b0 = b0 + 1 - sum/dif*log(a1/a2)
          end if
        end if
        return
      end if
      if ((a1*a2).eq.0.d0) then
        if (a1.eq.a2) then
          b0 = del - log(abs(s)/amiu2) + 2
          if (s.gt.0.d0) b0 = b0 + (0.d0,1.d0)*pi
        else
          b0 = del - log((a1 + a2)*(a1 + a2)/amiu2) + 1 + f(s,a1,a2)
        end if
        return
      end if
      b0 = del - log(a1*a2/amiu2) + f(s,a1,a2)
      if (a1.ne.a2) then
        sum = a1*a1 + a2*a2
        dif = a1*a1 - a2*a2
        b0 = b0 + 1 - sum/dif*log(a1/a2)
      end if
      return
      end

      double complex function b1(s,a1,a2)
c     B1 Veltman function
      implicit double precision (a-h,o-z)
      double complex f,db1
      logical dbstat,infstat
      common/spence/ber(9),pi,pi6,eps
      common/renorm/del,amiu2,infstat
      common/bdif/dbstat
      if (dbstat) then
        b1 = db1(s,a1,a2)
        return
      end if
      if (infstat) then
        b1 = - del/2
        return
      end if
      if (s.eq.0.d0) then
        if (abs(a1) + abs(a2).eq.0.d0) stop 'b1 called for s=m1=m2=0'
        if (a1.eq.0.d0) then
          b1 = - del/2 + log(a2*a2/amiu2)/2 - 0.25d0
        else if (a2.eq.0.d0) then
          b1 = - del/2 + log(a1*a1/amiu2)/2 - 0.75d0
        else
          b1 = - del/2 + log(a1*a2/amiu2)/2
          if (a1.ne.a2) then
            dif = a1*a1 - a2*a2
            b1 = b1 - 0.75d0 - a2*a2/2/dif
     1         + (a1**4/dif/dif - 0.5d0)*log(a1/a2)
          end if
        end if
        return
      end if
      if ((a1*a2).eq.0.d0) then
        if (a1.eq.a2) then
          b1 = - (del - log(abs(s)/amiu2))/2 - 1
          if (s.gt.0.d0) b1 = b1 - (0.d0,1.d0)*pi/2
        else
          if (a1.eq.0.d0) then
            b1 = - (del - log(a2*a2/amiu2) + 1
     1         + (1 - a2*a2/s)*f(s,a1,a2))/2
          else
            b1 = - (del - log(a1*a1/amiu2) + 1
     1         + (1 + a1*a1/s)*f(s,a1,a2))/2
          end if
        end if
        return
      end if
      b1 = - (del - log(a1*a2/amiu2)
     1   + (s  + a1*a1 - a2*a2)/s*f(s,a1,a2))/2
      if (a1.ne.a2) then
        sum = a1*a1 + a2*a2
        dif = a1*a1 - a2*a2
        b1 = b1 - (1 - sum/dif*log(a1/a2))/2
      end if
      return
      end

      double complex function b22(s,a1,a2)
c     B22 Veltman function
      implicit double precision (a-h,o-z)
      double complex f,db22
      logical dbstat,infstat
      common/spence/ber(9),pi,pi6,eps
      common/renorm/del,amiu2,infstat
      common/bdif/dbstat
      if (dbstat) then
        b22 = db22(s,a1,a2)
        return
      end if
      if (infstat) then
        b22 = del*(3*a1*a1 + 3*a2*a2 - s)/12
        return
      end if
      if (s.eq.0.d0) then
        if (abs(a1) + abs(a2).eq.0.d0) then
          b22 = (0.d0,0.d0)
        else if (a1*a2.eq.0.d0) then
          as = a1 + a2
          b22 = as*as/4*(del - log(as*as/amiu2) + 1.5d0)
        else
          b22 = (a1*a1 + a2*a2)/4*(del - log(a1*a2/amiu2) + 1)
          if (a1.ne.a2) then
            sum = a1*a1 + a2*a2
            dif = a1*a1 - a2*a2
            b22 = b22 + sum/8 - (dif + 2*a1*a1*a2*a2/dif)*log(a1/a2)/4
          end if
        end if
        return
      end if
      if ((a1*a2).eq.0.d0) then
        if (a1.eq.a2) then
          b22 = - s/12*(del - log(abs(s)/amiu2)) - 2*s/9
          if (s.gt.0.d0) b22 = b22 - (0.d0,1.d0)*pi*s/12
        else
          as = a1 + a2
          b22 = ((3*as*as - s)*(5.d0/3 + del - log(as*as/amiu2))
     1        - (s - as*as)**2/s*f(s,a1,a2))/12
        end if
        return
      end if
      if (a1.eq.a2) then
        b22 = ((a1*a1 - s/6)*(del - log(a1*a1/amiu2)) - s/9 + a1*a1
     1      - (s - 4*a1*a1)/6*f(s,a1,a2))/2
      else
        sum = a1*a1 + a2*a2
        dif = a1*a1 - a2*a2
        b22 = ((3*sum - s)*(5.d0/3 + del - log(a1*a2/amiu2))
     1        + (s*sum - 8*a1*a1*a2*a2 - 3*dif*dif)/dif*log(a1/a2)
     2        - (s*s - 2*s*sum + dif**2)/s*f(s,a1,a2))/12
      end if
      return
      end

      double complex function x21(s,a1,a2)
c     X21 Veltman function - auxiliary combination
      implicit double precision (a-h,o-z)
      double complex b0
      logical infstat
      common/renorm/del,amiu2,infstat
      if (infstat) then
        x21 = del*s/3
        return
      end if
      sum = s + a2*a2 - a1*a1
      x21 = ((s + sum)*a(a1) - sum*a(a2))/3/s
     1    + (a2*a2 - sum*sum/s)*b0(s,a1,a2)/3 + (a1*a1 + a2*a2 - s/3)/6
      return
      end

      double complex function li2(x)
c     Complex dilogarithm
      implicit double precision (a-h,o-z)
      double complex x,y,z,sign
      common/spence/ber(9),pi,pi6,eps
      external init_spence
      if (x.eq.(0.d0,0.d0)) then
        li2 = (0.d0,0.d0)
        return
      else if (x.eq.(1.d0,0.d0)) then
        li2 = pi6
        return
      else if (x.eq.(-1.d0,0.d0)) then
        li2 = - pi6/2
        return
      else if (x.eq.(0.5d0,0.d0)) then
        li2 = pi6/2 - log(2.d0)**2/2
        return
      end if
      sign = (1.d0,0.d0)
      y = x
      li2 = (0.d0,0.d0)
      if ((abs(y).gt.1)) then
        li2 = - pi6 - log(-y)**2/2
        y = 1/y
        sign = - sign
      end if
      if (dble(y).gt.0.5d0) then
        li2 = li2 + sign*(pi6 - log(y)*log(1 - y))
        y = 1 - y
        sign = - sign
      end if
      y = - log(1 - y)
      sign = sign*y
      li2 = li2 + sign*(1 - y/4)
      y = y*y
      z = (1.d0,0.d0)
      den = 1
      do 10 i=1,9
        z = z*y
c     avoid underflow
        if (abs(z).lt.1.d-80) z=(0.d0,0.d0)
        den = 2*i*(2*i + 1)*den
10      li2 = li2 + sign*ber(i)*z/den
      return
      end

      double precision function theta(x)
      implicit double precision (a-h,o-z)
      if (x.lt.0.d0) then
         theta = 0
      else
         theta = 1
      end if
      return
      end

      block data init_spence
      implicit double precision (a-h,o-z)
      logical dbstat,infstat
      double complex c0_0
      common/cmem/c0_0,p2_0,q2_0,pq_0,a1_0,a2_0,a3_0,eps_c0
      common/spence/ber(9),pi,pi6,eps
      common/renorm/del,amiu2,infstat
      common/bdif/dbstat
c     Bernoulli numbers
      data ber/1.6666666666667d-1,-3.33333333333333d-2,
     1         2.3809538095381d-2,-3.33333333333333d-2,
     2         7.5757575757576d-2,-2.53113553113553d-1,
     3         1.1666666666667d0 ,-7.09215686274510d0,
     4         54.9711779448622d0/
c     pi6 = pi^2/6
c     eps: variable in loop integrals denominators: p^2 - m^2 - i*eps 
      data pi,pi6,eps/3.1415926535897932d0,1.6449340668482264d0,1.d-5/
      data dbstat/.false./
      data del,amiu2,infstat/0.d0,8.3152513d3,.false./
      data c0_0/(0.d0,0.d0)/
      data p2_0,q2_0,pq_0,a1_0,a2_0,a3_0/6*-1.d-5/
      data eps_c0/1.d-8/
      end

      subroutine db_stat(dbstat)
      logical dbstat,dbst
      common/bdif/dbst
      dbst = dbstat
      return
      end

      subroutine inf_stat(istat)
      implicit double precision (a-h,o-z)
      logical istat,infstat
      common/renorm/del,amiu2,infstat
      infstat = istat
      return
      end



