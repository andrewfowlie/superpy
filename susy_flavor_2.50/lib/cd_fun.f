c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: CD_FUN.F
c     Released: 15:11:1995 (J.R.)
c     Corrected: 20.01.2000 (J.R.)
c     Behaviour around poles in denominators regularized 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains three- and four-point loop integrals      c
c     at zero external momentum. The parameter amiu2 in the        c
c     common/renorm/ helps to check the cancelation of the         c
c     divergences by counterterms. The amiu2 parameter in the      c
c     common/renorm/ is the 't Hooft mass unit which also should   c
c     cancel in the final expressions.                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function cp0(am,bm,cm)
c     Three point scalar function (1 in the numerator)
      implicit double precision (a-h,o-z)
      dimension x(3),ind(3)
      x(1) = am*am
      x(2) = bm*bm
      x(3) = cm*cm
      call sort_arr(x,ind,3)
      a = x(1)
      b = x(2)
      c = x(3)
      if (a.eq.0) then
        if (b.eq.0) stop '2 or more masses=0 in cp0'
        if (b.eq.c) then
          cp0 = - 1/b
        else 
          cp0 = log(b/c)/(c - b)
        end if
        return
      end if
      if (a.eq.c) then
        cp0 = - 1.d0/2/a
      else if (a.eq.b) then
        cp0 =  (1 - log(c/a)/(1 - a/c))/(c - a)
      else if (b.eq.c) then
        cp0 =  (1 - log(a/c)/(1 - c/a))/(a - c)
      else
        cp0 = - (log(b/a)/(1 - a/b) - log(c/a)/(1 - a/c))/(b - c)
      end if
      return
      end

      double precision function cp1(am,bm,cm)
c     Three point scalar function (k^2 in the numerator)
      implicit double precision (a-h,o-z)
      logical infstat
      dimension x(3),ind(3)
      common/renorm/del,amiu2,infstat
      external init_spence
      x(1) = am*am
      x(2) = bm*bm
      x(3) = cm*cm
      call sort_arr(x,ind,3)
      a = x(1)
      b = x(2)
      c = x(3)
      if (a.eq.0) then
        if (c.eq.0) stop '3x0 masses in cp1'
        cp1 = del + log(amiu2/c)
        if (b.eq.0) then
          cp1 = cp1 + 1
          return
        end if
        if (b.ne.c) cp1 = cp1 + 1 + log(c/b)/(1 - c/b)
        return
      end if
      if (a.eq.c) then
        cp1 = del + log(amiu2/c) - 0.5d0
      else if (a.eq.b) then
        cp1 = del + log(amiu2/a) + (1 - log(c/a)/(1 - a/c))/(1 - a/c)
      else if (b.eq.c) then
        cp1 = del + log(amiu2/c) + (1 - log(a/c)/(1 - c/a))/(1 - c/a)
      else
        cp1 = del + log(amiu2/a) + 1 - log(b/a)/(a/b - 1)/(c/b - 1)
     1      - log(c/a)/(a/c - 1)/(b/c - 1)
      end if
      return
      end

      double precision function dp0(am,bm,cm,dm)
c     Four point scalar function (1 in the numerator)
      implicit double precision (a-h,o-z)
      dimension x(4),ind(4)
      x(1) = am*am
      x(2) = bm*bm
      x(3) = cm*cm
      x(4) = dm*dm
      call sort_arr(x,ind,4)
      a = x(1)
      b = x(2)
      c = x(3)
      d = x(4)
      if (a.eq.0) then
c     Mass combination 00cd
        if (b.eq.0) stop '2 or more masses=0 in dp0'
        if (b.eq.d) then
c     Mass combination 0bbb
          dp0 = 1.d0/2/b/b
        else if (b.eq.c) then
c     Mass combination 0bbd
          dp0 = (1/b - log(d/b)/(d - b))/(d - b)
        else if (c.eq.d) then
c     Mass combination 0bcc
          dp0 = (1/c - log(b/c)/(b - c))/(b - c)
        else
c     Mass combination 0bcd
          dp0 = (log(c/b)/(c - b) - log(d/b)/(d - b))/(d - c)
        end if
      return
      end if
      if (a.eq.d) then
c     Mass combination aaaa
        dp0 = 1.d0/6/a/a
      else if (a.eq.c) then
c     Mass combination aaad
        dp0 = ((1 + d/a)/2 - log(d/a)/(1 - a/d))/(d - a)/(d - a)
      else if (a.eq.b) then
        if (c.eq.d) then
c     Mass combination aacc
          dp0 = - (2 + (a + c)/(a - c)*log(c/a))/(c - a)/(c - a)
        else
c     Mass combination aacd
          dp0 = - 1/(d - a)/(c - a) - log(c/a)/(1 - d/c)/(c - a)/(c - a)
     1                              - log(d/a)/(1 - c/d)/(d - a)/(d - a)
        end if
      else
        if (b.eq.d) then
c     Mass combination abbb
          dp0 = ((1 + a/b)/2 - log(a/b)/(1 - b/a))/(a - b)/(a - b)
        else if (b.eq.c) then
c     Mass combination abbd
          dp0 = - 1/(d - b)/(a - b) - log(a/b)/(1 - d/a)/(a - b)/(a - b)
     1                              - log(d/b)/(1 - a/d)/(d - b)/(d - b)
        else if (c.eq.d) then
c     Mass combination abcc
          dp0 = - 1/(b - c)/(a - c) - log(a/c)/(1 - b/a)/(a - c)/(a - c)
     1                              - log(b/c)/(1 - a/b)/(b - c)/(b - c)
        else 
c     Mass combination abcd
          dp0 = - log(b/a)/(1 - a/b)/(b - c)/(b - d)
     1          - log(c/a)/(1 - a/c)/(c - b)/(c - d)
     2          - log(d/a)/(1 - a/d)/(d - b)/(d - c)
        end if
      end if
      return
      end

      double precision function dp1(am,bm,cm,dm)
c     Four point scalar function (k^2 in the numerator)
      implicit double precision (a-h,o-z)
      dimension x(4),ind(4)
      x(1) = am*am
      x(2) = bm*bm
      x(3) = cm*cm
      x(4) = dm*dm
      call sort_arr(x,ind,4)
      a = x(1)
      b = x(2)
      c = x(3)
      d = x(4)
      if (a.eq.0) then
c     Mass combination 000d
        if (c.eq.0) stop '3 or more masses=0 in dp1'
        if (b.eq.0) then
          if (c.eq.d) then
c     Mass combination 00cc
            dp1 = - 1/c
          else
c     Mass combination 00cd
            dp1 = - log(d/c)/(d - c)
          end if
        else if (b.eq.d) then
c     Mass combination 0bbb
          dp1 = - 1.d0/2/b
        else if (b.eq.c) then
c     Mass combination 0bbd
          dp1 = (1 - log(d/b)/(1 - b/d))/(d - b)
        else if (c.eq.d) then
c     Mass combination 0bcc
          dp1 = (1 - log(b/c)/(1 - c/b))/(b - c)
        else
c     Mass combination 0bcd
          dp1 = (log(b/d)/(1 - d/b) - log(c/d)/(1 - d/c))/(c - b)
        end if
      return
      end if
      if (a.eq.d) then
c     Mass combination aaaa
        dp1 = - 1.d0/3/a
      else if (a.eq.c) then
c     Mass combination aaad
        dp1 = ((3*d - a)/2 - d*log(d/a)/(1 - a/d))/(d - a)/(d - a)
      else if (a.eq.b) then
        if (c.eq.d) then
c     Mass combination aacc
          dp1 = - (a + c + 2*c*log(c/a)/(1 - c/a))/(c - a)/(c - a)
        else
c     Mass combination aacd
          dp1 = - a/(d - a)/(c - a)
     1          - c*log(c/a)/(1 - d/c)/(c - a)/(c - a)
     2          - d*log(d/a)/(1 - c/d)/(d - a)/(d - a)
        end if
      else
        if (b.eq.d) then
c     Mass combination abbb
          dp1 = ((3*a - b)/2 - a*log(a/b)/(1 - b/a))/(a - b)/(a - b)
        else if (b.eq.c) then
c     Mass combination abbd
          dp1 = - b/(d - b)/(a - b)
     1          - a*log(a/b)/(1 - d/a)/(a - b)/(a - b)
     2          - d*log(d/b)/(1 - a/d)/(d - b)/(d - b)
        else if (c.eq.d) then
c     Mass combination abcc
          dp1 = - c/(b - c)/(a - c)
     1          - a*log(a/c)/(1 - b/a)/(a - c)/(a - c)
     2          - b*log(b/c)/(1 - a/b)/(b - c)/(b - c)
        else
c     Mass combination abcd
          dp1 = - b*log(b/a)/(1 - a/b)/(b - c)/(b - d)
     1          - c*log(c/a)/(1 - a/c)/(c - b)/(c - d)
     2          - d*log(d/a)/(1 - a/d)/(d - b)/(d - c)
        end if
      end if
      return
      end

      subroutine sort_arr(x,ind,n)
      implicit double precision (a-h,o-z)
      dimension x(n),ind(n)
      common/cp_acc/eps
      external init_cd_fun
      if (n.lt.2) return
c     initial ordering
      do i=1,n
         ind(i) = i
      end do
c     Sort array. Final ordering x(1) <= x(2) <= ... <= x(n)
      do i=1,n-1
        do j=1,n-i
          if (x(j).gt.x(j+1)) then
             tmp    = x(j+1)
             x(j+1) = x(j)
             x(j)   = tmp
             itmp     = ind(j+1)
             ind(j+1) = ind(j)
             ind(j)   = itmp
           end if
        end do
      end do
c     Replace very close masses by the equal ones, i.e.
c     x(i) = x(i+1) = ... = x(i+n) = (x(1) + ... x(n))/n
      i = 1
 10   if (abs(x(i+1) - x(i))/max(abs(x(i+1)),abs(x(i))).le.eps) then
         avg = x(i)
         j = 1
 20      avg = (j*avg + x(i+j))/(j+1)
         j = j + 1 
         if ((i+j).gt.n) goto 30 
         if (abs(x(i+j) - avg)/max(abs(x(i+j)),abs(avg)).le.eps) goto 20
 30      j = j - 1
         do k=0,j
            x(i+k) = avg
         end do
         i = i + j
      else
         i = i + 1
      end if      
      if (i.ge.n) return
      goto 10 
      end

      block data init_cd_fun
      implicit double precision (a-h,o-z)
      common/cp_acc/eps
      data eps/5.d-3/
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Three point scalar functions for b->s gamma decay           c
c     cp_ij(m1,m1,m2) = - c_ij(0,0,m1,m1,m2), ij=11,12,21,22,23   c
c     - c0(p,q-p,m1,m2,m3)  = cp0_1 + pq*cp0_2 + q^2*cp0_3        c
c                                   + p^2*cp0_4                   c
c     - c24(p,q-p,m1,m2,m3) = cp24_1 + pq*cp24_2 + q^2*cp24_3     c
c                                    + p^2*cp24_4                 c
c      Identities:                                                c
c      cp0_3 = - cp0_2                                            c
c      cp24_2 = - cp24_3 = cp23                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function cp0_1(am,bm)
      implicit double precision (a-h,o-z)
      common/cp_acc/eps
      a = am*am
      b = bm*bm
      if (abs(a - b)/max(abs(a),abs(b)).lt.eps) then
         avg = (a + b)/2
         a = avg
         b = avg
      end if
      if (a.eq.0) then
        stop 'first mass=0 in cp0_1'
      else if (b.eq.0) then
        cp0_1 = - 1/a
      else if (a.eq.b) then
        cp0_1 = - 0.5d0/a
      else
        d = a - b
        cp0_1 = - 1/d - b/d/d*log(b/a)
      end if
      return
      end

      double precision function cp0_2(am,bm)
      implicit double precision (a-h,o-z)
      common/cp_acc/eps
      a = am*am
      b = bm*bm
      if (abs(a - b)/max(abs(a),abs(b)).lt.eps) then
         avg = (a + b)/2
         a = avg
         b = avg
      end if
      if (a.eq.0) then
        stop 'first mass=0 in cp0_2'
      else if (b.eq.0) then
        cp0_2 = 0.5d0/a/a
      else if (a.eq.b) then
        cp0_2 = 1/12.d0/a/a
      else
        d = a - b
        cp0_2 = (a + 5*b)/2/d**3 + b*(2*a + b)/d**4*log(b/a)
      end if
      return
      end

      double precision function cp0_4(am,bm)
      implicit double precision (a-h,o-z)
      common/cp_acc/eps
      a = am*am
      b = bm*bm
      if (abs(a - b)/max(abs(a),abs(b)).lt.eps) then
         avg = (a + b)/2
         a = avg
         b = avg
      end if
      if (a.eq.0) then
        stop 'first mass=0 in cp0_4'
      else if (b.eq.0) then
        cp0_4 = - 1/3.d0/a/a
      else if (a.eq.b) then
        cp0_4 = - 1/12.d0/a/a
      else
        d = a - b
        cp0_4 = - (2*a*a + 5*a*b - b*b)/6/a/d**3 - a*b/d**4*log(b/a)
      end if
      return
      end

      double precision function cp11(am,bm)
      implicit double precision (a-h,o-z)
      common/cp_acc/eps
      a = am*am
      b = bm*bm
      if (abs(a - b)/max(abs(a),abs(b)).lt.eps) then
         avg = (a + b)/2
         a = avg
         b = avg
      end if
      if (a.eq.0) then
        stop 'first mass=0 in cp11'
      else if (b.eq.0) then
        cp11 = - 1/4.d0/a
      else if (a.eq.b) then
        cp11 = - 1/6.d0/a
      else
        d = a - b
        cp11 = - (a - 3*b)/4/d/d + b*b/2/d**3*log(b/a)
      end if
      return
      end

      double precision function cp12(am,bm)
      implicit double precision (a-h,o-z)
      common/cp_acc/eps
      a = am*am
      b = bm*bm
      if (abs(a - b)/max(abs(a),abs(b)).lt.eps) then
         avg = (a + b)/2
         a = avg
         b = avg
      end if
      if (a.eq.0) then
        if (b.eq.0) stop '2 masses=0 in cp12'
        cp12 = - 0.5d0/b
      else if (b.eq.0) then
        cp12 = - 0.5d0/a
      else if (a.eq.b) then
        cp12 = - 1/6.d0/a
      else
        d = a - b
        cp12 = - (a + b)/2/d/d - a*b/d**3*log(b/a)
      end if
      return
      end

      double precision function cp21(am,bm)
      implicit double precision (a-h,o-z)
      common/cp_acc/eps
      a = am*am
      b = bm*bm
      if (abs(a - b)/max(abs(a),abs(b)).lt.eps) then
         avg = (a + b)/2
         a = avg
         b = avg
      end if
      if (a.eq.0) then
        stop 'first mass=0 in cp21'
      else if (b.eq.0) then
        cp21 = - 1/9.d0/a
      else if (a.eq.b) then
        cp21 = - 1/12.d0/a
      else
        d = a - b
        cp21 = - (2*a*a - 7*a*b + 11*b*b)/18/d**3 - b**3/3/d**4*log(b/a)
      end if
      return
      end

      double precision function cp22(am,bm)
c     CAUTION: cp22 = am*am*cp0_4
      implicit double precision (a-h,o-z)
      common/cp_acc/eps
      a = am*am
      b = bm*bm
      if (abs(a - b)/max(abs(a),abs(b)).lt.eps) then
         avg = (a + b)/2
         a = avg
         b = avg
      end if
      if (a.eq.0) then
        if (b.eq.0) stop '2 masses=0 in cp22'
        cp22 = - 1/6.d0/b
      else if (b.eq.0) then
        cp22 = - 1/3.d0/a
      else if (a.eq.b) then
        cp22 = - 1/12.d0/a
      else
        d = a - b
        cp22 = - (2*a*a + 5*a*b - b*b)/6/d**3 - a*a*b/d**4*log(b/a)
      end if
      return
      end

      double precision function cp23(am,bm)
      implicit double precision (a-h,o-z)
      common/cp_acc/eps
      a = am*am
      b = bm*bm
      if (abs(a - b)/max(abs(a),abs(b)).lt.eps) then
         avg = (a + b)/2
         a = avg
         b = avg
      end if
      if (a.eq.0) then
        if (b.eq.0) stop '2 masses=0 in cp23'
        cp23 = - 1/6.d0/b
      else if (b.eq.0) then
        cp23 = - 1/12.d0/a
      else if (a.eq.b) then
        cp23 = - 1/24.d0/a
      else
        d = a - b
        cp23 = - (a*a - 5*a*b - 2*b*b)/12/d**3 + a*b*b/2/d**4*log(b/a)
      end if
      return
      end

      double precision function cp24_4(am,bm)
      implicit double precision (a-h,o-z)
      common/cp_acc/eps
      a = am*am
      b = bm*bm
      if (abs(a - b)/max(abs(a),abs(b)).lt.eps) then
         avg = (a + b)/2
         a = avg
         b = avg
      end if
      if (a.eq.0) then
        stop 'first mass=0 in cp24_4'
      else if (b.eq.0) then
        cp24_4 = 5/72.d0/a
      else if (a.eq.b) then
        cp24_4 =  1/24.d0/a
      else
        d = a - b
        cp24_4 = (5*a*a - 22*a*b + 5*b*b)/72/d**3
     1         - b*b*(3*a - b)/12/d**4*log(b/a)
      end if
      return
      end




