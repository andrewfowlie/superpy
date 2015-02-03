c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: C_FUN.F
c     Revised: 29: 9:1992 (P.Ch.)
c     Corrected code for s3 for a=b=0
c     Revised: 23: 6:1999 (J.R.)
c     Preventing overflows in c0 for alfa=0,1
c     Revised: 06: 6:2013 (J.R.)
c     Added moment to list of C-functions arguments; simplified
c     resetting after argument change - now automatic and inside c0 only
c     Revised: 17.11.2013 (J.R.)
c     New code with for c0 written to prevent loss of accuracy due to
c     strong cancellations. Error in c0 phase corrected (coming from
c     misprint in original PV paper!)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains three-point standard loop integrals         c
c     c0(p2,q2,pq;m1,m2,m3), c11, c12, c21, c22, c23, c24.           c
c     For definitions of these functions see:                        c
c     A.Axelrod@Nucl.Phys.B209(1982)p.349 (compare: Chankowski,      c
c     Pokorski,Rosiek, Nucl.Phys.B423(1994)p.437 available as        c
c     hep-ph/9303309.                                                c
c     The convention for external momenta and internal masses are:   c
c                                                                    c
c                       ________ p (outgoing)                        c
c                   / |                                              c
c            m1   /   |                                              c
c               /     |                                              c
c             /       |                                              c
c     -------         | m2              = c0(p2,q2,pq;m1,m2,m3)      c
c             \       |                                              c
c               \     |                                              c
c            m3   \   |                                              c
c                   \ |_________ q (outgoing)                        c
c                                                                    c
c                                                                    c
c                                                                    c
c     -i/(4 pi)^2 c0 = INT d^4k/(2 pi)^4                             c
c            1/(k^2-m2^2)/((k-q)^2-m3^2)/((k+p)^2-m1^2)              c
c                    = INT d^4k/(2 pi)^4                             c
c            1/(k^2-m1^2)/((k+p)^2-m2^2)/((k+p+q)^2-m3^2)            c
c                                                                    c
c   p2 = p^2    q2 = q^2    pq = 1/2 ((p + q)^2 - p^2 - q^2)         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function c11(p2,q2,pq,a1,a2,a3)
      implicit double precision (a-h,o-z)
      double complex b0,c0
      logical infstat
      common/renorm/del,amiu2,infstat
      if (infstat) then
        c11 = (0.d0,0.d0)
        return
      end if
c     check momenta and mass scales. 
      det = p2*q2 - pq*pq
      scp = abs(det)**0.25d0
      scm = max(abs(a1),abs(a2),abs(a3))
c     call full or expanded version of c-fun code
      if ((p2.le.(a1+a2)**2).and.(q2.le.(a2+a3)**2)
     $     .and.((p2+2*pq+q2).le.(a1+a3)**2)
     $     .and.((scp/scm).le.0.2d0)) then 
         c11 = c11exp(p2,q2,pq,a1,a2,a3)
      else
         s1 = p2 + a1*a1 - a2*a2
         s2 = q2 + 2*pq + a2*a2 - a3*a3
         x = p2 + q2 + 2*pq
         c11 = q2*(b0(x,a1,a3) - b0(q2,a2,a3) 
     $        + s1*c0(p2,q2,pq,a1,a2,a3))
     $        - pq*(b0(p2,a1,a2) - b0(x,a1,a3) 
     $        + s2*c0(p2,q2,pq,a1,a2,a3))
         c11 = - c11/det/2
      end if
      return
      end

      double complex function c12(p2,q2,pq,a1,a2,a3)
      implicit double precision (a-h,o-z)
      double complex b0,c0
      logical infstat
      common/renorm/del,amiu2,infstat
      if (infstat) then
        c12 = (0.d0,0.d0)
        return
      end if
c     check momenta and mass scales. 
      det = p2*q2 - pq*pq
      scp = abs(det)**0.25d0
      scm = max(abs(a1),abs(a2),abs(a3))
c     call full or expanded version of c-fun code
      if ((p2.le.(a1+a2)**2).and.(q2.le.(a2+a3)**2)
     $     .and.((p2+2*pq+q2).le.(a1+a3)**2)
     $     .and.((scp/scm).le.0.2d0)) then 
         c12 = c12exp(p2,q2,pq,a1,a2,a3)
      else
         s1 = p2 + a1*a1 - a2*a2
         s2 = q2 + 2*pq + a2*a2 - a3*a3
         x = p2 + q2 + 2*pq
         c12 = - pq*(b0(x,a1,a3) - b0(q2,a2,a3) 
     $        + s1*c0(p2,q2,pq,a1,a2,a3))
     $        + p2*(b0(p2,a1,a2) - b0(x,a1,a3) 
     $        + s2*c0(p2,q2,pq,a1,a2,a3))
         c12 = - c12/det/2
      end if
      return
      end

      double complex function c21(p2,q2,pq,a1,a2,a3)
      implicit double precision (a-h,o-z)
      double complex b0,b1,c11,c24
      logical infstat
      common/renorm/del,amiu2,infstat
      if (infstat) then
        c21 = (0.d0,0.d0)
        return
      end if
      det = p2*q2 - pq*pq
      s1 = p2 + a1*a1 - a2*a2
      s2 = q2 + 2*pq + a2*a2 - a3*a3
      x = p2 + q2 + 2*pq
      c21 = q2*(b1(x,a1,a3) + b0(q2,a2,a3) + s1*c11(p2,q2,pq,a1,a2,a3)
     $     + 2*c24(p2,q2,pq,a1,a2,a3)) - pq*(b1(p2,a1,a2) - b1(x,a1,a3)
     $     + s2*c11(p2,q2,pq,a1,a2,a3))
      c21 = - c21/det/2
      return
      end

      double complex function c22(p2,q2,pq,a1,a2,a3)
      implicit double precision (a-h,o-z)
      double complex b1,c12,c24
      logical infstat
      common/renorm/del,amiu2,infstat
      if (infstat) then
        c22 = (0.d0,0.d0)
        return
      end if
      det = p2*q2 - pq*pq
      s1 = p2 + a1*a1 - a2*a2
      s2 = q2 + 2*pq + a2*a2 - a3*a3
      x = p2 + q2 + 2*pq
      c22 = - p2*(b1(x,a1,a3) - s2*c12(p2,q2,pq,a1,a2,a3) 
     $     - 2*c24(p2,q2,pq,a1,a2,a3)) - pq*(b1(x,a1,a3) 
     $     - b1(q2,a2,a3) + s1*c12(p2,q2,pq,a1,a2,a3)) 
      c22 = - c22/det/2
      return
      end

      double complex function c23(p2,q2,pq,a1,a2,a3)
      implicit double precision (a-h,o-z)
      double complex b0,b1,c11,c24
      logical infstat
      common/renorm/del,amiu2,infstat
      if (infstat) then
        c23 = (0.d0,0.d0)
        return
      end if
      det = p2*q2 - pq*pq
      s1 = p2 + a1*a1 - a2*a2
      s2 = q2 + 2*pq + a2*a2 - a3*a3
      x = p2 + q2 + 2*pq
      c23 = - pq*(b1(x,a1,a3) + b0(q2,a2,a3) + s1*c11(p2,q2,pq,a1,a2,a3)
     $    + 2*c24(p2,q2,pq,a1,a2,a3)) + p2*(b1(p2,a1,a2) - b1(x,a1,a3)
     $    + s2*c11(p2,q2,pq,a1,a2,a3))
      c23 = - c23/det/2
      return
      end

      double complex function c24(p2,q2,pq,a1,a2,a3)
      implicit double precision (a-h,o-z)
      double complex b0,c0,c11,c12
      logical infstat
      common/renorm/del,amiu2,infstat
      if (infstat) then
        c24 = - del/4
        return
      end if
      s1 = p2 + a1*a1 - a2*a2
      s2 = q2 + 2*pq + a2*a2 - a3*a3
      c24 = (s1*c11(p2,q2,pq,a1,a2,a3) + s2*c12(p2,q2,pq,a1,a2,a3) 
     $     + 2*a1*a1*c0(p2,q2,pq,a1,a2,a3) - b0(q2,a2,a3) - 1)/4
      return
      end

      double complex function c0(p2,q2,pq,a1,a2,a3)
      implicit double precision (a-h,m,o-z)
      logical infstat
      double complex c0_0,c0_body
      common/cmem/c0_0,p2_0,q2_0,pq_0,a1_0,a2_0,a3_0,eps_c0
      common/renorm/del,amiu2,infstat
      if (infstat) then
        c0 = (0.d0,0.d0)
        return
      end if
c     if arguments changed, calculate new c0 value; use old one otherwise 
      if (max(abs(p2-p2_0),abs(q2-q2_0),abs(pq-pq_0),abs(a1-a1_0),
     $     abs(a2-a2_0),abs(a3-a3_0)).ge.eps_c0) then
c     rescale momenta to maximal one to avoid cancellations
         r = max(abs(p2)**0.5d0,abs(q2)**0.5d0,abs(pq)**0.5d0,
     $        abs(a1),abs(a2),abs(a3))
         c0 = c0_body(p2/r/r,q2/r/r,pq/r/r,a1/r,a2/r,a3/r)/r/r
         p2_0 = p2
         q2_0 = q2
         pq_0 = pq
         a1_0 = a1
         a2_0 = a2
         a3_0 = a3
         c0_0 = c0
      else
         c0 = c0_0
      end if
      return
      end

      double complex function c0_body(p2,q2,pq,a1,a2,a3)
      implicit double precision (a-h,m,o-z)
      double complex s3,alfa,y0,y1,y2,y3,sf1,sf2,sf3,f
      common/spence/ber(9),pi,pi6,eps
      a = q2
      b = p2
      c = 2*pq
      d = a2*a2 - a3*a3 - q2
      e = a1*a1 - a2*a2 - p2 - 2*pq
      f = dcmplx(a3*a3,-eps)
      if (b.eq.0.d0) then
        alfa = - dcmplx(a/c,0.d0)
      else
        delta = c*c - 4*a*b
        if (delta.gt.0.d0) then
          alfa = dcmplx( - c - sign(sqrt(delta),c),0.d0)/2/b
        else
          alfa = dcmplx( - c, - sqrt( - delta))/2/b
        end if
      end if
      y0 = - (d + e*alfa)/(c + 2*b*alfa)
      y1 = y0 + alfa
      sf1 = s3(y1,b,c + e,a + d + f)
      if (abs(1 - alfa).ne.0.d0) then
         y2 = y0/(1 - alfa)
         sf2 = s3(y2,a + b + c,e + d,f)
      else 
         sf2 = (0.d0,0.d0)
      end if
      if (abs(alfa).ne.0.d0) then
         y3 = - y0/alfa
         sf3 = s3(y3,a,d,f)
      else 
         sf3 = (0.d0,0.d0)
      end if
      c0_body = (sf1 - sf2 + sf3)/(c + 2*b*alfa)
c     check if phase is just due to numerical error
      if ((p2.le.(a1+a2)**2).and.(q2.le.(a2+a3)**2)
     $     .and.((p2+2*pq+q2).le.(a1+a3)**2)) c0_body = dble(c0_body)
      return
      end

      double complex function s3(y0,a,b,c)
      implicit double precision (a-h,o-z)
      double complex c,eta,R,y0,y1,y2,z1,z2,delta
      if (a.eq.0.d0) then
        if (b.eq.0.d0) then
          s3 = (0.d0,0.d0)
          return
        end if
        z1 = b
        z2 = c/b
        y1 = - c/b
        s3 = R(y0,y1) + log(1 - 1/y0)*(eta(z1,z2) - eta(z1,y0 - y1))
        return
      end if
      delta = b*b - 4*a*c
      y1 = ( - b - sqrt(delta))/2/a
      y2 = ( - b + sqrt(delta))/2/a
      z1 = a
      z2 = c/a
      s3 = R(y0,y1) + R(y0,y2) - log(1 - 1/y0)*(eta(y0 - y1,y0 - y2)
     $   - eta( - y1, - y2) + eta(z1,z2) - eta(z1,y0*(y0 + b/a) + z2))
      return
      end

      double complex function R(z0,z1)
      implicit double precision (a-h,o-z)
      double complex eta,li2,z0,z1
      R = li2(z0/(z0 - z1)) - li2((z0 - 1)/(z0 - z1))
     $     + eta( - z1,1/(z0 - z1))*log(z0/(z0 - z1))
     $     - eta(1 - z1,1/(z0 - z1))*log((z0 - 1)/(z0 - z1))
      return
      end

      double complex function eta(z1,z2)
      implicit double precision (a-h,o-z)
      double complex z1,z2
      common/spence/ber(9),pi,pi6,eps
      zi1 = dimag(z1)
      zi2 = dimag(z2)
      zi12 = dimag(z1*z2)
      eta = dcmplx(0.d0,2*pi)*(theta( - zi1)*theta( - zi2)*theta(zi12)
     $     - theta(zi1)*theta(zi2)*theta( - zi12))
      return
      end

      double complex function c0_on(s,am)
c     c0 for p=q=0, pq=s/2 and equal masses in loop
      implicit double precision (a-h,m,o-z)
      logical infstat
      common/spence/ber(9),pi,pi6,eps
      common/renorm/del,amiu2,infstat
      if (infstat) then
        c0_on = (0.d0,0.d0)
        return
      end if
      if (s.eq.0) then
        c0_on = 0.5d0/am/am
        return
      end if
      x = 4*am*am/s
      if (x.le.0.d0) then
        stop 'c0_on undefined for s<0 or m=0'
      else if (x.lt.1.d0) then
        sq = sqrt(1 - x)
        c0_on = - dcmplx(log((1 + sq)/(1 - sq)),-pi)**2/2/s
        return
      else
        c0_on = 2/s*asin(sqrt(1/x))**2
      end if
      return
      end




