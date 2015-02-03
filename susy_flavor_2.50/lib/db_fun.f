c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: DB_FUN.F
c     Released: 28: 3:1994 (J.R.)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains derivatives of the two-point standard     c
c     loop integrals b0(s,m1,m2), b1(s,m1,m2), b22(s,m1,m2)        c
c     These functions are defined as in: A.Axelrod,                c
c     Nucl.Phys.B209(1982)p.349 except that they have              c
c     opposite signs (compare: Chankowski,Pokorski,Rosiek,         c
c     Nucl.Phys.B423(1994)p.437, available as hep-ph/9303309.      c
c     The parameter del in the common/renorm/ helps to check the   c
c     cancelation of the divergences by counterterms. The amiu2    c
c     parameter in the common/renorm/  is the 't Hooft mass        c
c     unit which also should cancel in the final expressions       c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function df(s,a1,a2)
c     F function derivative - general case
      implicit double precision (a-h,o-z)
      double complex dfeq,df0,dfneq
      if (s.eq.(a1 + a2)**2) stop 'DF called for s=(m1 + m2)^2'
      if (a1.eq.a2) then
        df = dfeq(s,a1)
        return
      else if ((a1*a2).eq.0.d0) then
        df = df0(s,a1 + a2)
        return
      else
        df = dfneq(s,a1,a2)
      end if
      return
      end

      double complex function dfeq(s,a1)
c     F function derivative - equal masses
      implicit double precision (a-h,o-z)
      double complex feq
      common/spence/ber(9),pi,pi6,eps
      if (a1.eq.0) then
        dfeq = - 1/s
        return
      end if
      if (abs(s/a1/a1).le.eps) then
        dfeq = (1 + s/5/a1/a1)/6/a1/a1
        return
      end if
      dfeq = (2*a1*a1/s*feq(s,a1) - 1)/(s - 4*a1*a1)
      return
      end

      double complex function df0(s,a1)
c     F function derivative - one mass zero
      implicit double precision (a-h,o-z)
      double complex f0
      common/spence/ber(9),pi,pi6,eps
      if (abs(s/a1/a1).le.eps) then
        df0 = (0.5d0 + s/3/a1/a1)/a1/a1
        return
      end if
      df0 = (a1*a1/s*f0(s,a1) - 1)/(s - a1*a1)
      return
      end

      double complex function dfneq(s,a1,a2)
c     F function derivative - non equal masses
      implicit double precision (a-h,o-z)
      double complex fneq
      common/spence/ber(9),pi,pi6,eps
      sum  = a1*a1 + a2*a2
      dif  = a1*a1 - a2*a2
      prod = a1*a1*a2*a2
      al   = log(a1/a2)
      if (abs(s/sum).le.eps) then
        f1 = (sum*dif - 4*prod*al)/2
        f2 = (3*sum*sum - 2*dif*dif - 12*sum*prod*al/dif)/3/dif
        dfneq = (f1 + s*f2)/dif**3
        return
      end if
      if (s.eq.(a1 - a2)**2) then
        dfneq = ((a1 + a2)/(a1 - a2)*al - 2)/s
        return
      end if
      dfneq = ((sum - dif*dif/s)*fneq(s,a1,a2) + sum - s
     1      - 4*prod/dif*al)/(s*s - 2*s*sum + dif*dif)
      return
      end

      double complex function db0(s,a1,a2)
c     B0 Veltman function derivative
      implicit double precision (a-h,o-z)
      double complex df
      logical infstat
      common/renorm/del,amiu2,infstat
      if (infstat) then
        db0 = 0
        return
      end if
      db0 = df(s,a1,a2)
      return
      end

      double complex function db1(s,a1,a2)
c     B1 Veltman function derivative
      implicit double precision (a-h,o-z)
      double complex f,df
      logical infstat
      common/renorm/del,amiu2,infstat
      if (infstat) then
        db1 = 0
        return
      end if
      if (s.eq.0.d0) then
        if (abs(a1) + abs(a2).eq.0.d0) then
          stop 'DB1 called for s=m1=m2=0'
        else if (a1.eq.0.d0) then
          db1 = - 1.d0/6/a2/a2
        else if (a2.eq.0.d0) then
          db1 = - 1.d0/3/a1/a1
        else if (a1.eq.a2) then
          db1 = - 1.d0/12/a1/a2
        else
          dif = a1*a1 - a2*a2
          db1 = - (3*a1*a1*(dif*(a1*a1 + a2*a2)
     1        - 4*a1*a1*a2*a2*log(a1/a2)) - dif**3)/6/dif**4
        end if
        return
      end if
      if (a1.eq.a2) then
        db1 = - df(s,a1,a2)/2
      else
        dif = a1*a1 - a2*a2
        db1 = - ((s + dif)*df(s,a1,a2) - dif/s*f(s,a1,a2))/2/s
      end if
      return
      end

      double complex function db22(s,a1,a2)
c     B22 Veltman function derivative
      implicit double precision (a-h,o-z)
      double complex f,df
      logical infstat
      common/spence/ber(9),pi,pi6,eps
      common/renorm/del,amiu2,infstat
      if (infstat) then
        db22 = - del/12
        return
      end if
      if (s.eq.0.d0) then
        if (abs(a1) + abs(a2).eq.0.d0) then
          stop 'DB22 called for s=m1=m2=0'
        else if (a1*a2.eq.0.d0) then
          db22 = - (del - log((a1 + a2)*(a1 + a2)/amiu2) + 5.d0/6)/12
        else
          db22 = - (del - log(a1*a2/amiu2))/12
          if (a1.ne.a2) then
            sum  = a1*a1 + a2*a2
            dif  = a1*a1 - a2*a2
            prod = a1*a1*a2*a2
            al   = log(a1/a2)
            db22 = db22 - (5.d0/3 - sum/dif*al + (12*sum/dif*prod*al
     1           - 3*sum*sum - 2*dif*dif)/6/dif/dif)/12
          end if
        end if
        return
      end if
      if (s.eq.(a1 + a2)**2) then
        if (a1*a2.eq.0.d0) then
          db22 = - (del - log(s/amiu2) + 5.d0/3)/12
        else
          db22 = - (del - log(a1*a2/amiu2) + 8.d0/3)/12
          if (a1.ne.a2) then
            sum  = a1*a1 + a2*a2
            dif  = a1*a1 - a2*a2
            prod = a1*a1*a2*a2
            al   = log(a1/a2)
            db22 = db22 - ((8*prod - s*sum)/dif/s*al
     1           - (a1 - a2)*(a1 - a2)/s)/12
          end if
        end if
        return
      end if
      if ((a1*a2).eq.0.d0) then
        if (a1.eq.a2) then
          db22 = - (del - log(abs(s)/amiu2) + 5.d0/3)/12
          if (s.gt.0.d0) db22 = db22 - (0.d0,1.d0)*pi/12
        else
          as = a1 + a2
          db22 = - (del - log(as*as/amiu2) + 5.d0/3
     1         + (1 - as**4/s/s)*f(s,a1,a2)
     2         + (s - as*as)**2/s*df(s,a1,a2))/12
        end if
        return
      end if
      if (a1.eq.a2) then
        db22 = - (del - log(a1*a1/amiu2) + 2.d0/3 + f(s,a1,a2)
     1       + (s - 4*a1*a2)*df(s,a1,a2))/12
      else
        sum = a1*a1 + a2*a2
        dif = a1*a1 - a2*a2
        db22 = - (del - log(a1*a2/amiu2) + 5.d0/3 - sum/dif*log(a1/a2)
     1         + (1 - dif*dif/s/s)*f(s,a1,a2)
     2         + (s*s - 2*s*sum + dif**2)/s*df(s,a1,a2))/12
      end if
      return
      end




