c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor
      
c     FILENAME: ROMBINT.F
c     Released: 19:03:1999(J.R.)
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      Simple numerical integration procedure: Romberg method       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function romberg(fun,a,b,errin,errout,niter)
c      a,b:      integration limits
c      fun:      integrated function
c      errin:    required relative error of integration
c      errout:   real (achieved) integration error
c      niter:    maximal number of iterations (maximum 20)
c      aint:     result of integration
      parameter(maxint=20)
      implicit double precision (a-h,o-z)
      logical rstat
      dimension w(maxint)
      common/romb_coeff/q(maxint),rstat
      external fun,romb_stat
      if (.not.rstat) call init_romberg
      n = min(maxint,niter)
      oldint = (b-a)/4.d0*(fun(a) + 2*fun((a+b)/2) + fun(b))
      w(1) = oldint
      do i=2,n
        call trapez(a,b,i,fun,pint,oldint)
        oldint = pint
        call richards(w,pint,i,errout,aint)
        romberg = aint
        if (errout.le.errin) return
      end do
      return
      end
      
      subroutine trapez(a,b,i,fun,pint,oldint)
      implicit double precision (a-h,o-z)
      n = 2**i
      h = (b - a)/dble(n)
      pint = 0.d0
      do k = 1,n-1,2
        x = a + h*dble(k)
        pint = pint + fun(x)
      end do
      pint = oldint/2 + h*pint
      return
      end
      
      subroutine richards(w,pint,i,err,aint)
      parameter(maxint=20)
      implicit double precision (a-h,o-z)
      logical rstat
      dimension w(maxint)
      common/romb_coeff/q(maxint),rstat
      a = w(1)
      b = a
      w(1) = pint
      do n=1,i-1
        if (n.lt.i-1) b = w(n+1)
        w(n+1) = w(n) + q(n)*(w(n) - a)
        a = b
      end do
      aint = w(i)
      err = abs(w(i-1) - a)/(abs(aint) + 1.d-7)
      return
      end
      
      subroutine init_romberg
      parameter(maxint=20)
      implicit double precision (a-h,o-z)
      logical rstat
      common/romb_coeff/q(maxint),rstat
      rstat = .true.
      do i=1,maxint-1
        q(i)=1/(4.d0**i-1)
      end do
      return
      end
      
      block data romb_stat
      parameter(maxint=20)
      implicit double precision (a-h,o-z)
      logical rstat
      common/romb_coeff/q(maxint),rstat
      data rstat/.false./
      end
         


