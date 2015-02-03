c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor
      
c     FILENAME: CDM_G.F
c     Released: 12.01.1999 (J.R.)
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expressions for CDM of gluons             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
c     Definition of VEGAS integrand
      double precision function cdm_g_vegas(x)
      implicit double precision (a-h,o-z)
      dimension x(3)
      common/cdm_g_aux/z1,z2,z3,errin,eps,npoint,itmax
      cdm_g_vegas = x(1)*(1 - x(1))*x(2)/2
     $    *(x(2)*(1 - x(1)) + z3*x(1)*(1 - x(1))*(1 - x(2))
     $    - 2*x(1)*x(2)*(z1*x(3) + z2*(1 - x(3))))
     $    *((1 - 2*x(1) + 8.d0*x(1)*x(1)/9)*(1 - x(2))*(1 - x(2)) 
     $    + x(2)*x(2))
     $    /(x(2)*(1 - x(1)) + z3*x(1)*(1 - x(1))*(1 - x(2)) 
     $    + x(1)*x(2)*(z1*x(3) + z2*(1 - x(3))))**4
      return
      end
      
      double precision function cdm_g()
c     Full gluon CDM (only stops/gluino contribution included)
      implicit double precision (a-h,o-z)
      double complex zu0,zd0
      double complex zu
      double complex gm2,gm3
      double precision cdm_g_vegas
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/cdm_g_aux/z1,z2,z3,errin,eps,npoint,itmax
      common/vegas_result/avgi,sd,ti,tsi,nnew
      common/gmass/gm1,gm2,gm3
      common/fmass/em(3),um(3),dm(3)
      common/debug_4q/ih,ic,in,ig
      external cdm_g_vegas,init_cdm_g
      if (ig.ne.1) then  
         cdm_g = 0.d0
         return
      end if
      if (init_alpha_susy) call init_alpha_s_susy
c     identify left and right stop
      zlmax = 0
      zrmax = 0
      do i=1,6
        if (abs(zu(3,i)).gt.zlmax) then
          zlmax = abs(zu(3,i))
          il = i
        end if
        if (abs(zu(6,i)).gt.zrmax) then
          zrmax = abs(zu(6,i))
          ir = i
        end if
      end do
      z1 = (sum(il)/gm1)**2
      z2 = (sum(ir)/gm1)**2
      z3 = (um(3)/gm1)**2
      call vegas(cdm_g_vegas,errin,3,npoint,itmax,nprn,igraph)
      cdm_g = 3.d0*pi/8*um(3)*(g3u/2/pi)**5/gm1/gm1/gm1
     $     *(z1 - z2)*dimag(zu(6,ir)*dconjg(zu(3,ir)))*avgi
      return
      end
      
      double precision function cdm_g_stop()
c     Full gluon CDM (only stops contribution for now...)
c     separate stop mass matrix diagonalization
      implicit double precision (a-h,o-z)
      double complex gm2,gm3
      double complex zst(2,2),h,tmp(2,2)
      double complex lms,rms,ums,dms,qms
      double complex ls,ks,ds,es,us,ws
      double complex yl,yu,yd
      double precision au(2,2),bu(2,2),zu1(2,2),zu2(2,2),work(4),stm(2)
      double precision cdm_g_vegas
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/msoft/lms(3,3),rms(3,3),ums(3,3),dms(3,3),qms(3,3)
      common/hpar/hm1,hm2,hs,h
      common/vev/v1,v2
      common/yukawa/yl(3),yu(3),yd(3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/cdm_g_aux/z1,z2,z3,errin,eps,npoint,itmax
      common/vegas_result/avgi,sd,ti,tsi,nnew
      common/gmass/gm1,gm2,gm3
      common/fmass/em(3),um(3),dm(3)
      common/debug_4q/ih,ic,in,ig
      external cdm_g_vegas,init_cdm_g
      if (ig.ne.1) then  
         cdm_g_stop = 0.d0
         return
      end if
      vmin = v1*v1 - v2*v2
c     Stop mass matrix initialization
      tmp(1,1) = dconjg(qms(3,3)) - e2*vmin*(1-4*ct2)/24/sct2 + um(3)**2
      tmp(2,2) = dconjg(ums(3,3)) + e2/6/ct2*vmin + um(3)**2
      tmp(1,2) = - dconjg(v2*us(3,3) + v1*ws(3,3))/sq2 - v1*yu(3)/sq2*h
      tmp(2,1) = dconjg(tmp(1,2))
      do i=1,2
        do j=1,2
          au(i,j) = dble(tmp(i,j))
          bu(i,j) = dimag(tmp(i,j))
        end do
      end do
      call eisch1(2,2,au,bu,stm,zu1,zu2,ierr,work)
      do i=1,2
        if (stm(i).le.0.d0) return
        stm(i) = sign(sqrt(abs(stm(i))),stm(i))
        do j=1,2
          zst(i,j) = dcmplx(zu1(i,j),zu2(i,j))
        end do
      end do
      z1 = (stm(1)/gm1)**2
      z2 = (stm(2)/gm1)**2
      z3 = (um(3)/gm1)**2
      call vegas(cdm_g_vegas,errin,3,npoint,itmax,nprn,igraph)
      cdm_g_stop = 3.d0*pi/8*um(3)*(alfas(zm)/pi)**2.5d0/gm1/gm1/gm1
     $    *(z1 - z2)*dimag(zst(2,2)*dconjg(zst(1,2)))*avgi
      return
      end
      
      block data init_cdm_g
      implicit double precision (a-h,o-z)
      common/cdm_g_aux/z1,z2,z3,errin,eps,npoint,itmax
      common/vegas_random/iseed
      data errin,eps,npoint,itmax/5.d-2,1.d-8,10000,100/
      data iseed/1/
      end
      
