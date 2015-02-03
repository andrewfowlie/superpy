c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: LL_GAMMA.F
c     Released: 11:02:2012 (J.R.)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expressions for the coefficients of the   c
c     Hamiltonian for the l^J -> l^I + gamma decay, e.g.           c
c     mu -> e gamma at the MZ energy scale                         c
c     General form of the Hamiltonian is:                          c
c     -iH = SUM_{i=1}^5 (cfl^i_L H^i_L + cfr^i_R H^i_R)            c
c     Small SM-like W or Higgs corrections from U_MNS not included c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine ll_gam_c(i,j,cfl,cfr)
c     chargino-sneutrino contributions
      implicit double precision (a-h,o-z)
      double complex cfl(5),cfr(5)
      double complex ai,aj,bi,bj
      double complex zpos,zneg,zv,zl
      double complex yl,yu,yd
      common/yukawa/yl(3),yu(3),yd(3)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      qf = 1
      do k=1,3
        do l=1,2
          ai = e/st*zv(i,k)*dconjg(zpos(1,l))
          bi = dconjg(yl(i))*zv(i,k)*zneg(2,l)
          aj = e/st*dconjg(zv(j,k))*zpos(1,l) 
          bj = yl(j)*dconjg(zv(j,k)*zneg(2,l))
          call dd_ffs(ai,aj,bi,bj,qf,fcm(l),vm(k),cfl,cfr)
        end do
      end do
      return
      end

      subroutine ll_gam_n(i,j,cfl,cfr)
c     neutralino-slepton contributions
      implicit double precision (a-h,o-z)
      double complex cfl(5),cfr(5)
      double complex ai,aj,bi,bj
      double complex vl_lln,vr_lln,zn,zv,zl
      common/neut/fnm(4),zn(4,4)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      qs = - 1
      do k=1,6
        do l=1,4
          ai = dconjg(vl_lln(i,k,l))
          bi = dconjg(vr_lln(i,k,l))
          aj = vl_lln(j,k,l)
          bj = vr_lln(j,k,l)
          call dd_ssf(ai,aj,bi,bj,qs,fnm(l),slm(k),cfl,cfr)
        end do
      end do
      return
      end

      subroutine ll_gam(i,j,cfl,cfr)
c     Full coefficients
      implicit double precision (a-h,o-z)
      double complex cfl(5),cfr(5)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      do k=1,5
        cfl(k) = (0.d0,0.d0)
        cfr(k) = (0.d0,0.d0)
      end do
      if (ic.eq.1) call ll_gam_c(i,j,cfl,cfr)
      if (in.eq.1) call ll_gam_n(i,j,cfl,cfr)
      do k=1,5
        cfl(k) = cfl(k)/16/pi/pi
        cfr(k) = cfr(k)/16/pi/pi
      end do
      return
      end


