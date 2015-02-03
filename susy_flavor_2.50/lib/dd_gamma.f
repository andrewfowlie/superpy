c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: DD_GAMMA.F
c     Released: 25:03:1996 (J.R.)
c     Revised: 15.10.2013(J.R.)
c     signs in gluino diagrams corrected
c     Revised: 25.10.2013(J.R.)
c     proper resummation of chiral corrections in W/Higgs diagrams done

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expressions for the coefficients of the   c
c     Hamiltonian for the d^J -> d^I + gamma decay, e.g.           c
c     b -> s gamma at the MZ energy scale                          c
c     General form of the Hamiltonian is:                          c
c     -iH = SUM_{i=1}^5 (A^i_L H^i_L + A^i_R H^i_R)                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dd_gam_w(i,j,cfl,cfr)
c     W and u quark contributions (SM like)
      implicit double precision (a-h,o-z)
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      double complex cfl(5),cfr(5)
      double complex ai,aj
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      qf = - 2.d0/3
      qv = 1
      do k=1,3
        ai = e/sq2/st*dconjg(ckm_phys(k,i))
        aj = e/sq2/st*ckm_phys(k,j)
        call dd_ffv(ai,aj,qf,umu(k),wm,cfl,cfr)
        call dd_vvf(ai,aj,qv,umu(k),wm,cfl,cfr)
      end do
      return
      end

      subroutine dd_gam_h(i,j,cfl,cfr)
c     Higgs and Goldstone contributions
      implicit double precision (a-h,o-z)
      double complex cfl(5),cfr(5)
      double complex ai,aj,bi,bj
      double complex yh_eff_l,yh_eff_r
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      qf = - 2.d0/3
      qs = - 1
      do k=1,3
        do l=1,2
           ai = dconjg(yh_eff_l(i,k,l)) ! zh(2,l)*yu(k)*dconjg(ckm_phys(k,i))
           bi = dconjg(yh_eff_r(i,k,l)) ! - zh(1,l)*dconjg(yd(i)*ckm_phys(k,i))
           aj = yh_eff_l(j,k,l) ! zh(2,l)*dconjg(yu(k))*ckm_phys(k,j)
           bj = yh_eff_r(j,k,l) ! - zh(1,l)*yd(j)*ckm_phys(k,j)
           call dd_ffs(ai,aj,bi,bj,qf,umu(k),cm(l),cfl,cfr)
           call dd_ssf(ai,aj,bi,bj,qs,umu(k),cm(l),cfl,cfr)
        end do
      end do
      return
      end

      subroutine dd_gam_c(i,j,cfl,cfr)
c     chargino/up squark contributions
      implicit double precision (a-h,o-z)
      double complex cfl(5),cfr(5)
      double complex ai,aj,bi,bj
      double complex vl_duc,vr_duc,zpos,zneg,zu,zd
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      qf = 1
      qs = 2.d0/3
      do k=1,6
        do l=1,2
          ai = - dconjg(vl_duc(i,k,l))
          bi = - dconjg(vr_duc(i,k,l))
          aj = - vl_duc(j,k,l)
          bj = - vr_duc(j,k,l)
          call dd_ffs(ai,aj,bi,bj,qf,fcm(l),sum(k),cfl,cfr)
          call dd_ssf(ai,aj,bi,bj,qs,fcm(l),sum(k),cfl,cfr)
        end do
      end do
      return
      end

      subroutine dd_gam_n(i,j,cfl,cfr)
c     neutralino-down squark contributions
      implicit double precision (a-h,o-z)
      double complex cfl(5),cfr(5)
      double complex ai,aj,bi,bj
      double complex vl_ddn,vr_ddn,zn,zu,zd
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      qs = - 1.d0/3
      do k=1,6
        do l=1,4
          ai = - dconjg(vl_ddn(i,k,l))
          bi = - dconjg(vr_ddn(i,k,l))
          aj = - vl_ddn(j,k,l)
          bj = - vr_ddn(j,k,l)
          call dd_ssf(ai,aj,bi,bj,qs,fnm(l),sdm(k),cfl,cfr)
        end do
      end do
      return
      end

      subroutine dd_gam_g(i,j,cfl,cfr)
c     gluino-down squark contributions
      implicit double precision (a-h,o-z)
      double complex cfl(5),cfr(5)
      double complex ai,aj,bi,bj
      double complex zu0,zd0,gm2,gm3
      double complex zd
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      qs = - 8/9.d0*g3d*g3d 
      do k=1,6
         ai = - dconjg(zd(i,k))
         bi =   dconjg(zd(i+3,k))
         aj = - zd(j,k)
         bj =   zd(j+3,k)
         call dd_ssf(ai,aj,bi,bj,qs,gm1,sdm(k),cfl,cfr)
      end do
      return
      end

      subroutine dd_gam(i,j,cfl,cfr)
c     Full coefficients
      implicit double precision (a-h,o-z)
      double complex cfl(5),cfr(5)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      do k=1,5
        cfl(k) = (0.d0,0.d0)
        cfr(k) = (0.d0,0.d0)
      end do
      if (ih.eq.1) then
         call dd_gam_w(i,j,cfl,cfr)
         call dd_gam_h(i,j,cfl,cfr)
      end if
      if (ic.eq.1) call dd_gam_c(i,j,cfl,cfr)
      if (in.eq.1) call dd_gam_n(i,j,cfl,cfr)
      if (ig.eq.1) call dd_gam_g(i,j,cfl,cfr)
      do k=1,5
        cfl(k) = cfl(k)/16/pi/pi
        cfr(k) = cfr(k)/16/pi/pi
      end do
      return
      end


