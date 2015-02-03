c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: DD_MIX.F
c     Released: 25:03:1995 (J.R.)
c     Revised: 25:09:1999 (J.R.,P.Ch.)
c     Errors in non-SM formfactors corrected.
c     Revised: 21.02.2001 (J.R.)
c     alpha_s(msusy) used in gluino contribution calculations.
c     Revised: 02.04.2001 (J.R.)
c     Effective chargino couplings introduced.
c     Revised: 14.03.2002 (J.R.)
c     Double neutral Higgs penguins added.
c     Effective charged Higgs couplings introduced.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expression for the coefficient of the     c
c     Hamiltonian for the B-Bbar and the K-Kbar mixing (at the     c
c     quark level).                                                c
c     General form of the Hamiltonian is:                          c
c     H = A^V_LL H^V_LL + A^V_RR H^V_RR + A^V_LR H^V_LR            c
c       + A^S_LL H^S_LL + A^S_RR H^S_RR + A^S_LR H^S_LR            c
c       + A^T_L H^T_L   + A^T_R H^T_R                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^V_LL                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dd_vll_sm(i,j)
c     Full A^V_LL SM formfactor
      implicit double precision (a-h,o-z)
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      common/debug_4q/ih,ic,in,ig
      dd_vll_sm = (0.d0,0.d0)
      if (ih.ne.1) return
      do m=1,3
         do n=1,3
            dd_vll_sm = dd_vll_sm + e2*e2/4/st2/st2
     $           * ckm_phys(m,i)*ckm_phys(n,i)
     $           * dconjg(ckm_phys(m,j)*ckm_phys(n,j))
     $           * ((1 + (umu(m)*umu(n)/wm2)**2/4)/2
     $           * dp1(wm,wm,umu(m),umu(n)) 
     $           - (umu(m)*umu(n)/wm)**2*dp0(wm,wm,umu(m),umu(n)))
         end do
      end do
      return
      end
      
      double complex function dd_vll_hg(i,j)
c     Charged Higgs contributions
      implicit double precision (a-h,o-z)
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      double complex yh_eff_l
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      dd_vll_hg = (0.d0,0.d0)
      do m=1,3
         do n=1,3
            dd_vll_hg = dd_vll_hg - e2/st2/2*yh_eff_l(i,m,1)
     $           * dconjg(yh_eff_l(j,n,1)*ckm_phys(m,j))*ckm_phys(n,i)
     $           * umu(m)*umu(n)*dp0(wm,cm(1),umu(m),umu(n))
            dd_vll_hg = dd_vll_hg + yh_eff_l(i,m,1)*yh_eff_l(i,n,1)
     $           * dconjg(yh_eff_l(j,n,1)*yh_eff_l(j,m,1))/8
     $           * dp1(cm(1),cm(1),umu(m),umu(n))
     $           + yh_eff_l(i,m,1)*yh_eff_l(i,n,2)
     $           * dconjg(yh_eff_l(j,n,1)*yh_eff_l(j,m,2))/4
     $           * dp1(wm,cm(1),umu(m),umu(n))
         end do
      end do
      return
      end

      double complex function dd_vll_c(i,j)
c     chargino contributions
      implicit double precision (a-h,o-z)
      double complex vl_duc,zpos,zneg,zu,zd
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      dd_vll_c = (0.d0,0.d0)
      do m=1,2
         do n=1,2
            do k=1,6
               do l=1,6
                  dd_vll_c = dd_vll_c 
     $                 + vl_duc(i,k,m)*vl_duc(i,l,n)
     $                 * dconjg(vl_duc(j,k,n)*vl_duc(j,l,m))/8
     $                 * dp1(fcm(m),fcm(n),sum(k),sum(l))
               end do
            end do
         end do
      end do
      return
      end

      double complex function dd_vll_n(i,j)
c     neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_ddn,zn,zu,zd
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      dd_vll_n = (0.d0,0.d0)
      do m=1,4
         do n=1,4
            do k=1,6
               do l=1,6
                  dd_vll_n = dd_vll_n 
     $                 + vl_ddn(i,k,m)*vl_ddn(i,l,n)
     $                 * dconjg(vl_ddn(j,k,n)*vl_ddn(j,l,m))/8
     $                 * dp1(fnm(m),fnm(n),sdm(k),sdm(l))
                  dd_vll_n = dd_vll_n 
     $                 + vl_ddn(i,k,m)*vl_ddn(i,l,m)
     $                 * dconjg(vl_ddn(j,l,n)*vl_ddn(j,k,n))/4
     $                 * fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sdm(k),sdm(l))
               end do
            end do
         end do
      end do
      return
      end

      double complex function dd_vll_ng(i,j)
c     gluino-neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_ddn,zn,zu0,zd0,gm2,gm3
      double complex zd
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/neut/fnm(4),zn(4,4)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      dd_vll_ng = (0.d0,0.d0)
      do m=1,4
         do k=1,6
            do l=1,6
               dd_vll_ng = dd_vll_ng + g3d/6*vl_ddn(i,k,m)*zd(i,l)
     $              * dconjg(vl_ddn(j,l,m)*zd(j,k))
     $              * dp1(gm1,fnm(m),sdm(k),sdm(l))
               dd_vll_ng = dd_vll_ng 
     $              + g3d/6*(vl_ddn(i,k,m)*vl_ddn(i,l,m)
     $              * dconjg(zd(j,k)*zd(j,l)) + zd(i,k)*zd(i,l)
     $              * dconjg(vl_ddn(j,k,m)*vl_ddn(j,l,m)))
     $              * gm1*fnm(m)*dp0(gm1,fnm(m),sdm(k),sdm(l))
            end do
         end do
      end do
      return
      end
      
      double complex function dd_vll_g(i,j)
c     gluino contributions
      implicit double precision (a-h,o-z)
      double complex zu0,zd0,gm2,gm3
      double complex zd
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      dd_vll_g = (0.d0,0.d0)
      do k=1,6
         do l=1,6
            dd_vll_g = dd_vll_g 
     $           + g3d*g3d/36*zd(i,k)*zd(i,l)*dconjg(zd(j,k)*zd(j,l))
     $           * (11*dp1(gm1,gm1,sdm(k),sdm(l))
     $           + 4*gm1*gm1*dp0(gm1,gm1,sdm(k),sdm(l)))
         end do
      end do
      return
      end
      
      double complex function dd_vll(i,j)
c     Full A^V_LL formfactor
      implicit double precision (a-h,o-z)
c      double complex dd_vll_sm
      double complex dd_vll_hg,dd_vll_c,dd_vll_n,dd_vll_ng,dd_vll_g
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      dd_vll = (0.d0,0.d0)
      if (ih.eq.1)  dd_vll = dd_vll + dd_vll_hg(i,j) 
      if (ic.eq.1)  dd_vll = dd_vll + dd_vll_c(i,j) 
      if (in.eq.1)  dd_vll = dd_vll + dd_vll_n(i,j) 
      if ((in.eq.1).and.(ig.eq.1)) dd_vll = dd_vll + dd_vll_ng(i,j) 
      if (ig.eq.1)  dd_vll = dd_vll + dd_vll_g(i,j) 
      dd_vll = - dd_vll/16/pi/pi
c     SM contribution has to be treated separately from the point of 
c     view of QCD evolution factors. Therefore we do not add it here.
c     Perhaps the same should be done for Higgs contribution dd_vll_hg?
c     In order to calculate the "full" high-scale dd_vll, uncomment:
c     dd_vll = dd_vll - dd_vll_sm(i,j)/16/pi/pi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^V_RR                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dd_vrr_hg(i,j)
c     Higgs and gauge contributions
      implicit double precision (a-h,o-z)
      double complex yh_eff_r
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
      dd_vrr_hg = (0.d0,0.d0)
      do m=1,3
         do n=1,3
            do k=1,2
               do l=1,2
                  dd_vrr_hg = dd_vrr_hg 
     $                 + yh_eff_r(i,m,k)*yh_eff_r(i,n,l)
     $                 * dconjg(yh_eff_r(j,n,k)*yh_eff_r(j,m,l))
     $                 * dp1(cm(k),cm(l),umu(m),umu(n))
               end do
            end do
         end do
      end do
      dd_vrr_hg = dd_vrr_hg/8
      return
      end

      double complex function dd_vrr_c(i,j)
c     chargino contributions
      implicit double precision (a-h,o-z)
      double complex vr_duc,zpos,zneg,zu,zd
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      dd_vrr_c = (0.d0,0.d0)
      do m=1,2
         do n=1,2
            do k=1,6
               do l=1,6
                  dd_vrr_c = dd_vrr_c 
     $                 + vr_duc(i,k,m)*vr_duc(i,l,n)
     $                 * dconjg(vr_duc(j,k,n)*vr_duc(j,l,m))/8
     $                 * dp1(fcm(m),fcm(n),sum(k),sum(l))
               end do
            end do
         end do
      end do
      return
      end

      double complex function dd_vrr_n(i,j)
c     neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vr_ddn,zn,zu,zd
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      dd_vrr_n = (0.d0,0.d0)
      do m=1,4
         do n=1,4
            do k=1,6
               do l=1,6
                  dd_vrr_n = dd_vrr_n 
     $                 + vr_ddn(i,k,m)*vr_ddn(i,l,n)
     $                 * dconjg(vr_ddn(j,k,n)*vr_ddn(j,l,m))/8
     $                 * dp1(fnm(m),fnm(n),sdm(k),sdm(l))
                  dd_vrr_n = dd_vrr_n 
     $                 + vr_ddn(i,k,m)*vr_ddn(i,l,m)
     $                 * dconjg(vr_ddn(j,l,n)*vr_ddn(j,k,n))/4
     $                 * fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sdm(k),sdm(l))
               end do
            end do
         end do
      end do
      return
      end

      double complex function dd_vrr_ng(i,j)
c     gluino-neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vr_ddn,zn,zu0,zd0,gm2,gm3
      double complex zd
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/neut/fnm(4),zn(4,4)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      dd_vrr_ng = (0.d0,0.d0)
      do m=1,4
         do k=1,6
            do l=1,6
               dd_vrr_ng = dd_vrr_ng + g3d/6*vr_ddn(i,k,m)*zd(i+3,l)
     $              * dconjg(vr_ddn(j,l,m)*zd(j+3,k))
     $              * dp1(gm1,fnm(m),sdm(k),sdm(l))
               dd_vrr_ng = dd_vrr_ng 
     $              + g3d/6*(vr_ddn(i,k,m)*vr_ddn(i,l,m)
     $              * dconjg(zd(j+3,k)*zd(j+3,l)) + zd(i+3,k)*zd(i+3,l)
     $              * dconjg(vr_ddn(j,k,m)*vr_ddn(j,l,m)))
     $              * gm1*fnm(m)*dp0(gm1,fnm(m),sdm(k),sdm(l))
            end do
         end do
      end do
      return
      end

      double complex function dd_vrr_g(i,j)
c      gluino contributions
      implicit double precision (a-h,o-z)
      double complex zu0,zd0,gm2,gm3
      double complex zd
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      dd_vrr_g = (0.d0,0.d0)
      do k=1,6
         do l=1,6
            dd_vrr_g = dd_vrr_g + g3d*g3d/36*zd(i+3,k)*zd(i+3,l)
     $           * dconjg(zd(j+3,k)*zd(j+3,l))
     $           * (11*dp1(gm1,gm1,sdm(k),sdm(l))
     $           + 4*gm1*gm1*dp0(gm1,gm1,sdm(k),sdm(l)))
         end do
      end do
      return
      end
      
      double complex function dd_vrr(i,j)
c      Full A^V_RR formfactor
      implicit double precision (a-h,o-z)
      double complex dd_vrr_hg,dd_vrr_c,dd_vrr_n,dd_vrr_ng,dd_vrr_g
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      dd_vrr = (0.d0,0.d0)
      if (ih.eq.1)  dd_vrr = dd_vrr + dd_vrr_hg(i,j) 
      if (ic.eq.1)  dd_vrr = dd_vrr + dd_vrr_c(i,j) 
      if (in.eq.1)  dd_vrr = dd_vrr + dd_vrr_n(i,j) 
      if ((in.eq.1).and.(ig.eq.1)) dd_vrr = dd_vrr + dd_vrr_ng(i,j) 
      if (ig.eq.1)  dd_vrr = dd_vrr + dd_vrr_g(i,j) 
      dd_vrr = - dd_vrr/16/pi/pi
      return
      end
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^V_LR                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dd_vlr_hg(i,j)
c      Higgs and gauge contributions
      implicit double precision (a-h,o-z)
      double complex yh_eff_l,yh_eff_r
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
      dd_vlr_hg = (0.d0,0.d0)
      do m=1,3
         do n=1,3
            do k=1,2
               do l=1,2
                  dd_vlr_hg = dd_vlr_hg
     $                 + yh_eff_l(i,m,k)*yh_eff_r(i,n,l) 
     $                 * dconjg(yh_eff_l(j,m,l)*yh_eff_r(j,n,k))
     $                 * dp1(cm(k),cm(l),umu(m),umu(n))
               end do
            end do
         end do
      end do
      dd_vlr_hg = dd_vlr_hg/4
      return
      end

      double complex function dd_vlr_c(i,j)
c      chargino contributions
      implicit double precision (a-h,o-z)
      double complex vl_duc,vr_duc,zpos,zneg,zu,zd
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      dd_vlr_c = (0.d0,0.d0)
      do m=1,2
         do n=1,2
            do k=1,6
               do l=1,6
                  dd_vlr_c = dd_vlr_c 
     $                 - vl_duc(i,k,m)*vr_duc(i,l,n)
     $                 * dconjg(vl_duc(j,k,n)*vr_duc(j,l,m))/2
     $                 * fcm(m)*fcm(n)*dp0(fcm(m),fcm(n),sum(k),sum(l))
            end do
          end do
        end do
      end do
      return
      end

      double complex function dd_vlr_n(i,j)
c     neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_ddn,vr_ddn,zn,zu,zd
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      dd_vlr_n = (0.d0,0.d0)
      do m=1,4
         do n=1,4
            do k=1,6
               do l=1,6
                  dd_vlr_n = dd_vlr_n 
     $                 - vl_ddn(i,k,m)*vr_ddn(i,l,n)
     $                 * dconjg(vl_ddn(j,k,n)*vr_ddn(j,l,m))/2
     $                 * fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sdm(k),sdm(l))
                  dd_vlr_n = dd_vlr_n 
     $                 - vl_ddn(i,k,m)*vr_ddn(i,l,m)
     $                 * dconjg(vl_ddn(j,k,n)*vr_ddn(j,l,n))/4
     $                 * dp1(fnm(m),fnm(n),sdm(k),sdm(l))
               end do
            end do
         end do
      end do
      return
      end

      double complex function dd_vlr_ng(i,j)
c      gluino-neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_ddn,vr_ddn,zn,zu0,zd0,gm2,gm3
      double complex zd
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/neut/fnm(4),zn(4,4)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      dd_vlr_ng = (0.d0,0.d0)
      do m=1,4
         do k=1,6
            do l=1,6
               dd_vlr_ng = dd_vlr_ng 
     $              + g3d/4*dp1(gm1,fnm(m),sdm(k),sdm(l))*(zd(i,l)
     $              *vr_ddn(i,k,m)*dconjg(zd(j,k)*vr_ddn(j,l,m)) 
     $              + zd(i+3,l)*vl_ddn(i,k,m)*dconjg(zd(j+3,k)
     $              * vl_ddn(j,l,m)))
     $              - g3d/6*gm1*fnm(m)*dp0(gm1,fnm(m),sdm(k),sdm(l))
     $              * (zd(i,l)*vr_ddn(i,k,m)
     $              * dconjg(zd(j+3,k)*vl_ddn(j,l,m))
     $              + zd(i+3,l)*vl_ddn(i,k,m)
     $              * dconjg(zd(j,k)*vr_ddn(j,l,m)))
               dd_vlr_ng = dd_vlr_ng 
     $              - g3d/12*dp1(gm1,fnm(m),sdm(k),sdm(l))
     $              * (3*(zd(i,l)*zd(i+3,k)
     $              * dconjg(vl_ddn(j,k,m)*vr_ddn(j,l,m))
     $              + vl_ddn(i,k,m)*vr_ddn(i,l,m)
     $              * dconjg(zd(j,l)*zd(j+3,k)))
     $              + (zd(i+3,k)*zd(i,l)
     $              * dconjg(vl_ddn(j,l,m)*vr_ddn(j,k,m))
     $              + vr_ddn(i,k,m)*vl_ddn(i,l,m)
     $              * dconjg(zd(j,l)*zd(j+3,k))))
            end do
         end do
      end do
      return
      end

      double complex function dd_vlr_g(i,j)
c      gluino contributions
      implicit double precision (a-h,o-z)
      double complex zu0,zd0,gm2,gm3
      double complex zd
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      dd_vlr_g = (0.d0,0.d0)
      do k=1,6
         do l=1,6
            dd_vlr_g = dd_vlr_g 
     $           - gm1*gm1*g3d*g3d/18.d0*zd(i,k)*zd(i+3,l)
     $           * dconjg(zd(j,k)*zd(j+3,l))
     $           *dp0(gm1,gm1,sdm(k),sdm(l))
            dd_vlr_g = dd_vlr_g + 5*g3d*g3d/36.d0*zd(i,k)*zd(i+3,l)
     $           * dconjg(3*zd(j+3,k)*zd(j,l) -  2*zd(j,k)*zd(j+3,l))
     $           * dp1(gm1,gm1,sdm(k),sdm(l))
         end do
      end do
      return
      end

      double complex function dd_vlr(i,j)
c     Full A^V_LR formfactor
      implicit double precision (a-h,o-z)
      double complex dd_vlr_hg,dd_vlr_c,dd_vlr_n,dd_vlr_ng,dd_vlr_g
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      dd_vlr = (0.d0,0.d0)
      if (ih.eq.1)  dd_vlr = dd_vlr + dd_vlr_hg(i,j) 
      if (ic.eq.1)  dd_vlr = dd_vlr + dd_vlr_c(i,j) 
      if (in.eq.1)  dd_vlr = dd_vlr + dd_vlr_n(i,j) 
      if ((in.eq.1).and.(ig.eq.1)) dd_vlr = dd_vlr + dd_vlr_ng(i,j) 
      if (ig.eq.1)  dd_vlr = dd_vlr + dd_vlr_g(i,j) 
      dd_vlr = - dd_vlr/16/pi/pi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^S_LL                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dd_sll_hg(i,j)
c     Higgs and gauge contributions
      implicit double precision (a-h,o-z)
      double complex yh_eff_l,yh_eff_r
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
      dd_sll_hg = (0.d0,0.d0)
      do m=1,3
         do n=1,3
            do k=1,2
               do l=1,2
                  dd_sll_hg = dd_sll_hg 
     $                 + umu(m)*umu(n)*yh_eff_l(i,m,k)*yh_eff_l(i,n,l)
     $                 * dconjg(yh_eff_r(j,m,l)*yh_eff_r(j,n,k))
     $                 * dp0(cm(k),cm(l),umu(m),umu(n))
               end do
            end do
         end do
      end do
      dd_sll_hg = dd_sll_hg/2
      return
      end

      double complex function dd_sll_c(i,j)
c      chargino contributions
      implicit double precision (a-h,o-z)
      double complex vl_duc,vr_duc,zpos,zneg,zu,zd
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dd_sll_c = (0.d0,0.d0)
      do m=1,2
         do n=1,2
            do k=1,6
               do l=1,6
                  dd_sll_c = dd_sll_c 
     $                 - vl_duc(i,k,m)*vl_duc(i,l,n)
     $                 * dconjg(vr_duc(j,k,n)*vr_duc(j,l,m))/4
     $                 * fcm(m)*fcm(n)*dp0(fcm(m),fcm(n),sum(k),sum(l))
               end do
            end do
         end do
      end do
      return
      end

      double complex function dd_sll_n(i,j)
c      neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_ddn,vr_ddn,zn,zu,zd
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dd_sll_n = (0.d0,0.d0)
      do m=1,4
         do n=1,4
            do k=1,6
               do l=1,6
                  dd_sll_n = dd_sll_n 
     $                 - vl_ddn(i,k,m)*vl_ddn(i,l,n)
     $                 * dconjg(vr_ddn(j,k,n)*vr_ddn(j,l,m))/4
     $                 * fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sdm(k),sdm(l))
                  dd_sll_n = dd_sll_n 
     $                 - vl_ddn(i,k,m)*vl_ddn(i,l,m)
     $                 * dconjg(vr_ddn(j,k,n)*vr_ddn(j,l,n))/4
     $                 * fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sdm(k),sdm(l))
               end do
            end do
         end do
      end do
      return
      end

      double complex function dd_sll_ng(i,j)
c     gluino-neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_ddn,vr_ddn,zn,zu0,zd0,gm2,gm3
      double complex zd
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/neut/fnm(4),zn(4,4)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      dd_sll_ng = (0.d0,0.d0)
      do m=1,4
         do k=1,6
            do l=1,6
               dd_sll_ng = dd_sll_ng 
     $              - 7*g3d/6.d0*zd(i,l)*vl_ddn(i,k,m)
     $              * dconjg(vr_ddn(j,l,m)*zd(j+3,k))
     $              * gm1*fnm(m)*dp0(gm1,fnm(m),sdm(k),sdm(l))
               dd_sll_ng = dd_sll_ng - g3d/6*gm1*fnm(m)
     $              * dp0(gm1,fnm(m),sdm(k),sdm(l))
     $              * (zd(i,k)*zd(i,l)
     $              * dconjg(vr_ddn(j,k,m)*vr_ddn(j,l,m))
     $              + vl_ddn(i,k,m)*vl_ddn(i,l,m)
     $              * dconjg(zd(j+3,k)*zd(j+3,l))) 
            end do
         end do
      end do
      return
      end

      double complex function dd_sll_g(i,j)
c     gluino contributions
      implicit double precision (a-h,o-z)
      double complex zu0,zd0,gm2,gm3
      double complex zd
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      dd_sll_g = (0.d0,0.d0)
      do k=1,6
         do l=1,6
            dd_sll_g = dd_sll_g + 37*g3d*g3d/36.d0*zd(i,k)*zd(i,l)
     $           * dconjg(zd(j+3,k)*zd(j+3,l))
     $           * gm1*gm1*dp0(gm1,gm1,sdm(k),sdm(l))
         end do
      end do
      return
      end
      
      double complex function dd_sll(i,j)
c     Full A^S_LL formfactor
      implicit double precision (a-h,o-z)
      double complex dd_sll_hg,dd_sll_c,dd_sll_n,dd_sll_ng,dd_sll_g
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      dd_sll = (0.d0,0.d0)
      if (ih.eq.1)  dd_sll = dd_sll + dd_sll_hg(i,j) 
      if (ic.eq.1)  dd_sll = dd_sll + dd_sll_c(i,j) 
      if (in.eq.1)  dd_sll = dd_sll + dd_sll_n(i,j) 
      if ((in.eq.1).and.(ig.eq.1)) dd_sll = dd_sll + dd_sll_ng(i,j) 
      if (ig.eq.1)  dd_sll = dd_sll + dd_sll_g(i,j) 
      dd_sll = - dd_sll/16/pi/pi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^S_RR                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dd_srr_hg(i,j)
c     Higgs and gauge contributions
      implicit double precision (a-h,o-z)
      double complex yh_eff_l,yh_eff_r
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
      dd_srr_hg = (0.d0,0.d0)
      do m=1,3
         do n=1,3
            do k=1,2
               do l=1,2
                  dd_srr_hg = dd_srr_hg 
     $                 + umu(m)*umu(n)*yh_eff_r(i,m,k)*yh_eff_r(i,n,l)
     $                 * dconjg(yh_eff_l(j,m,l)*yh_eff_l(j,n,k))
     $                 * dp0(cm(k),cm(l),umu(m),umu(n))
               end do
            end do
         end do
      end do
      dd_srr_hg = dd_srr_hg/2
      return
      end

      double complex function dd_srr_c(i,j)
c     chargino contributions
      implicit double precision (a-h,o-z)
      double complex vl_duc,vr_duc,zpos,zneg,zu,zd
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dd_srr_c = (0.d0,0.d0)
      do m=1,2
         do n=1,2
            do k=1,6
               do l=1,6
                  dd_srr_c = dd_srr_c 
     $                 - vr_duc(i,k,m)*vr_duc(i,l,n)
     $                 * dconjg(vl_duc(j,k,n)*vl_duc(j,l,m))/4
     $                 * fcm(m)*fcm(n)*dp0(fcm(m),fcm(n),sum(k),sum(l))
               end do
            end do
         end do
      end do
      return
      end

      double complex function dd_srr_n(i,j)
c     neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_ddn,vr_ddn,zn,zu,zd
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dd_srr_n = (0.d0,0.d0)
      do m=1,4
         do n=1,4
            do k=1,6
               do l=1,6
                  dd_srr_n = dd_srr_n 
     $                 - vr_ddn(i,k,m)*vr_ddn(i,l,n)
     $                 * dconjg(vl_ddn(j,k,n)*vl_ddn(j,l,m))/4
     $                 * fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sdm(k),sdm(l))
                  dd_srr_n = dd_srr_n 
     $                 - vr_ddn(i,k,m)*vr_ddn(i,l,m)
     $                 * dconjg(vl_ddn(j,k,n)*vl_ddn(j,l,n))/4
     $                 * fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sdm(k),sdm(l))
               end do
            end do
         end do
      end do
      return
      end

      double complex function dd_srr_ng(i,j)
c     gluino-neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_ddn,vr_ddn,zn,zu0,zd0,gm2,gm3
      double complex zd
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/neut/fnm(4),zn(4,4)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      dd_srr_ng = (0.d0,0.d0)
      do m=1,4
         do k=1,6
            do l=1,6
               dd_srr_ng = dd_srr_ng 
     $              - 7*g3d/6.d0*zd(i+3,l)*vr_ddn(i,k,m)
     $              * dconjg(vl_ddn(j,l,m)*zd(j,k))
     $              * gm1*fnm(m)*dp0(gm1,fnm(m),sdm(k),sdm(l))
               dd_srr_ng = dd_srr_ng - g3d/6*gm1*fnm(m)
     $              * dp0(gm1,fnm(m),sdm(k),sdm(l))
     $              *(zd(i+3,k)*zd(i+3,l)
     $              * dconjg(vl_ddn(j,k,m)*vl_ddn(j,l,m))
     $              + vr_ddn(i,k,m)*vr_ddn(i,l,m)
     $              * dconjg(zd(j,k)*zd(j,l)))
            end do
         end do
      end do
      return
      end

      double complex function dd_srr_g(i,j)
c      gluino contributions
      implicit double precision (a-h,o-z)
      double complex zu0,zd0,gm2,gm3
      double complex zd
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      dd_srr_g = (0.d0,0.d0)
      do k=1,6
         do l=1,6
            dd_srr_g = dd_srr_g + 37*g3d*g3d/36.d0*zd(i+3,k)*zd(i+3,l)
     $           *dconjg(zd(j,k)*zd(j,l))
     $           *gm1*gm1*dp0(gm1,gm1,sdm(k),sdm(l))
         end do
      end do
      return
      end

      double complex function dd_srr(i,j)
c     Full A^S_RR formfactor
      implicit double precision (a-h,o-z)
      double complex dd_srr_hg,dd_srr_c,dd_srr_n,dd_srr_ng,dd_srr_g
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      dd_srr = (0.d0,0.d0)
      if (ih.eq.1)  dd_srr = dd_srr + dd_srr_hg(i,j) 
      if (ic.eq.1)  dd_srr = dd_srr + dd_srr_c(i,j) 
      if (in.eq.1)  dd_srr = dd_srr + dd_srr_n(i,j) 
      if ((in.eq.1).and.(ig.eq.1)) dd_srr = dd_srr + dd_srr_ng(i,j) 
      if (ig.eq.1)  dd_srr = dd_srr + dd_srr_g(i,j) 
      dd_srr = - dd_srr/16/pi/pi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^S_LR                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dd_slr_hg(i,j)
c     Higgs and gauge contributions
      implicit double precision (a-h,o-z)
      double complex yh_eff_l,yh_eff_r
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
      dd_slr_hg = (0.d0,0.d0)
      do m=1,3
         do n=1,3
            do k=1,2
               dd_slr_hg = dd_slr_hg - e2/2/st2*yh_eff_r(i,m,k)
     $              *ckm_phys(n,i)*dconjg(ckm_phys(m,j)*yh_eff_r(j,n,k))
     $              * dp1(wm,cm(k),umu(m),umu(n))
               do l=1,2
                  dd_slr_hg = dd_slr_hg
     $                 + umu(m)*umu(n)*yh_eff_l(i,m,k)*yh_eff_r(i,n,l)
     $                 * dconjg(yh_eff_l(j,n,k)*yh_eff_r(j,m,l))
     $                 * dp0(cm(k),cm(l),umu(m),umu(n))
               end do
            end do
         end do
      end do
      return
      end
     
      double complex function dd_slr_c(i,j)
c     chargino contributions
      implicit double precision (a-h,o-z)
      double complex vl_duc,vr_duc,zpos,zneg,zu,zd
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dd_slr_c = (0.d0,0.d0)
      do m=1,2
         do n=1,2
            do k=1,6
               do l=1,6
                  dd_slr_c = dd_slr_c 
     $                 - vl_duc(i,k,m)*vr_duc(i,l,n)
     $                 * dconjg(vl_duc(j,l,m)*vr_duc(j,k,n))/2
     $                 * dp1(fcm(m),fcm(n),sum(k),sum(l))
               end do
            end do
         end do
      end do
      return
      end

      double complex function dd_slr_n(i,j)
c     neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_ddn,vr_ddn,zn,zu,zd
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dd_slr_n = (0.d0,0.d0)
      do m=1,4
         do n=1,4
            do k=1,6
               do l=1,6
                  dd_slr_n = dd_slr_n 
     $                 - vl_ddn(i,k,m)*vr_ddn(i,l,n)
     $                 * dconjg(vl_ddn(j,l,m)*vr_ddn(j,k,n))/2
     $                 * dp1(fnm(m),fnm(n),sdm(k),sdm(l))
                  dd_slr_n = dd_slr_n 
     $                 - vl_ddn(i,k,m)*vr_ddn(i,l,m)
     $                 * dconjg(vl_ddn(j,l,n)*vr_ddn(j,k,n))/2
     $                 * dp1(fnm(m),fnm(n),sdm(k),sdm(l))
               end do
            end do
         end do
      end do
      return
      end

      double complex function dd_slr_ng(i,j)
c     gluino-neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_ddn,vr_ddn,zn,zu0,zd0,gm2,gm3
      double complex zd
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/neut/fnm(4),zn(4,4)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      dd_slr_ng = (0.d0,0.d0)
      do m=1,4
         do k=1,6
            do l=1,6
               dd_slr_ng = dd_slr_ng 
     $              + g3d/6*dp1(gm1,fnm(m),sdm(k),sdm(l))
     $              * (zd(i,l)*vr_ddn(i,k,m)
     $              * dconjg(zd(j,k)*vr_ddn(j,l,m)) 
     $              + zd(i+3,l)*vl_ddn(i,k,m)
     $              * dconjg(zd(j+3,k)*vl_ddn(j,l,m)))
     $              - g3d*gm1*fnm(m)*dp0(gm1,fnm(m),sdm(k),sdm(l))
     $              * (zd(i,l)*vr_ddn(i,k,m)
     $              * dconjg(zd(j+3,k)*vl_ddn(j,l,m))
     $              + zd(i+3,l)*vl_ddn(i,k,m)
     $              * dconjg(zd(j,k)*vr_ddn(j,l,m)))
               dd_slr_ng = dd_slr_ng 
     $              - g3d/6*dp1(gm1,fnm(m),sdm(k),sdm(l))
     $              * ((zd(i,l)*zd(i+3,k)
     $              * dconjg(vl_ddn(j,k,m)*vr_ddn(j,l,m))
     $              + vl_ddn(i,k,m)*vr_ddn(i,l,m)
     $              * dconjg(zd(j,l)*zd(j+3,k)))
     $              + 3*(zd(i+3,k)*zd(i,l)
     $              * dconjg(vl_ddn(j,l,m)*vr_ddn(j,k,m))
     $              + vr_ddn(i,k,m)*vl_ddn(i,l,m)
     $              * dconjg(zd(j,l)*zd(j+3,k))))
            end do
         end do
      end do
      return
      end

      double complex function dd_slr_g(i,j)
c      gluino contributions
      implicit double precision (a-h,o-z)
      double complex zu0,zd0,gm2,gm3
      double complex zd
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      dd_slr_g = (0.d0,0.d0)
      do k=1,6
         do l=1,6
            dd_slr_g = dd_slr_g + 7*g3d*g3d/3.d0*zd(i,k)*zd(i+3,l)
     $           * gm1*gm1*dconjg(zd(j,k)*zd(j+3,l))
     $           * dp0(gm1,gm1,sdm(k),sdm(l))
            dd_slr_g = dd_slr_g - g3d*g3d/18*zd(i,k)*zd(i+3,l)
     $           * dconjg(6*zd(j,k)*zd(j+3,l) + 11*zd(j+3,k)*zd(j,l))
     $           * dp1(gm1,gm1,sdm(k),sdm(l))
         end do
      end do
      return
      end

      double complex function dd_slr(i,j)
c     Full A^S_LR formfactor
      implicit double precision (a-h,o-z)
      double complex dd_slr_hg,dd_slr_c,dd_slr_n,dd_slr_ng,dd_slr_g
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      dd_slr = (0.d0,0.d0)
      if (ih.eq.1)  dd_slr = dd_slr + dd_slr_hg(i,j) 
      if (ic.eq.1)  dd_slr = dd_slr + dd_slr_c(i,j) 
      if (in.eq.1)  dd_slr = dd_slr + dd_slr_n(i,j) 
      if ((in.eq.1).and.(ig.eq.1)) dd_slr = dd_slr + dd_slr_ng(i,j) 
      if (ig.eq.1)  dd_slr = dd_slr + dd_slr_g(i,j) 
      dd_slr = - dd_slr/16/pi/pi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^T_L                                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dd_tl_c(i,j)
c     chargino contributions
      implicit double precision (a-h,o-z)
      double complex vl_duc,vr_duc,zpos,zneg,zu,zd
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dd_tl_c = (0.d0,0.d0)
      do m=1,2
         do n=1,2
            do k=1,6
               do l=1,6
                  dd_tl_c = dd_tl_c 
     $                 - vl_duc(i,k,m)*vl_duc(i,l,n)
     $                 * dconjg(vr_duc(j,k,n)*vr_duc(j,l,m))/16
     $                 * fcm(m)*fcm(n)*dp0(fcm(m),fcm(n),sum(k),sum(l))
               end do
            end do
         end do
      end do
      return
      end

      double complex function dd_tl_n(i,j)
c     neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_ddn,vr_ddn,zn,zu,zd
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dd_tl_n = (0.d0,0.d0)
      do m=1,4
         do n=1,4
            do k=1,6
               do l=1,6
                  dd_tl_n = dd_tl_n 
     $                 - vl_ddn(i,k,m)*vl_ddn(i,l,n)
     $                 * dconjg(vr_ddn(j,k,n)*vr_ddn(j,l,m))/16
     $                 * fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sdm(k),sdm(l))
                  dd_tl_n = dd_tl_n 
     $                 + vl_ddn(i,k,m)*vl_ddn(i,l,m)
     $                 * dconjg(vr_ddn(j,l,n)*vr_ddn(j,k,n))/16
     $                 * fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sdm(k),sdm(l))
               end do
            end do
         end do
      end do
      return
      end

      double complex function dd_tl_ng(i,j)
c     gluino-neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_ddn,vr_ddn,zn,zu0,zd0,gm2,gm3
      double complex zd
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/neut/fnm(4),zn(4,4)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      dd_tl_ng = (0.d0,0.d0)
      do m=1,4
         do k=1,6
            do l=1,6
               dd_tl_ng = dd_tl_ng - g3d/24*zd(i,l)*vl_ddn(i,k,m)
     $              * dconjg(vr_ddn(j,l,m)*zd(j+3,k))
     $              * gm1*fnm(m)*dp0(gm1,fnm(m),sdm(k),sdm(l))
               dd_tl_ng = dd_tl_ng + g3d/24*gm1*fnm(m)
     $              * dp0(gm1,fnm(m),sdm(k),sdm(l))
     $              * (zd(i,k)*zd(i,l)
     $              * dconjg(vr_ddn(j,k,m)*vr_ddn(j,l,m))
     $              + vl_ddn(i,k,m)*vl_ddn(i,l,m)
     $              * dconjg(zd(j+3,k)*zd(j+3,l)))
            end do
         end do
      end do
      return
      end

      double complex function dd_tl_g(i,j)
c     gluino contributions
      implicit double precision (a-h,o-z)
      double complex zu0,zd0,gm2,gm3
      double complex zd
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      dd_tl_g = (0.d0,0.d0)
      do k=1,6
         do l=1,6
            dd_tl_g = dd_tl_g + g3d*g3d/48*zd(i,k)*zd(i,l)
     $           * dconjg(zd(j+3,k)*zd(j+3,l))
     $           * gm1*gm1*dp0(gm1,gm1,sdm(k),sdm(l))
         end do
      end do
      return
      end

      double complex function dd_tl(i,j)
c     Full A^T_L formfactor
      implicit double precision (a-h,o-z)
      double complex dd_tl_c,dd_tl_n,dd_tl_ng,dd_tl_g
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      dd_tl = (0.d0,0.d0)
      if (ic.eq.1)  dd_tl = dd_tl + dd_tl_c(i,j) 
      if (in.eq.1)  dd_tl = dd_tl + dd_tl_n(i,j) 
      if ((in.eq.1).and.(ig.eq.1)) dd_tl = dd_tl + dd_tl_ng(i,j) 
      if (ig.eq.1)  dd_tl = dd_tl + dd_tl_g(i,j) 
      dd_tl = - dd_tl/16/pi/pi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^T_R                                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dd_tr_c(i,j)
c     chargino contributions
      implicit double precision (a-h,o-z)
      double complex vl_duc,vr_duc,zpos,zneg,zu,zd
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dd_tr_c = (0.d0,0.d0)
      do m=1,2
         do n=1,2
            do k=1,6
               do l=1,6
                  dd_tr_c = dd_tr_c 
     $                 - vr_duc(i,k,m)*vr_duc(i,l,n)
     $                 * dconjg(vl_duc(j,k,n)*vl_duc(j,l,m))/16
     $                 * fcm(m)*fcm(n)*dp0(fcm(m),fcm(n),sum(k),sum(l))
               end do
            end do
         end do
      end do
      return
      end

      double complex function dd_tr_n(i,j)
c     neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_ddn,vr_ddn,zn,zu,zd
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dd_tr_n = (0.d0,0.d0)
      do m=1,4
         do n=1,4
            do k=1,6
               do l=1,6
                  dd_tr_n = dd_tr_n 
     $                 - vr_ddn(i,k,m)*vr_ddn(i,l,n)
     $                 * dconjg(vl_ddn(j,k,n)*vl_ddn(j,l,m))/16
     $                 * fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sdm(k),sdm(l))
                  dd_tr_n = dd_tr_n 
     $                 + vr_ddn(i,k,m)*vr_ddn(i,l,m)
     $                 * dconjg(vl_ddn(j,l,n)*vl_ddn(j,k,n))/16
     $                 * fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sdm(k),sdm(l))
               end do
            end do
         end do
      end do
      return
      end

      double complex function dd_tr_ng(i,j)
c     gluino-neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_ddn,vr_ddn,zn,zu0,zd0,gm2,gm3
      double complex zd
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/neut/fnm(4),zn(4,4)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      dd_tr_ng = (0.d0,0.d0)
      do m=1,4
         do k=1,6
            do l=1,6
               dd_tr_ng = dd_tr_ng - g3d/24*zd(i+3,l)*vr_ddn(i,k,m)
     $              * dconjg(vl_ddn(j,l,m)*zd(j,k))
     $              * gm1*fnm(m)*dp0(gm1,fnm(m),sdm(k),sdm(l))
               dd_tr_ng = dd_tr_ng + g3d/24*gm1*fnm(m)
     $              * dp0(gm1,fnm(m),sdm(k),sdm(l))
     $              * (zd(i+3,k)*zd(i+3,l)
     $              * dconjg(vl_ddn(j,k,m)*vl_ddn(j,l,m))
     $              + vr_ddn(i,k,m)*vr_ddn(i,l,m)
     $              * dconjg(zd(j,k)*zd(j,l)))
            end do
         end do
      end do
      return
      end

      double complex function dd_tr_g(i,j)
c     gluino contributions
      implicit double precision (a-h,o-z)
      double complex zu0,zd0,gm2,gm3
      double complex zd
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      dd_tr_g = (0.d0,0.d0)
      do k=1,6
        do l=1,6
           dd_tr_g = dd_tr_g + g3d*g3d/48*zd(i+3,k)*zd(i+3,l)
     $          * dconjg(zd(j,k)*zd(j,l))
     $          * gm1*gm1*dp0(gm1,gm1,sdm(k),sdm(l))
        end do
      end do
      return
      end
      
      double complex function dd_tr(i,j)
c     Full A^T_R formfactor
      implicit double precision (a-h,o-z)
      double complex dd_tr_c,dd_tr_n,dd_tr_ng,dd_tr_g
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      dd_tr = (0.d0,0.d0)
      if (ic.eq.1)  dd_tr = dd_tr + dd_tr_c(i,j) 
      if (in.eq.1)  dd_tr = dd_tr + dd_tr_n(i,j) 
      if ((in.eq.1).and.(ig.eq.1)) dd_tr = dd_tr + dd_tr_ng(i,j) 
      if (ig.eq.1)  dd_tr = dd_tr + dd_tr_g(i,j) 
      dd_tr = - dd_tr/16/pi/pi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Wilson coefficients generated by the effective neutral       c
c     Higgs Yukawa couplings                                       c
c     General form of the Hamiltonian is:                          c
c     H = A^S_LL H^S_LL + A^S_RR H^S_RR + A^S_LR H^S_LR            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      Scalar-LL,RR,LR formfactor to delta_F=2 mixing (KK,BB)
c     Y_dL^IJk = ysd(I,J,k)^* etc. 

      double complex function dd_sll_yuk(i,j)
      implicit double precision (a-h,o-z)
      double complex yhl,yhr,ysu,ypu,ysd,ypd
      logical init_yukawa_eff
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/yukawa_eff/yhl(3,3,2),yhr(3,3,2),ysu(3,3,2),ypu(3,3,2),
     $     ysd(3,3,2),ypd(3,3,2),init_yukawa_eff
      if (init_yukawa_eff) call yukawa_eff_init
      dd_sll_yuk = (0.d0,0.d0)
      do k=1,2
         dd_sll_yuk = dd_sll_yuk + (ysd(i,j,k)/rm(k))**2
     $        - (ypd(i,j,k)/pm(k))**2
      end do
      dd_sll_yuk = dconjg(dd_sll_yuk)/4
      return
      end

      double complex function dd_srr_yuk(i,j)
      implicit double precision (a-h,o-z)
      double complex yhl,yhr,ysu,ypu,ysd,ypd
      logical init_yukawa_eff
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/yukawa_eff/yhl(3,3,2),yhr(3,3,2),ysu(3,3,2),ypu(3,3,2),
     $     ysd(3,3,2),ypd(3,3,2),init_yukawa_eff
      if (init_yukawa_eff) call yukawa_eff_init
      dd_srr_yuk = (0.d0,0.d0)
      do k=1,2
         dd_srr_yuk = dd_srr_yuk + (ysd(j,i,k)/rm(k))**2
     $        - (ypd(j,i,k)/pm(k))**2
      end do
      dd_srr_yuk = dd_srr_yuk/4
      return
      end
      
      double complex function dd_slr_yuk(i,j)
      implicit double precision (a-h,o-z)
      double complex yhl,yhr,ysu,ypu,ysd,ypd
      logical init_yukawa_eff
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/yukawa_eff/yhl(3,3,2),yhr(3,3,2),ysu(3,3,2),ypu(3,3,2),
     $     ysd(3,3,2),ypd(3,3,2),init_yukawa_eff
      if (init_yukawa_eff) call yukawa_eff_init
      dd_slr_yuk = (0.d0,0.d0)
      do k=1,2
         dd_slr_yuk = dd_slr_yuk 
     $        + ysd(j,i,k)*dconjg(ysd(i,j,k))/rm(k)**2
     $        + ypd(j,i,k)*dconjg(ypd(i,j,k))/pm(k)**2
      end do
      dd_slr_yuk = dd_slr_yuk/2
      return
      end

