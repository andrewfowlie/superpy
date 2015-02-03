c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM} 
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor
 
c     FILENAME: UU_MIX.F
c     Released: 25:11:1999 (J.R.,P.Ch.)
c     Revised: 21:02:2001 (J.R.)
c     alpha_s(msusy) in gluino contribution calculation used

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expression for the coefficient of the     c
c     Hamiltonian for the D-Dbar (at the quark level).             c
c     General form of the Hamiltonian is:                          c
c     H = A^V_LL H^V_LL + A^V_RR H^V_RR + A^V_LR H^V_LR            c
c       + A^S_LL H^S_LL + A^S_RR H^S_RR + A^S_LR H^S_LR            c
c       + A^T_L H^T_L   + A^T_R H^T_R                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^V_LL                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function uu_vll_sm(i,j)
c     Full A^V_LL SM formfactor
      implicit double precision (a-h,o-z)
      double complex ckm
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/fmass/em(3),um(3),dm(3)
      common/km_mat/ckm(3,3)
      uu_vll_sm = (0.d0,0.d0)
      do m=1,3
        do n=1,3
          uu_vll_sm = uu_vll_sm + e2*e2/4/st2/st2
     $        * dconjg(ckm(m,j)*ckm(n,j))*ckm(m,i)*ckm(n,i)
     $        * ((1 + (dm(m)*dm(n)/wm2)**2/4)/2*dp1(wm,wm,dm(m),dm(n)) 
     $        - (dm(m)*dm(n)/wm)**2*dp0(wm,wm,dm(m),dm(n)))
        end do
      end do
      return
      end

      double complex function uu_vll_hg(i,j)
c     Charged Higgs contributions
      implicit double precision (a-h,o-z)
      double complex ckm
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/km_mat/ckm(3,3)
      common/hangle/ca,sa,cb,sb
      uu_vll_hg = (0.d0,0.d0)
      tb2 = sb*sb/cb/cb 
      do m=1,3
        do n=1,3
          a = 0
          do l=1,2
            a = a + l*tb2*(zh(1,l)/cb)**2*dp1(cm(1),cm(l),dm(m),dm(n))
          end do
          uu_vll_hg = uu_vll_hg + e2*e2/4/st2/st2*(dm(m)*dm(n)/wm2)**2
     $        * dconjg(ckm(m,j)*ckm(n,j))*ckm(m,i)*ckm(n,i)
     $        * (a/8 - tb2*wm2*dp0(wm,cm(1),dm(m),dm(n)))
        end do
      end do
      return
      end

      double complex function uu_vll_c(i,j)
c     chargino contributions
      implicit double precision (a-h,o-z)
      double complex vl_udc,zpos,zneg,zu,zd
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      uu_vll_c = (0.d0,0.d0)
      do m=1,2
        do n=1,2
          do k=1,6
            do l=1,6
              uu_vll_c = uu_vll_c + vl_udc(i,k,m)*vl_udc(i,l,n)
     $            * dconjg(vl_udc(j,k,n)*vl_udc(j,l,m))/8
     $            * dp1(fcm(m),fcm(n),sdm(k),sdm(l))
            end do
          end do
        end do
      end do
      return
      end

      double complex function uu_vll_n(i,j)
c     neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_uun,zn,zu,zd
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      uu_vll_n = (0.d0,0.d0)
      do m=1,4
        do n=1,4
          do k=1,6
            do l=1,6
              uu_vll_n = uu_vll_n + vl_uun(i,k,m)*vl_uun(i,l,n)
     $            * dconjg(vl_uun(j,k,n)*vl_uun(j,l,m))/8
     $            * dp1(fnm(m),fnm(n),sum(k),sum(l))
              uu_vll_n = uu_vll_n + vl_uun(i,k,m)*vl_uun(i,l,m)
     $            * dconjg(vl_uun(j,l,n)*vl_uun(j,k,n))/4
     $            * fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sum(k),sum(l))
            end do
          end do
        end do
      end do
      return
      end

      double complex function uu_vll_ng(i,j)
c     gluino-neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_uun,zn,zu0,zd0,gm2,gm3
      double complex zu
      logical init_alpha_susy
      common/alpha_s_susy/g3,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/neut/fnm(4),zn(4,4)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      uu_vll_ng = (0.d0,0.d0)
      do m=1,4
        do k=1,6
          do l=1,6
            uu_vll_ng = uu_vll_ng + g3/6*vl_uun(i,k,m)*zu(j,k)
     $          * dconjg(vl_uun(j,l,m)*zu(i,l))
     $          * dp1(gm1,fnm(m),sum(k),sum(l))
            uu_vll_ng = uu_vll_ng 
     $          + g3/6*(vl_uun(i,k,m)*vl_uun(i,l,m)*zu(j,k)*zu(j,l) 
     $          + dconjg(zu(i,k)*zu(i,l)*vl_uun(j,k,m)*vl_uun(j,l,m)))
     $          * gm1*fnm(m)*dp0(gm1,fnm(m),sum(k),sum(l))
          end do
        end do
      end do
      return
      end

      double complex function uu_vll_g(i,j)
c     gluino contributions
      implicit double precision (a-h,o-z)
      double complex zu0,zd0,gm2,gm3
      double complex zu
      logical init_alpha_susy
      common/alpha_s_susy/g3,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      uu_vll_g = (0.d0,0.d0)
      do k=1,6
        do l=1,6
          uu_vll_g = uu_vll_g + g3*g3/36*zu(j,k)*zu(j,l)
     $        * dconjg(zu(i,k)*zu(i,l))
     $        * (11*dp1(gm1,gm1,sum(k),sum(l))
     $        + 4*gm1*gm1*dp0(gm1,gm1,sum(k),sum(l)))        
         end do
      end do
      return
      end
      
      double complex function uu_vll(i,j)
c     Full A^V_LL formfactor
      implicit double precision (a-h,o-z)
      double complex uu_vll_sm,uu_vll_hg,uu_vll_c,uu_vll_n,
     $    uu_vll_ng,uu_vll_g
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      uu_vll = (0.d0,0.d0)
      if (ih.eq.1) then 
c     SM contribution should be treated separately from the point of 
c     view of QCD evolution factors. We add it here temporarily,
c     improve that later?
c     Perhaps the same should be done for Higgs contribution uu_vll_hg?
         uu_vll = uu_vll + uu_vll_sm(i,j)
         uu_vll = uu_vll + uu_vll_hg(i,j)
      end if 
      if (in.eq.1) uu_vll = uu_vll + uu_vll_n(i,j) 
      if (ic.eq.1) uu_vll = uu_vll + uu_vll_c(i,j) 
      if (ig.eq.1) uu_vll = uu_vll + uu_vll_g(i,j) 
      if ((ig.eq.1).and.(in.eq.1)) uu_vll = uu_vll + uu_vll_ng(i,j) 
      uu_vll = - uu_vll/16/pi/pi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^V_RR                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function uu_vrr_hg(i,j)
c     Higgs and gauge contributions
      implicit double precision (a-h,o-z)
      double complex ckm
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/km_mat/ckm(3,3)
      common/hangle/ca,sa,cb,sb
      uu_vrr_hg = (0.d0,0.d0)
      do m=1,3
        do n=1,3
          do k=1,2
            do l=1,2
              uu_vrr_hg = uu_vrr_hg
     $            + dconjg(ckm(m,j)*ckm(n,j))*ckm(m,i)*ckm(n,i)
     $            * (zh(2,k)*zh(2,l))**2*dp1(cm(k),cm(l),dm(m),dm(n))
            end do
          end do
        end do
      end do
      uu_vrr_hg = (e2/st2*um(i)*um(j)/wm2/sb/sb)**2*uu_vrr_hg/32
      return
      end

      double complex function uu_vrr_c(i,j)
c     chargino contributions
      implicit double precision (a-h,o-z)
      double complex vr_udc,zpos,zneg,zu,zd
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      uu_vrr_c = (0.d0,0.d0)
      do m=1,2
        do n=1,2
          do k=1,6
            do l=1,6
              uu_vrr_c = uu_vrr_c + vr_udc(i,k,m)*vr_udc(i,l,n)
     $            * dconjg(vr_udc(j,k,n)*vr_udc(j,l,m))/8
     $            * dp1(fcm(m),fcm(n),sdm(k),sdm(l))
            end do
          end do
        end do
      end do
      return
      end

      double complex function uu_vrr_n(i,j)
c     neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vr_uun,zn,zu,zd
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      uu_vrr_n = (0.d0,0.d0)
      do m=1,4
        do n=1,4
          do k=1,6
            do l=1,6
              uu_vrr_n = uu_vrr_n + vr_uun(i,k,m)*vr_uun(i,l,n)
     $            * dconjg(vr_uun(j,k,n)*vr_uun(j,l,m))/8
     $            * dp1(fnm(m),fnm(n),sum(k),sum(l))
              uu_vrr_n = uu_vrr_n + vr_uun(i,k,m)*vr_uun(i,l,m)
     $            * dconjg(vr_uun(j,l,n)*vr_uun(j,k,n))/4
     $            * fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sum(k),sum(l))
            end do
          end do
        end do
      end do
      return
      end

      double complex function uu_vrr_ng(i,j)
c     gluino-neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vr_uun,zn,zu0,zd0,gm2,gm3
      double complex zu
      logical init_alpha_susy
      common/alpha_s_susy/g3,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/neut/fnm(4),zn(4,4)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      uu_vrr_ng = (0.d0,0.d0)
      do m=1,4
        do k=1,6
          do l=1,6
            uu_vrr_ng = uu_vrr_ng + g3/6*vr_uun(i,k,m)*zu(j+3,k)
     $          * dconjg(vr_uun(j,l,m)*zu(i+3,l))
     $          * dp1(gm1,fnm(m),sum(k),sum(l))
            uu_vrr_ng = uu_vrr_ng 
     $          + g3/6*(vr_uun(i,k,m)*vr_uun(i,l,m)*zu(j+3,k)*zu(j+3,l) 
     $          + dconjg(zu(i+3,k)*zu(i+3,l)
     $          *vr_uun(j,k,m)*vr_uun(j,l,m)))
     $          * gm1*fnm(m)*dp0(gm1,fnm(m),sum(k),sum(l))
          end do
        end do
      end do
      return
      end

      double complex function uu_vrr_g(i,j)
c      gluino contributions
      implicit double precision (a-h,o-z)
      double complex zu0,zd0,gm2,gm3
      double complex zu
      logical init_alpha_susy
      common/alpha_s_susy/g3,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      uu_vrr_g = (0.d0,0.d0)
      do k=1,6
        do l=1,6
          uu_vrr_g = uu_vrr_g + g3*g3/36*zu(j+3,k)*zu(j+3,l)
     $        * dconjg(zu(i+3,k)*zu(i+3,l))
     $        * (11*dp1(gm1,gm1,sum(k),sum(l))
     $        + 4*gm1*gm1*dp0(gm1,gm1,sum(k),sum(l)))
        end do
      end do
      return
      end

      double complex function uu_vrr(i,j)
c      Full A^V_RR formfactor
      implicit double precision (a-h,o-z)
      double complex uu_vrr_hg,uu_vrr_c,uu_vrr_n,uu_vrr_ng,uu_vrr_g
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      uu_vrr = (0.d0,0.d0)
      if (ih.eq.1) uu_vrr = uu_vrr + uu_vrr_hg(i,j)
      if (in.eq.1) uu_vrr = uu_vrr + uu_vrr_n(i,j) 
      if (ic.eq.1) uu_vrr = uu_vrr + uu_vrr_c(i,j) 
      if (ig.eq.1) uu_vrr = uu_vrr + uu_vrr_g(i,j) 
      if ((ig.eq.1).and.(in.eq.1)) uu_vrr = uu_vrr + uu_vrr_ng(i,j) 
      uu_vrr = - uu_vrr/16/pi/pi
      return
      end
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^V_LR                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function uu_vlr_hg(i,j)
c     Higgs and gauge contributions
      implicit double precision (a-h,o-z)
      double complex ckm
      double complex yl,yu,yd
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/yukawa/yl(3),yu(3),yd(3)
      common/km_mat/ckm(3,3)
      uu_vlr_hg = (0.d0,0.d0)
      do m=1,3
        do n=1,3
          do k=1,2
            do l=1,2
              uu_vlr_hg = uu_vlr_hg
     $            + dconjg(ckm(m,j)*ckm(n,j))*ckm(m,i)*ckm(n,i)
     $            * zh(1,k)*zh(1,l)*zh(2,k)*zh(2,l)
     $            * yd(m)*yd(m)*dp1(cm(k),cm(l),dm(m),dm(n))
            end do
          end do
        end do
      end do
      uu_vlr_hg = yu(i)*yu(j)/4*uu_vlr_hg
      return
      end

      double complex function uu_vlr_c(i,j)
c      chargino contributions
      implicit double precision (a-h,o-z)
      double complex vl_udc,vr_udc,zpos,zneg,zu,zd
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      uu_vlr_c = (0.d0,0.d0)
      do m=1,2
        do n=1,2
          do k=1,6
            do l=1,6
              uu_vlr_c = uu_vlr_c - vl_udc(i,k,m)*vr_udc(i,l,n)
     $            * dconjg(vl_udc(j,k,n)*vr_udc(j,l,m))/2
     $            * fcm(m)*fcm(n)*dp0(fcm(m),fcm(n),sdm(k),sdm(l))
            end do
          end do
        end do
      end do
      return
      end

      double complex function uu_vlr_n(i,j)
c     neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_uun,vr_uun,zn,zu,zd
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      uu_vlr_n = (0.d0,0.d0)
      do m=1,4
        do n=1,4
          do k=1,6
            do l=1,6
              uu_vlr_n = uu_vlr_n - vl_uun(i,k,m)*vr_uun(i,l,m)
     $            * dconjg(vl_uun(j,k,n)*vr_uun(j,l,n))/4
     $            * dp1(fnm(m),fnm(n),sum(k),sum(l))
              uu_vlr_n = uu_vlr_n - vl_uun(i,k,m)*vr_uun(i,l,n)
     $            * dconjg(vl_uun(j,k,n)*vr_uun(j,l,m))/2
     $            * fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sum(k),sum(l))
            end do
          end do
        end do
      end do
      return
      end

      double complex function uu_vlr_ng(i,j)
c      gluino-neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_uun,vr_uun,zn,zu0,zd0,gm2,gm3
      double complex zu
      logical init_alpha_susy
      common/alpha_s_susy/g3,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/neut/fnm(4),zn(4,4)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      uu_vlr_ng = (0.d0,0.d0)
      do m=1,4
        do k=1,6
          do l=1,6
            uu_vlr_ng = uu_vlr_ng + g3/4*dp1(gm1,fnm(m),sum(k),sum(l))
     $        * (zu(j,k)*vr_uun(i,k,m)*dconjg(zu(i,l)*vr_uun(j,l,m)) 
     $        + zu(j+3,k)*vl_uun(i,k,m)*dconjg(zu(i+3,l)*vl_uun(j,l,m)))
     $        - g3/6*gm1*fnm(m)*dp0(gm1,fnm(m),sum(k),sum(l))
     $        * (zu(j+3,k)*vr_uun(i,k,m)*dconjg(zu(i,l)*vl_uun(j,l,m))
     $        + zu(j,k)*vl_uun(i,k,m)*dconjg(zu(i+3,l)*vr_uun(j,l,m)))
            uu_vlr_ng = uu_vlr_ng -g3/12*dp1(gm1,fnm(m),sum(k),sum(l))
     $        *(3*(dconjg(zu(i,l)*zu(i+3,k)*vl_uun(j,k,m)*vr_uun(j,l,m))
     $        + vl_uun(i,k,m)*vr_uun(i,l,m)*zu(j,l)*zu(j+3,k))
     $        + (dconjg(zu(i+3,k)*zu(i,l)*vl_uun(j,l,m)*vr_uun(j,k,m))
     $        + vr_uun(i,k,m)*vl_uun(i,l,m)*zu(j,l)*zu(j+3,k)))
          end do
        end do
      end do
      return
      end

      double complex function uu_vlr_g(i,j)
c      gluino contributions
      implicit double precision (a-h,o-z)
      double complex zu0,zd0,gm2,gm3
      double complex zu
      logical init_alpha_susy
      common/alpha_s_susy/g3,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      uu_vlr_g = (0.d0,0.d0)
      do k=1,6
        do l=1,6
          uu_vlr_g = uu_vlr_g - gm1*gm1*g3*g3/18*zu(j,k)*zu(j+3,l)
     $        * dconjg(zu(i,k)*zu(i+3,l))*dp0(gm1,gm1,sum(k),sum(l))
          uu_vlr_g = uu_vlr_g + 5*g3*g3/36.d0*dconjg(zu(i,k)*zu(i+3,l))
     $        * (3*zu(j+3,k)*zu(j,l) -  2*zu(j,k)*zu(j+3,l))
     $        * dp1(gm1,gm1,sum(k),sum(l))
        end do
      end do
      return
      end

      double complex function uu_vlr(i,j)
c     Full A^V_LR formfactor
      implicit double precision (a-h,o-z)
      double complex uu_vlr_hg,uu_vlr_c,uu_vlr_n,uu_vlr_ng,uu_vlr_g
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      uu_vlr = (0.d0,0.d0)
      if (ih.eq.1) uu_vlr = uu_vlr + uu_vlr_hg(i,j)
      if (in.eq.1) uu_vlr = uu_vlr + uu_vlr_n(i,j) 
      if (ic.eq.1) uu_vlr = uu_vlr + uu_vlr_c(i,j) 
      if (ig.eq.1) uu_vlr = uu_vlr + uu_vlr_g(i,j) 
      if ((ig.eq.1).and.(in.eq.1)) uu_vlr = uu_vlr + uu_vlr_ng(i,j) 
      uu_vlr = - uu_vlr/16/pi/pi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^S_LL                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function uu_sll_hg(i,j)
c     Higgs and gauge contributions
      implicit double precision (a-h,o-z)
      double complex ckm
      double complex yl,yu,yd
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/yukawa/yl(3),yu(3),yd(3)
      common/km_mat/ckm(3,3)
      uu_sll_hg = (0.d0,0.d0)
      do m=1,3
        do n=1,3
          do k=1,2
            do l=1,2
              uu_sll_hg = uu_sll_hg + yd(m)*yd(n)*dm(m)*dm(n)
     $            * dconjg(ckm(m,j)*ckm(n,j))*ckm(m,i)*ckm(n,i)
     $            * zh(1,k)*zh(1,l)*zh(2,k)*zh(2,l)
     $            * dp0(cm(k),cm(l),dm(m),dm(n))
            end do
          end do
        end do
      end do
      uu_sll_hg = yu(j)*yu(j)/2*uu_sll_hg
      return
      end

      double complex function uu_sll_c(i,j)
c      chargino contributions
      implicit double precision (a-h,o-z)
      double complex vl_udc,vr_udc,zpos,zneg,zu,zd
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      uu_sll_c = (0.d0,0.d0)
      do m=1,2
        do n=1,2
          do k=1,6
            do l=1,6
              uu_sll_c = uu_sll_c - vl_udc(i,k,m)*vl_udc(i,l,n)
     $            * dconjg(vr_udc(j,k,n)*vr_udc(j,l,m))/4
     $            * fcm(m)*fcm(n)*dp0(fcm(m),fcm(n),sdm(k),sdm(l))
            end do
          end do
        end do
      end do
      return
      end

      double complex function uu_sll_n(i,j)
c      neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_uun,vr_uun,zn,zu,zd
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      uu_sll_n = (0.d0,0.d0)
      do m=1,4
        do n=1,4
          do k=1,6
            do l=1,6
              uu_sll_n = uu_sll_n - vl_uun(i,k,m)*vl_uun(i,l,n)/4
     $            * dconjg(vr_uun(j,k,n)*vr_uun(j,l,m))
     $            * fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sum(k),sum(l))
              uu_sll_n = uu_sll_n - vl_uun(i,k,m)*vl_uun(i,l,m)/4
     $            * dconjg(vr_uun(j,k,n)*vr_uun(j,l,n))
     $            * fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sum(k),sum(l))
            end do
          end do
        end do
      end do
      return
      end

      double complex function uu_sll_ng(i,j)
c     gluino-neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_uun,vr_uun,zn,zu0,zd0,gm2,gm3
      double complex zu
      logical init_alpha_susy
      common/alpha_s_susy/g3,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/neut/fnm(4),zn(4,4)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      uu_sll_ng = (0.d0,0.d0)
      do m=1,4
        do k=1,6
          do l=1,6
            uu_sll_ng = uu_sll_ng - 7*g3/6.d0*zu(j+3,k)*vl_uun(i,k,m)
     $         * dconjg(vr_uun(j,l,m)*zu(i,l))
     $         * gm1*fnm(m)*dp0(gm1,fnm(m),sum(k),sum(l))
            uu_sll_ng = uu_sll_ng - g3/6*gm1*fnm(m)
     $         * dp0(gm1,fnm(m),sum(k),sum(l))
     $         * (dconjg(zu(i,k)*zu(i,l)*vr_uun(j,k,m)*vr_uun(j,l,m))
     $         + vl_uun(i,k,m)*vl_uun(i,l,m)*zu(j+3,k)*zu(j+3,l)) 
           end do
        end do
      end do
      return
      end

      double complex function uu_sll_g(i,j)
c     gluino contributions
      implicit double precision (a-h,o-z)
      double complex zu0,zd0,gm2,gm3
      double complex zu
      logical init_alpha_susy
      common/alpha_s_susy/g3,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      uu_sll_g = (0.d0,0.d0)
      do k=1,6
        do l=1,6
          uu_sll_g = uu_sll_g + 37*g3*g3/36.d0*dconjg(zu(i,k)*zu(i,l))
     $        * zu(j+3,k)*zu(j+3,l)
     $        * gm1*gm1*dp0(gm1,gm1,sum(k),sum(l))
        end do
      end do
      return
      end

      double complex function uu_sll(i,j)
c     Full A^S_LL formfactor
      implicit double precision (a-h,o-z)
      double complex uu_sll_hg,uu_sll_c,uu_sll_n,uu_sll_ng,uu_sll_g
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      uu_sll = (0.d0,0.d0)
      if (ih.eq.1) uu_sll = uu_sll + uu_sll_hg(i,j)
      if (in.eq.1) uu_sll = uu_sll + uu_sll_n(i,j) 
      if (ic.eq.1) uu_sll = uu_sll + uu_sll_c(i,j) 
      if (ig.eq.1) uu_sll = uu_sll + uu_sll_g(i,j) 
      if ((ig.eq.1).and.(in.eq.1)) uu_sll = uu_sll + uu_sll_ng(i,j) 
      uu_sll = - uu_sll/16/pi/pi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^S_RR                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function uu_srr_hg(i,j)
c     Higgs and gauge contributions
      implicit double precision (a-h,o-z)
      double complex ckm
      double complex yl,yu,yd
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/yukawa/yl(3),yu(3),yd(3)
      common/km_mat/ckm(3,3)
      uu_srr_hg = (0.d0,0.d0)
      do m=1,3
        do n=1,3
          do k=1,2
            do l=1,2
              uu_srr_hg = uu_srr_hg + yd(m)*yd(n)*dm(m)*dm(n)
     $            * dconjg(ckm(m,j)*ckm(n,j))*ckm(m,i)*ckm(n,i)
     $            * zh(1,k)*zh(1,l)*zh(2,k)*zh(2,l)
     $            * dp0(cm(k),cm(l),dm(m),dm(n))
            end do
          end do
        end do
      end do
      uu_srr_hg = yu(i)*yu(i)/2*uu_srr_hg
      return
      end

      double complex function uu_srr_c(i,j)
c     chargino contributions
      implicit double precision (a-h,o-z)
      double complex vl_udc,vr_udc,zpos,zneg,zu,zd
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      uu_srr_c = (0.d0,0.d0)
      do m=1,2
        do n=1,2
          do k=1,6
            do l=1,6
              uu_srr_c = uu_srr_c - vr_udc(i,k,m)*vr_udc(i,l,n)
     $            * dconjg(vl_udc(j,k,n)*vl_udc(j,l,m))/4
     $            * fcm(m)*fcm(n)*dp0(fcm(m),fcm(n),sdm(k),sdm(l))
            end do
          end do
        end do
      end do
      return
      end

      double complex function uu_srr_n(i,j)
c     neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_uun,vr_uun,zn,zu,zd
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      uu_srr_n = (0.d0,0.d0)
      do m=1,4
        do n=1,4
          do k=1,6
            do l=1,6
              uu_srr_n = uu_srr_n - vr_uun(i,k,m)*vr_uun(i,l,n)/4
     $            * dconjg(vl_uun(j,k,n)*vl_uun(j,l,m))
     $            * fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sum(k),sum(l))
              uu_srr_n = uu_srr_n - vr_uun(i,k,m)*vr_uun(i,l,m)/4
     $            * dconjg(vl_uun(j,k,n)*vl_uun(j,l,n))
     $            * fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sum(k),sum(l))
            end do
          end do
        end do
      end do
      return
      end

      double complex function uu_srr_ng(i,j)
c     gluino-neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_uun,vr_uun,zn,zu0,zd0,gm2,gm3
      double complex zu
      logical init_alpha_susy
      common/alpha_s_susy/g3,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/neut/fnm(4),zn(4,4)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      uu_srr_ng = (0.d0,0.d0)
      do m=1,4
        do k=1,6
          do l=1,6
            uu_srr_ng = uu_srr_ng - 7*g3/6.d0*zu(j,k)*vr_uun(i,k,m)
     $         * dconjg(vl_uun(j,l,m)*zu(i+3,l))
     $         * gm1*fnm(m)*dp0(gm1,fnm(m),sum(k),sum(l))
            uu_srr_ng = uu_srr_ng - g3/6*gm1*fnm(m)
     $         * dp0(gm1,fnm(m),sum(k),sum(l))
     $         *(dconjg(zu(i+3,k)*zu(i+3,l)*vl_uun(j,k,m)*vl_uun(j,l,m))
     $         + vr_uun(i,k,m)*vr_uun(i,l,m)*zu(j,k)*zu(j,l))
          end do
        end do
      end do
      return
      end

      double complex function uu_srr_g(i,j)
c      gluino contributions
      implicit double precision (a-h,o-z)
      double complex zu0,zd0,gm2,gm3
      double complex zu
      logical init_alpha_susy
      common/alpha_s_susy/g3,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      uu_srr_g = (0.d0,0.d0)
      do k=1,6
        do l=1,6
         uu_srr_g = uu_srr_g + 37*g3*g3/36.d0
     $          * dconjg(zu(i+3,k)*zu(i+3,l))*zu(j,k)*zu(j,l)
     $          * gm1*gm1*dp0(gm1,gm1,sum(k),sum(l))
        end do
      end do
      return
      end

      double complex function uu_srr(i,j)
c     Full A^S_RR formfactor
      implicit double precision (a-h,o-z)
      double complex uu_srr_hg,uu_srr_c,uu_srr_n,uu_srr_ng,uu_srr_g
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      uu_srr = (0.d0,0.d0)
      if (ih.eq.1) uu_srr = uu_srr + uu_srr_hg(i,j)
      if (in.eq.1) uu_srr = uu_srr + uu_srr_n(i,j) 
      if (ic.eq.1) uu_srr = uu_srr + uu_srr_c(i,j) 
      if (ig.eq.1) uu_srr = uu_srr + uu_srr_g(i,j) 
      if ((ig.eq.1).and.(in.eq.1)) uu_srr = uu_srr + uu_srr_ng(i,j) 
      uu_srr = - uu_srr/16/pi/pi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^S_LR                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function uu_slr_hg(i,j)
c     Higgs and gauge contributions
      implicit double precision (a-h,o-z)
      double complex ckm
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/yukawa/yl(3),yu(3),yd(3)
      common/km_mat/ckm(3,3)
      uu_slr_hg = (0.d0,0.d0)
      do m=1,3
        do n=1,3
          do k=1,2
            uu_slr_hg = uu_slr_hg - e2/2/st2*zh(2,k)*zh(2,k)
     $          * dconjg(ckm(m,j)*ckm(n,j))*ckm(m,i)*ckm(n,i)
     $          * dp1(wm,cm(k),dm(m),dm(n))
            do l=1,2
              uu_slr_hg = uu_slr_hg
     $            + dm(m)*dm(n)*yd(m)*yd(n)*(zh(2,k)*zh(1,l))**2
     $            * dconjg(ckm(m,j)*ckm(n,j))*ckm(m,i)*ckm(n,i)
     $            * dp0(cm(k),cm(l),dm(m),dm(n))
            end do
          end do
        end do
      end do
      uu_slr_hg = yu(i)*yu(j)*uu_slr_hg
      return
      end
      
      double complex function uu_slr_c(i,j)
c      chargino contributions
      implicit double precision (a-h,o-z)
      double complex vl_udc,vr_udc,zpos,zneg,zu,zd
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      uu_slr_c = (0.d0,0.d0)
      do m=1,2
        do n=1,2
          do k=1,6
            do l=1,6
              uu_slr_c = uu_slr_c - vl_udc(i,k,m)*vr_udc(i,l,n)
     4            * dconjg(vl_udc(j,l,m)*vr_udc(j,k,n))/2
     $            * dp1(fcm(m),fcm(n),sdm(k),sdm(l))
            end do
          end do
        end do
      end do
      return
      end

      double complex function uu_slr_n(i,j)
c     neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_uun,vr_uun,zn,zu,zd
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      uu_slr_n = (0.d0,0.d0)
      do m=1,4
        do n=1,4
          do k=1,6
            do l=1,6
              uu_slr_n = uu_slr_n - vl_uun(i,k,m)*vr_uun(i,l,n)
     $            * dconjg(vl_uun(j,l,m)*vr_uun(j,k,n))/2
     $            * dp1(fnm(m),fnm(n),sum(k),sum(l))
              uu_slr_n = uu_slr_n - vl_uun(i,k,m)*vr_uun(i,l,m)
     $            * dconjg(vl_uun(j,l,n)*vr_uun(j,k,n))/2
     $            * dp1(fnm(m),fnm(n),sum(k),sum(l))
            end do
          end do
        end do
      end do
      return
      end

      double complex function uu_slr_ng(i,j)
c     gluino-neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_uun,vr_uun,zn,zu0,zd0,gm2,gm3
      double complex zu
      logical init_alpha_susy
      common/alpha_s_susy/g3,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/neut/fnm(4),zn(4,4)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      uu_slr_ng = (0.d0,0.d0)
      do m=1,4
        do k=1,6
          do l=1,6
            uu_slr_ng = uu_slr_ng + g3/6*dp1(gm1,fnm(m),sum(k),sum(l))
     $        *(zu(j,k)*vr_uun(i,k,m)*dconjg(zu(i,l)*vr_uun(j,l,m)) 
     $        + zu(j+3,k)*vl_uun(i,k,m)*dconjg(zu(i+3,l)*vl_uun(j,l,m)))
     $        - g3*gm1*fnm(m)*dp0(gm1,fnm(m),sum(k),sum(l))
     $        *(zu(j+3,k)*vr_uun(i,k,m)*dconjg(zu(i,l)*vl_uun(j,l,m))
     $        + zu(j,k)*vl_uun(i,k,m)*dconjg(zu(i+3,l)*vr_uun(j,l,m)))
            uu_slr_ng = uu_slr_ng - g3/6*dp1(gm1,fnm(m),sum(k),sum(l))
     $        * ((dconjg(zu(i,l)*zu(i+3,k)*vl_uun(j,k,m)*vr_uun(j,l,m))
     $        + vl_uun(i,k,m)*vr_uun(i,l,m)*zu(j,l)*zu(j+3,k))
     $        + 3*(dconjg(zu(i+3,k)*zu(i,l)*vl_uun(j,l,m)*vr_uun(j,k,m))
     $        + vr_uun(i,k,m)*vl_uun(i,l,m)*zu(j,l)*zu(j+3,k)))
          end do
        end do
      end do
      return
      end

      double complex function uu_slr_g(i,j)
c      gluino contributions
      implicit double precision (a-h,o-z)
      double complex zu0,zd0,gm2,gm3
      double complex zu
      logical init_alpha_susy
      common/alpha_s_susy/g3,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      uu_slr_g = (0.d0,0.d0)
      do k=1,6
        do l=1,6
          uu_slr_g = uu_slr_g + 7*g3*g3/3.d0*zu(j,k)*zu(j+3,l)
     $        * dconjg(zu(i,k)*zu(i+3,l))
     $        * gm1*gm1*dp0(gm1,gm1,sum(k),sum(l))
          uu_slr_g = uu_slr_g - g3*g3/18*dconjg(zu(i,k)*zu(i+3,l))
     $        * (6*zu(j,k)*zu(j+3,l) + 11*zu(j+3,k)*zu(j,l))
     $        * dp1(gm1,gm1,sum(k),sum(l))
        end do
      end do
      return
      end

      double complex function uu_slr(i,j)
c     Full A^S_LR formfactor
      implicit double precision (a-h,o-z)
      double complex uu_slr_hg,uu_slr_c,uu_slr_n,uu_slr_ng,uu_slr_g
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      uu_slr = (0.d0,0.d0)
      if (ih.eq.1) uu_slr = uu_slr + uu_slr_hg(i,j)
      if (in.eq.1) uu_slr = uu_slr + uu_slr_n(i,j) 
      if (ic.eq.1) uu_slr = uu_slr + uu_slr_c(i,j) 
      if (ig.eq.1) uu_slr = uu_slr + uu_slr_g(i,j) 
      if ((ig.eq.1).and.(in.eq.1)) uu_slr = uu_slr + uu_slr_ng(i,j) 
      uu_slr = - uu_slr/16/pi/pi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^T_L                                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function uu_tl_c(i,j)
c     chargino contributions
      implicit double precision (a-h,o-z)
      double complex vl_udc,vr_udc,zpos,zneg,zu,zd
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      uu_tl_c = (0.d0,0.d0)
      do m=1,2
        do n=1,2
          do k=1,6
            do l=1,6
              uu_tl_c = uu_tl_c - vl_udc(i,k,m)*vl_udc(i,l,n)
     $            * dconjg(vr_udc(j,k,n)*vr_udc(j,l,m))/16
     $            * fcm(m)*fcm(n)*dp0(fcm(m),fcm(n),sdm(k),sdm(l))
            end do
          end do
        end do
      end do
      return
      end

      double complex function uu_tl_n(i,j)
c     neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_uun,vr_uun,zn,zu,zd
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      uu_tl_n = (0.d0,0.d0)
      do m=1,4
        do n=1,4
          do k=1,6
            do l=1,6
              uu_tl_n = uu_tl_n - vl_uun(i,k,m)*vl_uun(i,l,n)
     $            * dconjg(vr_uun(j,k,n)*vr_uun(j,l,m))/16
     $            * fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sum(k),sum(l))
              uu_tl_n = uu_tl_n + vl_uun(i,k,m)*vl_uun(i,l,m)
     $            * dconjg(vr_uun(j,l,n)*vr_uun(j,k,n))/16
     $            * fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sum(k),sum(l))
            end do
          end do
        end do
      end do
      return
      end

      double complex function uu_tl_ng(i,j)
c     gluino-neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_uun,vr_uun,zn,zu0,zd0,gm2,gm3
      double complex zu
      logical init_alpha_susy
      common/alpha_s_susy/g3,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/neut/fnm(4),zn(4,4)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      uu_tl_ng = (0.d0,0.d0)
      do m=1,4
        do k=1,6
          do l=1,6
            uu_tl_ng = uu_tl_ng - g3/24*zu(j+3,k)*vl_uun(i,k,m)
     $          * dconjg(vr_uun(j,l,m)*zu(i,l))
     $          * gm1*fnm(m)*dp0(gm1,fnm(m),sum(k),sum(l))
            uu_tl_ng = uu_tl_ng + g3/24*gm1*fnm(m)
     $         * dp0(gm1,fnm(m),sum(k),sum(l))
     $         * (dconjg(zu(i,k)*zu(i,l)*vr_uun(j,k,m)*vr_uun(j,l,m))
     $         + vl_uun(i,k,m)*vl_uun(i,l,m)*zu(j+3,k)*zu(j+3,l))
          end do
        end do
      end do
      return
      end

      double complex function uu_tl_g(i,j)
c     gluino contributions
      implicit double precision (a-h,o-z)
      double complex zu0,zd0,gm2,gm3
      double complex zu
      logical init_alpha_susy
      common/alpha_s_susy/g3,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      uu_tl_g = (0.d0,0.d0)
      do k=1,6
        do l=1,6
          uu_tl_g = uu_tl_g + g3*g3/48*zu(j+3,k)*zu(j+3,l)
     $        * dconjg(zu(i,k)*zu(i,l))
     $        * gm1*gm1*dp0(gm1,gm1,sum(k),sum(l))
        end do
      end do
      return
      end

      double complex function uu_tl(i,j)
c     Full A^T_L formfactor
      implicit double precision (a-h,o-z)
      double complex uu_tl_c,uu_tl_n,uu_tl_ng,uu_tl_g
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      uu_tl = (0.d0,0.d0)
      if (in.eq.1) uu_tl = uu_tl + uu_tl_n(i,j) 
      if (ic.eq.1) uu_tl = uu_tl + uu_tl_c(i,j) 
      if (ig.eq.1) uu_tl = uu_tl + uu_tl_g(i,j) 
      if ((ig.eq.1).and.(in.eq.1)) uu_tl = uu_tl + uu_tl_ng(i,j) 
      uu_tl = - uu_tl/16/pi/pi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^T_R                                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function uu_tr_c(i,j)
c     chargino contributions
      implicit double precision (a-h,o-z)
      double complex vl_udc,vr_udc,zpos,zneg,zu,zd
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      uu_tr_c = (0.d0,0.d0)
      do m=1,2
        do n=1,2
          do k=1,6
            do l=1,6
              uu_tr_c = uu_tr_c - vr_udc(i,k,m)*vr_udc(i,l,n)
     $            * dconjg(vl_udc(j,k,n)*vl_udc(j,l,m))/16
     $            * fcm(m)*fcm(n)*dp0(fcm(m),fcm(n),sdm(k),sdm(l))
            end do
          end do
        end do
      end do
      return
      end

      double complex function uu_tr_n(i,j)
c     neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_uun,vr_uun,zn,zu,zd
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      uu_tr_n = (0.d0,0.d0)
      do m=1,4
        do n=1,4
          do k=1,6
            do l=1,6
              uu_tr_n = uu_tr_n - vr_uun(i,k,m)*vr_uun(i,l,n)
     $            * dconjg(vl_uun(j,k,n)*vl_uun(j,l,m))/16
     $            * fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sum(k),sum(l))
              uu_tr_n = uu_tr_n + vr_uun(i,k,m)*vr_uun(i,l,m)
     $            * dconjg(vl_uun(j,l,n)*vl_uun(j,k,n))/16
     $            * fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sum(k),sum(l))
            end do
          end do
        end do
      end do
      return
      end

      double complex function uu_tr_ng(i,j)
c     gluino-neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_uun,vr_uun,zn,zu0,zd0,gm2,gm3
      double complex zu
      logical init_alpha_susy
      common/alpha_s_susy/g3,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/neut/fnm(4),zn(4,4)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      uu_tr_ng = (0.d0,0.d0)
      do m=1,4
        do k=1,6
          do l=1,6
            uu_tr_ng = uu_tr_ng - g3/24*zu(j,k)*vr_uun(i,k,m)
     $         * dconjg(vl_uun(j,l,m)*zu(i+3,l))
     $         * gm1*fnm(m)*dp0(gm1,fnm(m),sum(k),sum(l))
            uu_tr_ng = uu_tr_ng + g3/24*gm1*fnm(m)
     $         * dp0(gm1,fnm(m),sum(k),sum(l))
     $         *(dconjg(zu(i+3,k)*zu(i+3,l)*vl_uun(j,k,m)*vl_uun(j,l,m))
     $         + vr_uun(i,k,m)*vr_uun(i,l,m)*zu(j,k)*zu(j,l))
          end do
        end do
      end do
      return
      end

      double complex function uu_tr_g(i,j)
c     gluino contributions
      implicit double precision (a-h,o-z)
      double complex zu0,zd0,gm2,gm3
      double complex zu
      logical init_alpha_susy
      common/alpha_s_susy/g3,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      uu_tr_g = (0.d0,0.d0)
      do k=1,6
        do l=1,6
          uu_tr_g = uu_tr_g + g3*g3/48*zu(j,k)*zu(j,l)
     $        * dconjg(zu(i+3,k)*zu(i+3,l))
     $        * gm1*gm1*dp0(gm1,gm1,sum(k),sum(l))
        end do
      end do
      return
      end

      double complex function uu_tr(i,j)
c     Full A^T_R formfactor
      implicit double precision (a-h,o-z)
      double complex uu_tr_c,uu_tr_n,uu_tr_ng,uu_tr_g
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      uu_tr = (0.d0,0.d0)
      if (in.eq.1) uu_tr = uu_tr + uu_tr_n(i,j) 
      if (ic.eq.1) uu_tr = uu_tr + uu_tr_c(i,j) 
      if (ig.eq.1) uu_tr = uu_tr + uu_tr_g(i,j) 
      if ((ig.eq.1).and.(in.eq.1)) uu_tr = uu_tr + uu_tr_ng(i,j) 
      uu_tr = - uu_tr/16/pi/pi
      return
      end






