c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor
 
c     FILENAME: U_SELF.F
c     Released: 1: 4:1994 (P.Ch.)
c     Revised:  08:06:2013 (J.R.)
c     Cleaned, formatted and adapted for SUSY_FLAVOR

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expressions for up quarks self-energy          c
c     function and its renormalization.                                 c
c                                                                       c
c     The definition of the self energy as                              c
c     follows (arguments are s=p^2,i,j):                                c
c                                                                       c
c       p      ____    p                                                c
c             |    |          = i (uvl_sig G(p) P_L + uvr_sig G(p) P_R  c
c       ~~~~~~|____|~~~~~~       + usl_sig P_L + usr_sig P_R)           c
c      l_i             l_j                                              c
c                                                                       c
c                                                                       c
c     i and j are the flavors of incoming and outgoing quark            c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Left vector self-energy (proportional to p-slash P_L)             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Photon and gluon loops (also some constant and finite pieces of Z,
c     W loops) are included in overall QED/QCD factors on external
c     lines, so skipped here. They anyway vanish for flavor-violating
c     self-energies, i.e. for i<>j

c     flavor conserving contributions:

      double complex function uvl_sig_z(s,i,j)
c     Up quark-Z0 in loop
      implicit double precision (a-h,o-z)
      double complex b1
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/fmass/em(3),um(3),dm(3)
      uvl_sig_z = (0.d0,0.d0)
      if (i.ne.j) return
      uvl_sig_z = - e2/2/sct2*(1 - 8*st2/3.d0 + 16*st2*st2/9.d0)
     $     *b1(s,um(i),zm)
      return
      end

      double complex function uvl_sig_s(s,i,j)
c     Up quark + scalar in loop
      implicit double precision (a-h,o-z)
      double complex b1
      double complex yl,yd,yu
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/yukawa/yl(3),yd(3),yu(3)
      uvl_sig_s = (0.d0,0.d0)
      if (i.ne.j) return
      do k=1,2
         uvl_sig_s = uvl_sig_s 
     $        - (abs(yu(i))*zr(2,k))**2/2*b1(s,um(i),rm(k))
      end do
      return
      end

      double complex function uvl_sig_p(s,i,j)
c     Up quark + pseudoscalar in loop
      implicit double precision (a-h,o-z)
      double complex b1
      double complex yl,yd,yu
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/yukawa/yl(3),yd(3),yu(3)
      uvl_sig_p = (0.d0,0.d0)
      if (i.ne.j) return
      do k=1,2
         uvl_sig_p = uvl_sig_p 
     $        - (abs(yu(i))*zh(2,k))**2/2*b1(s,um(i),pm(k))
      end do
      return
      end

c     flavor violating contributions:

      double complex function uvl_sig_w(s,i,j)
c     Down quark-W in loop
      implicit double precision (a-h,o-z)
      double complex b1,ckm
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/fmass/em(3),um(3),dm(3)
      common/km_mat/ckm(3,3)
      uvl_sig_w = (0.d0,0.d0)
      do k=1,3
         uvl_sig_w = uvl_sig_w
     $        - e2/st2*ckm(k,i)*dconjg(ckm(k,j))*b1(s,dm(k),wm)
      end do
      return
      end

      double complex function uvl_sig_h(s,i,j)
c     Down quark+charged Higgs in loop
      implicit double precision (a-h,o-z)
      double complex b1,ckm
      double complex yl,yd,yu
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/km_mat/ckm(3,3)
      common/yukawa/yl(3),yd(3),yu(3)
      uvl_sig_h = (0.d0,0.d0)
      do l=1,3
         do k=1,2
            uvl_sig_h = uvl_sig_h - (zh(1,k)*abs(yd(l)))**2*ckm(l,i)
     $           *dconjg(ckm(l,j))*b1(s,dm(l),cm(k))
         end do
      end do
      return
      end

      double complex function uvl_sig_n(s,i,j)
c     Neutralino-up squark in loop
      implicit double precision (a-h,o-z)
      double complex b1
      double complex zu,zd,zn
      double complex vl_uun0
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      uvl_sig_n = (0.d0,0.d0)
      do k=1,6
         do l=1,4
            uvl_sig_n = uvl_sig_n - vl_uun0(i,k,l)
     $           *dconjg(vl_uun0(j,k,l))*b1(s,fnm(l),sum(k))
         end do
      end do
      return
      end

      double complex function uvl_sig_c(s,i,j)
c     Chargino-down squark in loop
      implicit double precision (a-h,o-z)
      double complex b1
      double complex zu,zd,zpos,zneg
      double complex vl_udc0
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      uvl_sig_c = (0.d0,0.d0)
      do k=1,6
         do l=1,2
            uvl_sig_c = uvl_sig_c - vl_udc0(i,k,l)
     $           *dconjg(vl_udc0(j,k,l))*b1(s,fcm(l),sdm(k))
         end do
      end do
      return
      end

      double complex function uvl_sig_gl(s,i,j)
c     gluino-up squark in loop
      implicit double precision (a-h,o-z)
      double complex b1
      double complex zu,zd
      double complex gm2,gm3
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      uvl_sig_gl = (0.d0,0.d0)
      do k=1,6
         uvl_sig_gl = uvl_sig_gl - 8*g3u*g3u*zu(j,k)*dconjg(zu(i,k))
     $        *b1(s,gm1,sum(k))/3.d0
      end do
      return
      end

      double complex function uvl_self(s,i,j)
c     Full bare up quark self-energy, left vector part
      implicit double precision (a-h,o-z)
      double complex uvl_sig_z,uvl_sig_w,uvl_sig_h,uvl_sig_s,uvl_sig_p,
     $     uvl_sig_n,uvl_sig_c,uvl_sig_gl
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      uvl_self = (0.d0,0.d0)
      if (ih.eq.1) then
         uvl_self = uvl_self + uvl_sig_z(s,i,j)
         uvl_self = uvl_self + uvl_sig_w(s,i,j)
         uvl_self = uvl_self + uvl_sig_h(s,i,j)
         uvl_self = uvl_self + uvl_sig_s(s,i,j)
         uvl_self = uvl_self + uvl_sig_p(s,i,j)
      end if
      if (in.eq.1) then
         uvl_self = uvl_self + uvl_sig_n(s,i,j)
      end if
      if (ic.eq.1) then
         uvl_self = uvl_self + uvl_sig_c(s,i,j)
      end if
      if (ig.eq.1) then
         uvl_self = uvl_self + uvl_sig_gl(s,i,j)
      end if
      uvl_self = uvl_self/16/pi/pi
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Right vector self-energy (proportional to p-slash P_R)            c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     flavor conserving contributions:

      double complex function uvr_sig_z(s,i,j)
c     Up quark-Z0 in loop
      implicit double precision (a-h,o-z)
      double complex b1
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/fmass/em(3),um(3),dm(3)
      uvr_sig_z = (0.d0,0.d0)
      if (i.ne.j) return
      uvr_sig_z = - 8*e2*st2/ct2/9.d0*b1(s,um(i),zm)
      return
      end

      double complex function uvr_sig_s(s,i,j)
c     Up quark + scalar in loop
      implicit double precision (a-h,o-z)
      double complex b1
      double complex yl,yd,yu
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/yukawa/yl(3),yd(3),yu(3)
      uvr_sig_s = (0.d0,0.d0)
      if (i.ne.j) return
      do k=1,2
         uvr_sig_s = uvr_sig_s 
     $        - (abs(yu(i))*zr(2,k))**2/2*b1(s,um(i),rm(k))
      end do
      return
      end

      double complex function uvr_sig_p(s,i,j)
c     Up quark + pseudoscalar in loop
      implicit double precision (a-h,o-z)
      double complex b1
      double complex yl,yd,yu
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/yukawa/yl(3),yd(3),yu(3)
      uvr_sig_p = (0.d0,0.d0)
      if (i.ne.j) return
      do k=1,2
         uvr_sig_p = uvr_sig_p 
     $        - (abs(yu(i))*zh(2,k))**2/2*b1(s,um(i),pm(k))
      end do
      return
      end

c     flavor violating contributions

      double complex function uvr_sig_h(s,i,j)
c     Down quark+charged Higgs in loop
      implicit double precision (a-h,o-z)
      double complex b1,ckm
      double complex yl,yd,yu
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/km_mat/ckm(3,3)
      common/yukawa/yl(3),yd(3),yu(3)
      uvr_sig_h = (0.d0,0.d0)
      do l=1,3
         do k=1,2
            uvr_sig_h = uvr_sig_h - dconjg(yu(j))*yu(i)*zh(2,k)**2
     $           *ckm(l,i)*dconjg(ckm(l,j))*b1(s,dm(l),cm(k))
         end do
      end do
      return
      end

      double complex function uvr_sig_n(s,i,j)
c     Neutralino-up squark in loop
      implicit double precision (a-h,o-z)
      double complex b1
      double complex zu,zd,zn
      double complex vr_uun0
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      uvr_sig_n = (0.d0,0.d0)
      do k=1,6
         do l=1,4
            uvr_sig_n = uvr_sig_n - vr_uun0(i,k,l)
     $           *dconjg(vr_uun0(j,k,l))*b1(s,fnm(l),sum(k))
         end do
      end do
      return
      end

      double complex function uvr_sig_c(s,i,j)
c     Chargino-down squark in loop
      implicit double precision (a-h,o-z)
      double complex b1
      double complex zu,zd,zpos,zneg
      double complex vr_udc0
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      uvr_sig_c = (0.d0,0.d0)
      do k=1,6
         do l=1,2
            uvr_sig_c = uvr_sig_c - vr_udc0(i,k,l)
     $           *dconjg(vr_udc0(j,k,l))*b1(s,fcm(l),sdm(k))
         end do
      end do
      return
      end

      double complex function uvr_sig_gl(s,i,j)
c     gluino-up squark in loop
      implicit double precision (a-h,o-z)
      double complex b1
      double complex zu,zd
      double complex gm2,gm3
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      uvr_sig_gl = (0.d0,0.d0)
      do k=1,6
         uvr_sig_gl = uvr_sig_gl - 8*g3u*g3u*zu(j+3,k)*dconjg(zu(i+3,k))
     $        *b1(s,gm1,sum(k))/3.d0
      end do
      return
      end

      double complex function uvr_self(s,i,j)
c     Full bare up quark self-energy, right vector part
      implicit double precision (a-h,o-z)
      double complex uvr_sig_z,uvr_sig_h,uvr_sig_s,uvr_sig_p,
     $     uvr_sig_n,uvr_sig_c,uvr_sig_gl
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      uvr_self = (0.d0,0.d0)
      if (ih.eq.1) then
         uvr_self = uvr_self + uvr_sig_z(s,i,j)
         uvr_self = uvr_self + uvr_sig_h(s,i,j)
         uvr_self = uvr_self + uvr_sig_s(s,i,j)
         uvr_self = uvr_self + uvr_sig_p(s,i,j)
      end if
      if (in.eq.1) then
         uvr_self = uvr_self + uvr_sig_n(s,i,j)
      end if
      if (ic.eq.1) then
         uvr_self = uvr_self + uvr_sig_c(s,i,j)
      end if
      if (ig.eq.1) then
         uvr_self = uvr_self + uvr_sig_gl(s,i,j)
      end if
      uvr_self = uvr_self/16/pi/pi
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Left scalar self-energy  (proportional to P_L)                    c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     flavor conserving contributions

      double complex function usl_sig_z(s,i,j)
c     Up quark-Z0 in loop
      implicit double precision (a-h,o-z)
      double complex b0
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/fmass/em(3),um(3),dm(3)
      usl_sig_z = (0.d0,0.d0)
      if (i.ne.j) return
      usl_sig_z = 4*e2/ct2/9.d0*um(i)*(3 - 4*st2)*b0(s,um(i),zm)
      return
      end

      double complex function usl_sig_s(s,i,j)
c     Up quark + scalar in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex yl,yd,yu
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/yukawa/yl(3),yd(3),yu(3)
      usl_sig_s = (0.d0,0.d0)
      if (i.ne.j) return
      do k=1,2
         usl_sig_s = usl_sig_s + um(i)/2*(abs(yu(i))*zr(2,k))**2
     $        *b0(s,um(i),rm(k))
      end do
      return
      end

      double complex function usl_sig_p(s,i,j)
c     Up quark + pseudoscalar in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex yl,yd,yu
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/yukawa/yl(3),yd(3),yu(3)
      usl_sig_p = (0.d0,0.d0)
      if (i.ne.j) return
      do k=1,2
         usl_sig_p = usl_sig_p - um(i)/2*(abs(yu(i))*zh(2,k))**2
     $        *b0(s,um(i),pm(k))
      end do
      return
      end

c     flavor violating contributions

      double complex function usl_sig_h(s,i,j)
c     Down quark + charged Higgs in loop
      implicit double precision (a-h,o-z)
      double complex b0,ckm
      double complex yl,yd,yu
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/km_mat/ckm(3,3)
      common/yukawa/yl(3),yd(3),yu(3)
      usl_sig_h = (0.d0,0.d0)
      do l=1,3
         do k=1,2
            usl_sig_h = usl_sig_h - zh(2,k)*zh(1,k)*dm(l)*yu(i)*yd(l)
     $           *ckm(l,i)*dconjg(ckm(l,j))*b0(s,dm(l),cm(k))
         end do
      end do
      return
      end

      double complex function usl_sig_n(s,i,j)
c     Neutralino-down squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd,zn
      double complex vl_uun0,vr_uun0
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      usl_sig_n = (0.d0,0.d0)
      do k=1,6
         do l=1,4
            usl_sig_n = usl_sig_n + fnm(l)*vl_uun0(i,k,l)
     $           *dconjg(vr_uun0(j,k,l))*b0(s,fnm(l),sum(k))
         end do
      end do
      return
      end

      double complex function usl_sig_c(s,i,j)
c     Chargino-down squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd,zpos,zneg
      double complex vl_udc0,vr_udc0
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      usl_sig_c = (0.d0,0.d0)
      do k=1,6
         do l=1,2
            usl_sig_c = usl_sig_c + fcm(l)*vl_udc0(i,k,l)
     $           *dconjg(vr_udc0(j,k,l))*b0(s,fcm(l),sdm(k))
         end do
      end do
      return
      end

      double complex function usl_sig_gl(s,i,j)
c     gluino-up squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd
      double complex gm2,gm3
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      usl_sig_gl = (0.d0,0.d0)
      do k=1,6
         usl_sig_gl = usl_sig_gl - 8*g3u*g3u*zu(j+3,k)*dconjg(zu(i,k))
     $        *gm1*b0(s,gm1,sum(k))/3.d0
      end do
      return
      end

      double complex function usl_self(s,i,j)
c     Full bare up quark self-energy, scalar part
      implicit double precision (a-h,o-z)
      double complex usl_sig_z,usl_sig_h,usl_sig_s,usl_sig_p,usl_sig_n,
     $     usl_sig_c,usl_sig_gl
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      usl_self = (0.d0,0.d0)
      if (ih.eq.1) then
         usl_self = usl_self + usl_sig_z(s,i,j)
         usl_self = usl_self + usl_sig_h(s,i,j)
         usl_self = usl_self + usl_sig_s(s,i,j)
         usl_self = usl_self + usl_sig_p(s,i,j)
      end if
      if (in.eq.1) then
         usl_self = usl_self + usl_sig_n(s,i,j)
      end if
      if (ic.eq.1) then
         usl_self = usl_self + usl_sig_c(s,i,j)
      end if
      if (ig.eq.1) then
         usl_self = usl_self + usl_sig_gl(s,i,j)
      end if
      usl_self = usl_self/16/pi/pi
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Scalar right self-energy (proportional to P_R)                    c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     flavor conserving contributions

      double complex function usr_sig_z(s,i,j)
c     Up quark-Z0 in loop
      implicit double precision (a-h,o-z)
      double complex b0
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/fmass/em(3),um(3),dm(3)
      usr_sig_z = (0.d0,0.d0)
      if (i.ne.j) return
      usr_sig_z = 4*e2/ct2/9.d0*um(i)*(3 - 4*st2)*b0(s,um(i),zm)
      return
      end

      double complex function usr_sig_s(s,i,j)
c     Up quark + scalar in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex yl,yd,yu
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/yukawa/yl(3),yd(3),yu(3)
      usr_sig_s = (0.d0,0.d0)
      if (i.ne.j) return
      do k=1,2
         usr_sig_s = usr_sig_s + um(i)/2*(abs(yu(i))*zr(2,k))**2
     $        *b0(s,um(i),rm(k))
      end do
      return
      end

      double complex function usr_sig_p(s,i,j)
c     Up quark + pseudoscalar in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex yl,yd,yu
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/yukawa/yl(3),yd(3),yu(3)
      usr_sig_p = (0.d0,0.d0)
      if (i.ne.j) return
      do k=1,2
         usr_sig_p = usr_sig_p - um(i)/2*(abs(yu(i))*zh(2,k))**2
     $        *b0(s,um(i),pm(k))
      end do
      return
      end

c     flavor violating contributions

      double complex function usr_sig_h(s,i,j)
c     Down quark + charged Higgs in loop
      implicit double precision (a-h,o-z)
      double complex b0,ckm
      double complex yl,yd,yu
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/km_mat/ckm(3,3)
      common/yukawa/yl(3),yd(3),yu(3)
      usr_sig_h = (0.d0,0.d0)
      do l=1,3
         do k=1,2
            usr_sig_h = usr_sig_h - zh(2,k)*zh(1,k)*dm(l)*dconjg(yu(j)
     $           *yd(l))*ckm(l,i)*dconjg(ckm(l,j))*b0(s,dm(l),cm(k))
         end do
      end do
      return
      end

      double complex function usr_sig_n(s,i,j)
c     Neutralino-down squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd,zn
      double complex vl_uun0,vr_uun0
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      usr_sig_n = (0.d0,0.d0)
      do k=1,6
         do l=1,4
            usr_sig_n = usr_sig_n + fnm(l)*vr_uun0(i,k,l)
     $           *dconjg(vl_uun0(j,k,l))*b0(s,fnm(l),sum(k))
         end do
      end do
      return
      end

      double complex function usr_sig_c(s,i,j)
c     Chargino-down squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd,zpos,zneg
      double complex vl_udc0,vr_udc0
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      usr_sig_c = (0.d0,0.d0)
      do k=1,6
         do l=1,2
            usr_sig_c = usr_sig_c + fcm(l)*vr_udc0(i,k,l)
     $           *dconjg(vl_udc0(j,k,l))*b0(s,fcm(l),sdm(k))
         end do
      end do
      return
      end

      double complex function usr_sig_gl(s,i,j)
c     gluino-up squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd
      double complex gm2,gm3
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      usr_sig_gl = (0.d0,0.d0)
      do k=1,6
         usr_sig_gl = usr_sig_gl - 8*g3u*g3u*zu(j,k)*dconjg(zu(i+3,k))
     $        *gm1*b0(s,gm1,sum(k))/3.d0
      end do
      return
      end

      double complex function usr_self(s,i,j)
c     Full bare up quark self-energy, scalar part
      implicit double precision (a-h,o-z)
      double complex usr_sig_z,usr_sig_h,usr_sig_s,usr_sig_p,usr_sig_n,
     $     usr_sig_c,usr_sig_gl
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      usr_self = (0.d0,0.d0)
      if (ih.eq.1) then
         usr_self = usr_self + usr_sig_z(s,i,j)
         usr_self = usr_self + usr_sig_h(s,i,j)
         usr_self = usr_self + usr_sig_s(s,i,j)
         usr_self = usr_self + usr_sig_p(s,i,j)
      end if
      if (in.eq.1) then
         usr_self = usr_self + usr_sig_n(s,i,j)
      end if
      if (ic.eq.1) then
         usr_self = usr_self + usr_sig_c(s,i,j)
      end if
      if (ig.eq.1) then
         usr_self = usr_self + usr_sig_gl(s,i,j)
      end if
      usr_self = usr_self/16/pi/pi
      return
      end

