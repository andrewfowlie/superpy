c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor
 
c     FILENAME: UU_GLUON.F
c     Released: 30:04:2014(J.R.)
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains gluon-up q-up q vertex formfactors           c
c                                                                     c
c     Outgoing gluon:       momentum q                                c
c     Incoming quark :      momentum pi                               c
c     Outgoing quark :      momentum pj                               c
c                                                                     c
c                                                                     c
c                            V3______\_____ gluon^mu                  c
c                            /|      /  q                             c
c                           / |                                       c 
c                       m1|/  |                                       c
c                         /~  |                                       c
c             u_I        /    | m1                                    c
c             _ _\_ _ _ /    /|\                                      c
c             pi /      \V1   |                                       c
c                        \    |                                       c
c                      m2 \|  |                                       c 
c                         ~\  |                                       c
c                           \ |                                       c
c                            \|______\______ u_J                      c
c                             V2     /  pj                            c
c                                                                     c
c     General form of the vertex (g_s extracted from VL, VR!          c
c                                                                     c
c     V = V_tree - g_s sigma^munu (V_L P_L + V_R P_R) q_nu            c
c                                                                     c
c     Momentum arguments in formfactors: q^2=0                        c
c                                        pi^2=um(I)^2                 c
c                                        pj^2=um(J)^2->0              c
c                                                                     c
c                                                                     c
c     vl_i,vr_i: LR couplings in the vertices (complex in general)    c
c                                                                     c
c     form:      complex output array containing formfactor values    c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine gluu_vert_ddw(i,j,form)
c     down-quark-down-quark-W in loop
      implicit double precision (a-h,o-z)
      double complex tmp(2),form(2),cv
      double complex cz,co,ci
      double complex ckm
      common/num/cz,co,ci,zero,one
      common/fmass/em(3),um(3),dm(3)
      common/km_mat/ckm(3,3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      do l=1,3
         cv = - e2/2/st2*ckm(l,i)*dconjg(ckm(l,j))
         call vvf_vert(dm(l),wm,1.d0,co,co,tmp)
         do kk=1,2
            form(kk) = form(kk) + cv*tmp(kk)
         end do
      end do
      return
      end
 
      subroutine gluu_vert_ddh(i,j,form)
c     Down-quark-down-quark-charged Higgs in loop
      implicit double precision (a-h,o-z)
      double complex tmp(2),form(2),cv
      double complex ckm
      double complex yl,yu,yd
      common/yukawa/yl(3),yu(3),yd(3)
      common/fmass/em(3),um(3),dm(3)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/km_mat/ckm(3,3)
      common/mssm_charged_higgs_min_index/mhmin
      do l=1,3
         cv = - ckm(l,i)*dconjg(ckm(l,j))
         do n=mhmin,2
            call vsf_vert(dm(l),cm(n),0.d0,1.d0, 
     $           - dconjg(yd(l))*zh(1,n),yu(i)*zh(2,n), 
     $           dconjg(yu(j))*zh(2,n), - yd(l)*zh(1,n),tmp)
            do kk=1,2
               form(kk) = form(kk) + cv*tmp(kk)
            end do
         end do
      end do
      return
      end
 
      subroutine gluu_vert_ddc(i,j,form)
c     Chargino-down squark-down squark in loop
      implicit double precision (a-h,o-z)
      double complex tmp(2),form(2)
      double complex vl_udc0,vr_udc0
      double complex zu,zd,zpos,zneg
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      do l=1,6
         do n=1,2
            call vsf_vert(fcm(n),sdm(l),0.d0,1.d0,vl_udc0(i,l,n),
     $           vr_udc0(i,l,n),dconjg(vr_udc0(j,l,n)), 
     $           dconjg(vl_udc0(j,l,n)),tmp)
            do kk=1,2
               form(kk) = form(kk) - tmp(kk)
            end do
         end do
      end do
      return
      end

      subroutine gluu_vert_uun(i,j,form)
c     Neutralino-up squark-up squark in loop
      implicit double precision (a-h,o-z)
      double complex tmp(2),form(2)
      double complex vl_uun0,vr_uun0
      double complex zu,zd,zn
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      do l=1,6
         do n=1,4
            call vsf_vert(fnm(n),sum(l),0.d0,1.d0,vl_uun0(i,l,n),
     $           vr_uun0(i,l,n),dconjg(vr_uun0(j,l,n)), 
     $           dconjg(vl_uun0(j,l,n)),tmp)
            do kk=1,2
               form(kk) = form(kk) - tmp(kk)
            end do
         end do
      end do
      return
      end

      subroutine gluu_vert_ug(i,j,form)
c     Gluino-up squark in loop (sum of 2 diagrams)
      implicit double precision (a-h,o-z)
      double complex tmp(2),form(2)
      double complex zu,zd
      double complex gm2,gm3
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/gmass/gm1,gm2,gm3
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      if (init_alpha_susy) call init_alpha_s_susy
      do l=1,6
         call vsf_vert(gm1,sum(l),-3.d0,1/3.d0, - dconjg(zu(i,l)),
     $        dconjg(zu(i+3,l)),zu(j+3,l), - zu(j,l),tmp)
         do kk=1,2
            form(kk) = form(kk) + g3u**2*tmp(kk)
         end do
      end do
      return
      end

      subroutine gluu_vert(i,j,form)
c     Full gluon-uu formfactors for on-shell momenta
      implicit double precision (a-h,o-z)
      double complex form(2)
      logical sqdiag,sldiag
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/qmass_pole/ump(3),dmp(3)
      common/vff_args/fmi
      common/debug_4q/ih,ic,in,ig
      if (j.ge.i) stop 'gluon-uu vertex can be called only for i>j!'
      fmi = ump(i)
      do l=1,2
         form(l) = (0.d0,0.d0)
      end do
      call set_ckm(0)
      call set_yukawa(0)
      if(sqdiag().or.sldiag()) 
     $     stop 'gluu_vert: problem with sfermion diagonalization?'
      if (ih.eq.1) then
         call gluu_vert_ddw(i,j,form)
         call gluu_vert_ddh(i,j,form)
      end if
      if (ic.eq.1) then
         call gluu_vert_ddc(i,j,form)
      end if
      if (in.eq.1) then
         call gluu_vert_uun(i,j,form)
      end if
      if (ig.eq.1) call gluu_vert_ug(i,j,form)
      call set_ckm(1)
      call set_yukawa(1)
      if(sqdiag().or.sldiag()) 
     $     stop 'gluu_vert: problem with sfermion diagonalization?'
      do l=1,2
         form(l) = form(l)/16/pi/pi
      end do
      return
      end

 
