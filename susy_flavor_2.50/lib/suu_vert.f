c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor
 
c     FILENAME: SUU_VERT.F
c     Released: 08:06:2013(J.R.)
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains scalar-up q-up q vertex formfactors          c
c                                                                     c
c     Outgoing scalar:      momentum p2                               c
c     Incoming quark :      momentum q                                c
c     Outgoing quark :      momentum p2                               c
c                                                                     c
c                                                                     c
c                            V1______/_____ u_J                       c
c                            /|      \  q                             c
c                   k+p1+p2 / |                                       c 
c                         |/  |                                       c
c                         /~  |                                       c
c             H_i   L2_l /    |                                       c
c             _ _/_ _ _ /    /|\ k                                    c
c             p2 \      \V2   |                                       c
c                        \    | L1_n                                  c
c                    L3_m \|  |                                       c 
c                         ~\  |                                       c
c                      k+p1 \ |                                       c
c                            \|______\______ u_K                      c
c                             V3     /  p1                            c
c                                                                     c
c     General form of the vertex                                      c
c                                                                     c
c     V = V_tree + i(V_L P_L + V_R P_R)                               c
c                                                                     c
c     Momentum arguments in formfactors: p1^2=um(K)^2                 c
c                                        q^2 =um(J)^2                 c
c                                        p2^2^2=hm(i)^2               c
c                                                                     c
c     Other arguments:                                                c
c     dm1,dm2,dm3:     masses of particles circulating in loop        c 
c                      on lines L1,L2,L3 respectively                 c
c                                                                     c
c     vl_i,vr_i: LR couplings in the vertices (complex in general)    c
c                                                                     c
c     form:      complex output array containing formfactor values    c
c                                                                     c
c     CAUTION: use those routines only for J<>K - some diagrams       c
c     contributing only to flavor diagonal Suu vertex (like photon,   c 
c     Z or gluon exchanges) are not implemented here                  c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine suu_vert_ddw(i,j,k,form)
c     down-quark-down-quark-W in loop
      implicit double precision (a-h,o-z)
      double complex tmp(2),form(2),cv
      double complex cz,co,ci
      double complex ckm
      common/num/cz,co,ci,zero,one
      common/fmass/em(3),um(3),dm(3)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/km_mat/ckm(3,3)
      common/vev/v1,v2
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      do l=1,3
         cv = - e2/2/st2*dm(l)*zr(1,i)/v1*ckm(l,j)*dconjg(ckm(l,k))
         call vff_svert(wm,dm(l),dm(l),co,co,tmp)
         do kk=1,2
            form(kk) = form(kk) + cv*tmp(kk)
         end do
      end do
      return
      end
 
      subroutine suu_vert_wwd(i,j,k,form)
c     W-W-down quark in loop
      implicit double precision (a-h,o-z)
      double complex tmp(2),form(2)
      double complex ckm
      common/fmass/em(3),um(3),dm(3)
      common/km_mat/ckm(3,3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      do n=1,3
         call fvv_svert(dm(n),wm,wm,tmp)
         do kk=1,2
            form(kk) = form(kk) + e2*e2/4/st2/st2*cr(i)*ckm(n,j)
     $           *dconjg(ckm(n,k))*tmp(kk)
         end do
      end do
      return
      end
 
      subroutine suu_vert_hwd(i,j,k,form)
c     Charged Higgs-W-down quark in loop
      implicit double precision (a-h,o-z)
      double complex tmp(2),form(2)
      double complex ckm
      double complex yl,yu,yd
      common/yukawa/yl(3),yu(3),yd(3)
      common/fmass/em(3),um(3),dm(3)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/km_mat/ckm(3,3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/mssm_charged_higgs_min_index/mhmin
      do l=mhmin,2
         do n=1,3
            call fsv_svert(dm(n),cm(l),wm, - dconjg(yd(n))*zh(1,l),
     $           yu(j)*zh(2,l),tmp)
            do kk=1,2
               form(kk) = form(kk) + e2/2/sq2/st2*am(i,l)*ckm(n,j)
     $              *dconjg(ckm(n,k))*tmp(kk)
            end do
         end do
      end do
      return
      end
 
      subroutine suu_vert_whd(i,j,k,form)
c     W-charged Higgs-down quark in loop
      implicit double precision (a-h,o-z)
      double complex tmp(2),form(2)
      double complex ckm
      double complex yl,yu,yd
      common/yukawa/yl(3),yu(3),yd(3)
      common/fmass/em(3),um(3),dm(3)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/km_mat/ckm(3,3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/mssm_charged_higgs_min_index/mhmin
      do l=mhmin,2
         do n=1,3
            call fvs_svert(dm(n),wm,cm(l),dconjg(yu(k))*zh(2,l),
     $           - yd(n)*zh(1,l),tmp)
            do kk=1,2
               form(kk) = form(kk) + e2/2/sq2/st2*am(i,l)*ckm(n,j)
     $              *dconjg(ckm(n,k))*tmp(kk)
            end do
         end do
      end do
      return
      end
 
      subroutine suu_vert_ddh(i,j,k,form)
c     Down-quark-down-quark-charged Higgs in loop
      implicit double precision (a-h,o-z)
      double complex tmp(2),form(2),cv
      double complex cz,co,ci
      double complex ckm
      double complex yl,yu,yd
      common/yukawa/yl(3),yu(3),yd(3)
      common/num/cz,co,ci,zero,one
      common/fmass/em(3),um(3),dm(3)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/km_mat/ckm(3,3)
      common/vev/v1,v2
      common/mssm_charged_higgs_min_index/mhmin
      do l=1,3
         cv = - dm(l)*zr(1,i)/v1*ckm(l,j)*dconjg(ckm(l,k))
         do n=mhmin,2
            call sff_svert(cm(n),dm(l),dm(l),co,co, 
     $           - dconjg(yd(l))*zh(1,n),yu(j)*zh(2,n), 
     $           dconjg(yu(k))*zh(2,n), - yd(l)*zh(1,n),tmp)
            do kk=1,2
               form(kk) = form(kk) + cv*tmp(kk)
            end do
         end do
      end do
      return
      end
 
      subroutine suu_vert_hhd(i,j,k,form)
c     Charged Higgs-charged Higgs-down quark in loop
      implicit double precision (a-h,o-z)
      double complex tmp(2),form(2),cv
      double complex ckm
      double complex yl,yu,yd
      common/fmass/em(3),um(3),dm(3)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/km_mat/ckm(3,3)
      common/yukawa/yl(3),yu(3),yd(3)
      common/mssm_charged_higgs_min_index/mhmin
      do l=mhmin,2
         do m=mhmin,2
            cv = v_hhs(l,m,i)
            do n=1,3
               call fss_svert(dm(n),cm(l),cm(m), 
     $              - dconjg(yd(n))*zh(1,l),yu(j)*zh(2,l), 
     $              dconjg(yu(k))*zh(2,m), - yd(n)*zh(1,m),tmp)
               do kk=1,2
                  form(kk) = form(kk) + cv*ckm(n,j)*dconjg(ckm(n,k))
     $                 *tmp(kk)
               end do
            end do
         end do
      end do
      return
      end
 
      subroutine suu_vert_ccd(i,j,k,form)
c     Chargino-chargino-down squark in loop
      implicit double precision (a-h,o-z)
      double complex tmp(2),form(2)
      double complex vl_ccs,vr_ccs,vl_udc0,vr_udc0
      double complex zu,zd,zpos,zneg
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      do l=1,2
         do m=1,2
            do n=1,6
               call sff_svert(sdm(n),fcm(l),fcm(m),vl_ccs(l,m,i),
     $              vr_ccs(l,m,i),vl_udc0(j,n,l),vr_udc0(j,n,l),
     $              dconjg(vr_udc0(k,n,m)),dconjg(vl_udc0(k,n,m)),tmp)
               do kk=1,2
                  form(kk) = form(kk) + tmp(kk)
               end do
            end do
         end do
      end do
      return
      end
 
      subroutine suu_vert_ddc(i,j,k,form)
c     Chargino-down squark-down squark in loop
      implicit double precision (a-h,o-z)
      double complex tmp(2),form(2),cv
      double complex vl_udc0,vr_udc0,v_dds
      double complex zu,zd,zpos,zneg
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      do l=1,6
         do m=1,6
            cv = v_dds(l,m,i)
            do n=1,2
               call fss_svert(fcm(n),sdm(l),sdm(m),vl_udc0(j,l,n),
     $              vr_udc0(j,l,n),dconjg(vr_udc0(k,m,n)),
     $              dconjg(vl_udc0(k,m,n)),tmp)
               do kk=1,2
                  form(kk) = form(kk) + cv*tmp(kk)
               end do
            end do
         end do
      end do
      return
      end
 
      subroutine suu_vert_nnu(i,j,k,form)
c     Up squark-neutralino-neutralino in loop
      implicit double precision (a-h,o-z)
      double complex tmp(2),form(2)
      double complex vl_nns,vr_nns,vl_uun0,vr_uun0
      double complex zu,zd,zn
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      do l=1,4
         do m=1,4
            do n=1,6
               call sff_svert(sum(n),fnm(l),fnm(m),vl_nns(l,m,i),
     $              vr_nns(l,m,i),vl_uun0(j,n,l),vr_uun0(j,n,l),
     $              dconjg(vr_uun0(k,n,m)),dconjg(vl_uun0(k,n,m)),tmp)
               do kk=1,2
                  form(kk) = form(kk) + tmp(kk)
               end do
            end do
         end do
      end do
      return
      end
 
      subroutine suu_vert_uun(i,j,k,form)
c     Neutralino-up squark-up squark in loop
      implicit double precision (a-h,o-z)
      double complex tmp(2),form(2),cv
      double complex vl_uun0,vr_uun0,v_uus
      double complex zu,zd,zn
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      do l=1,6
         do m=1,6
            cv = v_uus(m,l,i)
            do n=1,4
               call fss_svert(fnm(n),sum(l),sum(m),vl_uun0(j,l,n),
     $              vr_uun0(j,l,n),dconjg(vr_uun0(k,m,n)),
     $              dconjg(vl_uun0(k,m,n)),tmp)
               do kk=1,2
                  form(kk) = form(kk) + cv*tmp(kk)
               end do
            end do
         end do
      end do
      return
      end
 
      subroutine suu_vert_uug(i,j,k,form)
c     Gluino-up squark-up squark in loop
      implicit double precision (a-h,o-z)
      double complex tmp(2),form(2)
      double complex v_uus
      double complex zu,zd
      double complex gm2,gm3
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/gmass/gm1,gm2,gm3
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      if (init_alpha_susy) call init_alpha_s_susy
      do l=1,6
         do m=1,6
            call fss_svert(gm1,sum(l),sum(m), - dconjg(zu(j,l)),
     $           dconjg(zu(j+3,l)),zu(k+3,m), - zu(k,m),tmp)
            do kk=1,2
               form(kk) = form(kk) + 8*g3u*g3u/3.d0*v_uus(m,l,i)
     $              *tmp(kk)
            end do
         end do
      end do
      return
      end

      subroutine suu_triangle(i,j,k,form)
c     Full triangle Suu formfactors for on-shell momenta
      implicit double precision (a-h,o-z)
      double complex form(2)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      if (j.eq.k) stop 'Flavor diagonal Suu vertex not implemented!'
      do l=1,2
         form(l) = (0.d0,0.d0)
      end do
      if (ih.eq.1) then
         call suu_vert_wwd(i,j,k,form)
         call suu_vert_ddw(i,j,k,form)
         call suu_vert_hwd(i,j,k,form)
         call suu_vert_whd(i,j,k,form)
         call suu_vert_hhd(i,j,k,form)
         call suu_vert_ddh(i,j,k,form)
      end if
      if (ic.eq.1) then
         call suu_vert_ccd(i,j,k,form)
         call suu_vert_ddc(i,j,k,form)
      end if
      if (in.eq.1) then
         call suu_vert_nnu(i,j,k,form)
         call suu_vert_uun(i,j,k,form)
      end if
      if (ig.eq.1) call suu_vert_uug(i,j,k,form)
      do l=1,2
         form(l) = form(l)/16/pi/pi
      end do
      return
      end

      subroutine suu_vert(i,j,k,form)
c     Full Suu vertex (triangle + self energy) for on-shell momenta
      implicit double precision (a-h,o-z)
      double complex form(2)
      double complex uvl_self,uvr_self,usl_self,usr_self
      double complex uslj,usrj,uslk,usrk,uvlj,uvrj,uvlk,uvrk
      logical sqdiag,sldiag
      common/qmass_pole/ump(3),dmp(3)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/vev/v1,v2
      common/sff_args/sm,fm1,fm2
      if (j.eq.k) stop 'Flavor diagonal Suu vertex not implemented!'
      call set_ckm(0)
      call set_yukawa(0)
      if(sqdiag().or.sldiag()) 
     $     stop 'suu_vert: problem with sfermion diagonalization?'
      call set_2hdm(al,be,am,hm1,hm2,hmc,2) ! sets 2-loop neutral Higgs parameters
      sm = rm(i)                ! external scalar mass
      fm1 = ump(j)              ! incoming external fermion mass
      fm2 = ump(k)              ! outgoing external fermion mass
      call suu_triangle(i,j,k,form)
      fact = - zr(2,i)/v2/(fm1**2 - fm2**2)
      uslj = usl_self(fm1**2,j,k)
      uslk = usl_self(fm2**2,j,k)
      usrj = usr_self(fm1**2,j,k)
      usrk = usr_self(fm2**2,j,k)
      uvlj = uvl_self(fm1**2,j,k)
      uvlk = uvl_self(fm2**2,j,k)
      uvrj = uvr_self(fm1**2,j,k)
      uvrk = uvr_self(fm2**2,j,k)
      form(1) = form(1) + fact*((fm1*uslk + fm2*usrk + fm2*(fm1*uvlk +
     $     fm2*uvrk))*fm1 - (fm2*uslj + fm1*usrj + fm1*(fm1*uvlj + fm2
     $     *uvrj))*fm2)
      form(2) = form(2) + fact*((fm2*uslk + fm1*usrk + fm2*(fm2*uvlk +
     $     fm1*uvrk))*fm1 - (fm1*uslj + fm2*usrj + fm1*(fm2*uvlj + fm1
     $     *uvrj))*fm2)
      call set_2hdm(al,be,am,hm1,hm2,hmc,1) ! restores tree-level neutral Higgs parameters
      call set_ckm(1)
      call set_yukawa(1)
      if(sqdiag().or.sldiag()) 
     $     stop 'suu_vert: problem with sfermion diagonalization?'
      return
      end
 
