c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor
      
c     FILENAME: EDM_Q.F
c     Released: 25.03.1998 (J.R.)
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expressions for EDM of d and u quarks      c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     contributions to d quark EDM 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      double precision function edm_d_c(i)
c     chargino - up squark contributions
      implicit double precision (a-h,o-z)
      double complex vl_duc,vr_duc,zpos,zneg,zu,zd
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      edm_d_c= 0.d0
      do k=1,6
        do l=1,2
          edm_d_c  = edm_d_c 
     $        + fcm(l)*dimag(vl_duc(i,k,l)*dconjg(vr_duc(i,k,l)))
     $        *(cp11(fcm(l),sum(k)) + cp12(sum(k),fcm(l))/3)
        end do
      end do
      edm_d_c = e/8/pi/pi*edm_d_c
      return
      end
      
      double precision function edm_d_n(i)
c     neutralino-down squark contributions
      implicit double precision (a-h,o-z)
      double complex vl_ddn,vr_ddn,zn,zu,zd
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      edm_d_n= 0.d0
      do k=1,6
         do l=1,4
            edm_d_n = edm_d_n + fnm(l)*cp12(sdm(k),fnm(l))
     $           * dimag(vl_ddn(i,k,l)*dconjg(vr_ddn(i,k,l))) 
         end do
      end do
      edm_d_n = - e/48/pi/pi*edm_d_n
      return
      end
      
      double precision function edm_d_g(i)
c     gluino-down squark contributions
      implicit double precision (a-h,o-z)
      double complex zu0,zd0,gm2,gm3
      double complex zd
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      edm_d_g= 0.d0
      do k=1,6
        edm_d_g = edm_d_g 
     $      + dimag(zd(i,k)*dconjg(zd(i+3,k)))*cp12(sdm(k),gm1)
      end do
      edm_d_g = e*g3d*g3d/18/pi/pi*gm1*edm_d_g
      return
      end
      
      double precision function edm_d(i)
c     Full down quark EDM
      implicit double precision (a-h,o-z)
      common/debug_4q/ih,ic,in,ig
      edm_d = 0.d0
      if (ic.eq.1) edm_d = edm_d + edm_d_c(i)
      if (in.eq.1) edm_d = edm_d + edm_d_n(i)
      if (ig.eq.1) edm_d = edm_d + edm_d_g(i)
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     contributions to u quark EDM 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      double precision function edm_u_c(i)
c     chargino-down squark contributions
      implicit double precision (a-h,o-z)
      double complex vl_udc,vr_udc,zpos,zneg,zu,zd
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      edm_u_c= 0.d0
      do k=1,6
        do l=1,2
          edm_u_c  = edm_u_c 
     $        + fcm(l)*dimag(vl_udc(i,k,l)*dconjg(vr_udc(i,k,l)))
     $        * (cp11(fcm(l),sdm(k)) + cp12(sdm(k),fcm(l))/6)
        end do
      end do
      edm_u_c = - e/8/pi/pi*edm_u_c 
      return
      end
      
      double precision function edm_u_n(i)
c     neutralino-up squark contributions
      implicit double precision (a-h,o-z)
      double complex vl_uun,vr_uun,zn,zu,zd
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      edm_u_n= 0.d0
      do k=1,6
        do l=1,4
          edm_u_n = edm_u_n + fnm(l)*cp12(sum(k),fnm(l))
     $        * dimag(vl_uun(i,k,l)*dconjg(vr_uun(i,k,l)))
        end do
      end do
      edm_u_n = e/24/pi/pi*edm_u_n 
      return
      end
      
      double precision function edm_u_g(i)
c     gluino-up squark contributions
      implicit double precision (a-h,o-z)
      double complex zu0,zd0,gm2,gm3
      double complex zu
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      edm_u_g= 0.d0
      do k=1,6
        edm_u_g = edm_u_g 
     $      + dimag(zu(i,k)*dconjg(zu(i+3,k)))*cp12(sum(k),gm1)
      end do
      edm_u_g = e*g3u*g3u/9/pi/pi*gm1*edm_u_g 
      return
      end
      
      double precision function edm_u(i)
c     Full up quark EDM
      implicit double precision (a-h,o-z)
      common/debug_4q/ih,ic,in,ig
      edm_u = 0.d0
      if (ic.eq.1) edm_u = edm_u + edm_u_c(i)
      if (in.eq.1) edm_u = edm_u + edm_u_n(i)
      if (ig.eq.1) edm_u = edm_u + edm_u_g(i)
      return
      end

