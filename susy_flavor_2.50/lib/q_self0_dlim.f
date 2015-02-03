c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM}
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: Q_SELF0_DLIM.F
c     Released: 11: 8:2011(J.R.)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains initialization routines for the various       c
c     decompositions of the u and d quark self energies, calculated    c
c     at p^2 = 0 and in the decoupling limit only                      c
c     Conventions follow the arXiv:1103.4272 [hep-ph], i.e. comparing  c
c     to q_self0.f one has sig_dd(i,j) ~ dsl_sig(i,j)^*                c
c     and  sig_uu(i,j) ~ usr_sig(i,j)^*                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function eps_dd(i)
c     Diagonal "epsilon_d" terms from d-quark self energies in the
c     decoupling limit
      implicit double precision (a-h,o-z)
      double complex ls,ks,ds,es,us,ws
      double complex gm2,gm3,zu,zd,mu
      double complex tmp,tmp1,tmp2
      double complex ckm
      double complex yl,yu,yd
      logical init_alpha_susy
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      common/hpar/hm1,hm2,hs,mu
      common/yukawa/yl(3),yu(3),yd(3)
      common/ckm/ckm(3,3)
      common/vev/v1,v2
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/debug_4q/ih,ic,in,ig
      eps_dd = (0.d0,0.d0)
c     Neutralino+gluino contribution, terms proportional to Yukawa couplings
      if (in.eq.1) then
         do m=1,6
            eps_dd = eps_dd - e2/st2/2.d0*gm2*mu*abs(zd(i,m))**2
     $           * cp0(abs(mu),abs(gm2),sdm(m))
     $           - e2/ct2*gm3*mu/6.d0*(abs(zd(i,m))**2
     $           + 2*abs(zd(i+3,m))**2)*cp0(abs(mu),abs(gm3),sdm(m))
            do n=1,6
               eps_dd = eps_dd - abs(zd(i,m)*zd(i+3,n))**2*mu
     $              * e2/ct2/9.d0*gm3*cp0(abs(gm3),sdm(m),sdm(n))
            end do
         end do
      end if
      if (ig.eq.1) then
         do m=1,6
            do n=1,6
               eps_dd = eps_dd + abs(zd(i,m)*zd(i+3,n))**2*mu
     $              * 8/3.d0*g3d*g3d*gm1*cp0(gm1,sdm(m),sdm(n))
            end do
         end do
      end if
c     Chargino contribution
      if (ic.eq.1) then
         do m=1,6
            eps_dd = eps_dd - e2/st2*mu*gm2*abs(zd(i,m))**2
     $           * cp0(abs(mu),abs(gm2),sdm(m))
            do n=1,6
               tmp  = (0.d0,0.d0)
               tmp1 = (0.d0,0.d0)
               tmp2 = (0.d0,0.d0)
               do l=1,3
                  tmp  = tmp + yu(l)*dconjg(ckm(l,i)*zu(l+3,m))
                  tmp1 = tmp1 + ckm(l,i)*zu(l,n)
                  do k=1,3
                     tmp2 = tmp2 + zu(k+3,m)*(us(l,k) + v1/v2*ws(l,k))
     $                    * dconjg(zu(l,n))
                  end do
                  tmp2 = tmp2 + v1/v2*zu(l+3,m)*dconjg(mu*yu(l)*zu(l,n))
               end do
               eps_dd = eps_dd - mu*tmp*tmp1*tmp2
     $              * cp0(abs(mu),sum(m),sum(n))
            end do
         end do
      end if
      eps_dd = - eps_dd/16/pi/pi
      return
      end

      double complex function eps_fc()
c     "epsilon_FC" term, derived from chirally enhanced diagonal d-quark
c     self energy, part proportional to V_CKM in the decoupling limit
      implicit double precision (a-h,o-z)
      double complex zu,zd
      double complex ls,ks,ds,es,us,ws
      double complex mu
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/vev/v1,v2
      common/hpar/hm1,hm2,hs,mu
      common/yukawa/yl(3),yu(3),yd(3)
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/fmass/em(3),um(3),dm(3)
      common/debug_4q/ih,ic,in,ig
      eps_fc = (0.d0,0.d0)
c     Chargino contribution only
      if (ic.ne.1) return
      do m=1,6
         do n=1,6
            eps_fc = eps_fc + (v1*ws(3,3) + v1*dconjg(mu*yu(3))
     $           + v2*us(3,3))*abs(zu(6,m)*zu(3,n))**2
     $           * cp0(abs(mu),sum(m),sum(n))
         end do
      end do
      eps_fc = yd(3)/dm(3)*yu(3)/sq2/16/pi/pi*mu*eps_fc
      return
      end

      double complex function sig_dd_yb(i)
c     d and s self energy, part propotional to Yb, in the decoupling limit
      implicit double precision (a-h,o-z)
      double complex mu
      double complex gm2,gm3
      double complex zu,zd
      double complex yl,yu,yd
      logical init_alpha_susy
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      common/vev/v1,v2
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/hpar/hm1,hm2,hs,mu
      common/yukawa/yl(3),yu(3),yd(3)
      common/debug_4q/ih,ic,in,ig
      sig_dd_yb = (0.d0,0.d0)
      if (init_alpha_susy) call init_alpha_s_susy
c     Neutralino+gluino contribution
      do m=1,6
         do n=1,6
            if (in.eq.1) then
               sig_dd_yb = sig_dd_yb + yd(3)*v2/sq2*mu
     $              * dconjg(zd(6,n)*zd(i,m))*zd(3,m)*zd(i+3,n)
     $              * e2/ct2/9.d0*gm3*cp0(abs(gm3),sdm(m),sdm(n))
            end if
            if (ig.eq.1) then
               sig_dd_yb = sig_dd_yb - yd(3)*v2/sq2*mu
     $              * dconjg(zd(6,n)*zd(i,m))*zd(3,m)*zd(i+3,n)
     $              * 8/3.d0*g3d*g3d*gm1*cp0(gm1,sdm(m),sdm(n))
         end if
         end do
      end do
      sig_dd_yb = - sig_dd_yb/16/pi/pi
      return
      end

      double complex function sig_dd_y(i,j)
c     Chirally enhanced diagonal d-quark self energy, terms propotional
c     to down Yukawa, in the decoupling limit
      implicit double precision (a-h,o-z)
      double complex tmp,tmp1,tmp2
      double complex mu
      double complex gm2,gm3
      double complex zu,zd
      double complex ls,ks,ds,es,us,ws
      double complex yl,yu,yd
      double complex ckm
      logical init_alpha_susy
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      common/vev/v1,v2
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/hpar/hm1,hm2,hs,mu
      common/yukawa/yl(3),yu(3),yd(3)
      common/debug_4q/ih,ic,in,ig
      common/ckm/ckm(3,3)
      if (init_alpha_susy) call init_alpha_s_susy
      sig_dd_y = (0.d0,0.d0)
c     Neutralino+gluino contribution
      do m=1,6
         if (in.eq.1) then
            sig_dd_y = sig_dd_y + e2/st2*v2/sq2/2.d0*yd(j)
     $           * gm2*mu*zd(j,m)*dconjg(zd(i,m))
     $           * cp0(abs(mu),abs(gm2),sdm(m))
     $           + e2/ct2*v2/sq2*gm3*mu/6.d0
     $           * cp0(abs(mu),abs(gm3),sdm(m))
     $           * (dconjg(zd(i,m))*zd(j,m)*yd(j)
     $           + 2*dconjg(zd(i+3,m))*zd(j+3,m)*yd(i))
         end if
         do n=1,6
            tmp = (0.d0,0.d0)
            do l=1,3
               tmp = tmp + zd(l,m)*dconjg(zd(l+3,n))*yd(l)
            end do
            if (in.eq.1) then
               sig_dd_y = sig_dd_y + mu*v2/sq2*dconjg(zd(i,m))*zd(j+3,n)
     $              * e2/ct2/9.d0*gm3*cp0(abs(gm3),sdm(m),sdm(n))*tmp
            end if
            if (ig.eq.1) then
               sig_dd_y = sig_dd_y - mu*v2/sq2*dconjg(zd(i,m))*zd(j+3,n)
     $              * 8/3.d0*g3d*g3d*gm1*cp0(gm1,sdm(m),sdm(n))*tmp
            end if
         end do
      end do
c     Chargino contribution
      if (ic.eq.1) then
         do m=1,6
            sig_dd_y = sig_dd_y + ic*e2/st2*v2/sq2*yd(j)*gm2*mu
     $           * dconjg(zd(i,m))*zd(j,m)*cp0(abs(mu),abs(gm2),sdm(m))
            do n=1,6
               tmp  = (0.d0,0.d0)
               tmp1 = (0.d0,0.d0)
               tmp2 = (0.d0,0.d0)
               do l=1,3
                  tmp  = tmp + yu(l)*dconjg(ckm(l,i)*zu(l+3,m))
                  tmp1 = tmp1 + ckm(l,j)*zu(l,n)
                  do k=1,3
                     tmp2 = tmp2 + zu(k+3,m)*(v2*us(l,k) + v1*ws(l,k))
     $                    * dconjg(zu(l,n))
                  end do
                  tmp2 = tmp2 + v1*zu(l+3,m)*dconjg(mu*yu(l)*zu(l,n))
               end do
               sig_dd_y = sig_dd_y + yd(j)/sq2*mu*tmp*tmp1*tmp2*
     $              cp0(abs(mu),sum(m),sum(n))
            end do
         end do
      end if
      sig_dd_y = - sig_dd_y/16/pi/pi
      return
      end

      double complex function sig_dd_ny(i)
c     Chirally enhanced diagonal d-quark self energy, terms *not
c     propotional* to Yukawa_d, in the decoupling limit
      implicit double precision (a-h,o-z)
      double complex gm2,gm3
      double complex zu,zd
      double complex ls,ks,ds,es,us,ws
      double complex tmp
      logical init_alpha_susy
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      common/vev/v1,v2
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/debug_4q/ih,ic,in,ig
      if (init_alpha_susy) call init_alpha_s_susy
c     Neutralino+gluino contribution
      sig_dd_ny = (0.d0,0.d0)
      do m=1,6
         do n=1,6
            tmp = (0.d0,0.d0)
            do l=1,3
               do k=1,3
                  tmp = tmp + zd(l,m)*dconjg(zd(k+3,n)*(v1*ds(l,k)
     $                 - v2*es(l,k)))
               end do
            end do
            if (in.eq.1) then
               sig_dd_ny = sig_dd_ny + tmp/sq2*dconjg(zd(i,m))*zd(i+3,n)
     $              * e2/ct2/9.d0*gm3*cp0(abs(gm3),sdm(m),sdm(n))
            end if
            if (ig.eq.1) then
               sig_dd_ny = sig_dd_ny - tmp/sq2*dconjg(zd(i,m))*zd(i+3,n)
     $              * 8/3.d0*g3d*g3d*gm1*cp0(gm1,sdm(m),sdm(n))
            end if
         end do
      end do
      sig_dd_ny = - sig_dd_ny/16/pi/pi
      return
      end

      double complex function sig_dd_nckm(i,j)
c     Chirally enhanced d-quark self energy, terms *not propotional* to
c     CKM elements, in the decoupling limit
      implicit double precision (a-h,o-z)
      double complex mu
      double complex gm2,gm3
      double complex zu,zd
      double complex ls,ks,ds,es,us,ws
      double complex tmp
      double complex yl,yu,yd
      logical init_alpha_susy
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      common/vev/v1,v2
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/hpar/hm1,hm2,hs,mu
      common/yukawa/yl(3),yu(3),yd(3)
      common/debug_4q/ih,ic,in,ig
      if (init_alpha_susy) call init_alpha_s_susy
c     Neutralino+gluino contribution, terms without Yukawa couplings
      sig_dd_nckm = (0.d0,0.d0)
      do m=1,6
         do n=1,6
            tmp = (0.d0,0.d0)
            do l=1,3
               do k=1,3
                 tmp = tmp + zd(l,m)*dconjg(zd(k+3,n)*(v1*ds(l,k)
     $                 - v2*es(l,k)))
               end do
            end do
            if (in.eq.1) then
               sig_dd_nckm = sig_dd_nckm + tmp*dconjg(zd(i,m))*zd(j+3,n)
     $              * e2/ct2/sq2/9.d0*gm3*cp0(abs(gm3),sdm(m),sdm(n))
            end if
            if (ig.eq.1) then
               sig_dd_nckm = sig_dd_nckm - tmp*dconjg(zd(i,m))*zd(j+3,n)
     $              * 8/sq2/3.d0*g3d*g3d*gm1*cp0(gm1,sdm(m),sdm(n))
            end if
         end do
      end do
c     Neutralino+gluino contribution, terms proportional to Yukawa couplings
      do m=1,6
         if (in.eq.1) then
            sig_dd_nckm = sig_dd_nckm + e2/st2*v2/sq2/2.d0*yd(j)
     $           * gm2*mu*zd(j,m)*dconjg(zd(i,m))
     $           * cp0(abs(mu),abs(gm2),sdm(m))
     $           + e2/ct2*v2/sq2*gm3*mu/6.d0
     $           * cp0(abs(mu),abs(gm3),sdm(m))
     $           * (dconjg(zd(i,m))*zd(j,m)*yd(j)
     $           + 2*dconjg(zd(i+3,m))*zd(j+3,m)*yd(i))
         end if
         do n=1,6
            tmp = (0.d0,0.d0)
            do l=1,3
               tmp = tmp + zd(l,m)*dconjg(zd(l+3,n))*yd(l)
            end do
            if (in.eq.1) then
               sig_dd_nckm = sig_dd_nckm
     $              + tmp*v2/sq2*dconjg(zd(i,m))*zd(j+3,n)*mu
     $              * e2/ct2/9.d0*gm3*cp0(abs(gm3),sdm(m),sdm(n))
            end if
            if (ig.eq.1) then
               sig_dd_nckm = sig_dd_nckm
     $              - tmp*v2/sq2*dconjg(zd(i,m))*zd(j+3,n)*mu
     $              * 8/3.d0*g3d*g3d*gm1*cp0(gm1,sdm(m),sdm(n))
            end if
         end do
      end do
c     Chargino contribution
      if (ic.eq.1) then
         do m=1,6
            sig_dd_nckm = sig_dd_nckm + e2/st2*v2/sq2*yd(j)*gm2*mu
     $           * dconjg(zd(i,m))*zd(j,m)*cp0(abs(mu),abs(gm2),sdm(m))
         end do
      end if
      sig_dd_nckm = - sig_dd_nckm/16/pi/pi
      return
      end

      double complex function sig_dd_ckm(i,j)
c     Chirally enhanced d-quark self energy, terms propotional to CKM
c     elements, in the decoupling limit
      implicit double precision (a-h,o-z)
      double complex yl,yu,yd
      double complex mu,ckm,zu,zd
      double complex ls,ks,ds,es,us,ws
      double complex tmp,tmp1,tmp2
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/vev/v1,v2
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/hpar/hm1,hm2,hs,mu
      common/yukawa/yl(3),yu(3),yd(3)
      common/ckm/ckm(3,3)
      common/debug_4q/ih,ic,in,ig
      sig_dd_ckm = (0.d0,0.d0)
c     Chargino contribution only
      if (ic.ne.1) return
      do m=1,6
         do n=1,6
            tmp  = (0.d0,0.d0)
            tmp1 = (0.d0,0.d0)
            tmp2 = (0.d0,0.d0)
            do l=1,3
               tmp  = tmp + yu(l)*dconjg(ckm(l,i)*zu(l+3,m))
               tmp1 = tmp1 + ckm(l,j)*zu(l,n)
               do k=1,3
                  tmp2 = tmp2 + zu(k+3,m)*(v2*us(l,k) + v1*ws(l,k))
     $                 * dconjg(zu(l,n))
               end do
               tmp2 = tmp2 + v1*zu(l+3,m)*dconjg(mu*yu(l)*zu(l,n))
            end do
            sig_dd_ckm = sig_dd_ckm + ic/sq2*yd(j)*mu*tmp*tmp1*tmp2
     $           * cp0(abs(mu),sum(m),sum(n))
         end do
      end do
      sig_dd_ckm = - sig_dd_ckm/16/pi/pi
      return
      end

      double complex function sig_dd_nhol(i,j)
c     Chirally enhanced d-quark self energy, "non-holomorphic" terms, i.e.
c     propotional to v2 in the decoupling limit
      implicit double precision (a-h,o-z)
      double complex tmp,tmp1,tmp2
      double complex mu
      double complex gm2,gm3
      double complex zu,zd
      double complex ls,ks,ds,es,us,ws
      double complex yl,yu,yd
      double complex ckm
      logical init_alpha_susy
      common/ckm/ckm(3,3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      common/vev/v1,v2
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/hpar/hm1,hm2,hs,mu
      common/yukawa/yl(3),yu(3),yd(3)
      common/debug_4q/ih,ic,in,ig
      if (init_alpha_susy) call init_alpha_s_susy
c     Neutralino+gluino contribution, terms without Yukawa couplings
      sig_dd_nhol = (0.d0,0.d0)
      do m=1,6
         do n=1,6
            tmp = (0.d0,0.d0)
            do l=1,3
               do k=1,3
                  tmp = tmp - zd(l,m)*dconjg(zd(k+3,n)*es(l,k))
               end do
            end do
            if (in.eq.1) then
               sig_dd_nhol = sig_dd_nhol + tmp*dconjg(zd(i,m))*zd(j+3,n)
     $              * e2/ct2/9.d0*gm3*cp0(abs(gm3),sdm(m),sdm(n))
            end if
            if (ig.eq.1) then
               sig_dd_nhol = sig_dd_nhol - tmp*dconjg(zd(i,m))*zd(j+3,n)
     $              * 8/3.d0*g3d*g3d*gm1*cp0(gm1,sdm(m),sdm(n))
            end if
         end do
      end do
c     Neutralino+gluino contribution, terms proportional to Yukawa couplings
      do m=1,6
         if (in.eq.1) then
            sig_dd_nhol = sig_dd_nhol + e2/st2/2.d0*yd(j)*gm2*mu
     $           * zd(j,m)*dconjg(zd(i,m))*cp0(abs(mu),abs(gm2),sdm(m))
     $           + e2/ct2*gm3*mu/6.d0*cp0(abs(mu),abs(gm3),sdm(m))
     $           * (dconjg(zd(i,m))*zd(j,m)*yd(j)
     $           + 2*dconjg(zd(i+3,m))*zd(j+3,m)*yd(i))
         end if
         do n=1,6
            tmp = (0.d0,0.d0)
            do l=1,3
               tmp = tmp + zd(l,m)*dconjg(zd(l+3,n))*yd(l)
            end do
            if (in.eq.1) then
               sig_dd_nhol = sig_dd_nhol + tmp*dconjg(zd(i,m))*zd(j+3,n)
     $              * mu*e2/ct2/9.d0*gm3*cp0(abs(gm3),sdm(m),sdm(n))
            end if
            if (ig.eq.1) then
               sig_dd_nhol = sig_dd_nhol - tmp*dconjg(zd(i,m))*zd(j+3,n)
     $              * mu*8/3.d0*g3d*g3d*gm1*cp0(gm1,sdm(m),sdm(n))
            end if
         end do
      end do
c     Chargino contribution
      if (ic.eq.1) then
         do m=1,6
            sig_dd_nhol = sig_dd_nhol + e2/st2*yd(j)*gm2*mu*zd(j,m)
     $           * dconjg(zd(i,m))*cp0(abs(mu),abs(gm2),sdm(m))
            do n=1,6
               tmp  = (0.d0,0.d0)
               tmp1 = (0.d0,0.d0)
               tmp2 = (0.d0,0.d0)
               do l=1,3
                  tmp  = tmp + yu(l)*dconjg(ckm(l,i)*zu(l+3,m))
                  tmp1 = tmp1 + ckm(l,j)*zu(l,n)
                  do k=1,3
                     tmp2 = tmp2 + zu(k+3,m)*us(l,k)*dconjg(zu(l,n))
                  end do
               end do
               sig_dd_nhol = sig_dd_nhol + yd(j)*mu*tmp*tmp1*tmp2
     $              * cp0(abs(mu),sum(m),sum(n))
            end do
         end do
      end if
      sig_dd_nhol = - v2/sq2*sig_dd_nhol/16/pi/pi
      return
      end

      double complex function sig_dd_hol(i,j)
c     Chirally enhanced d-quark self energy, "holomorphic" terms, i.e.
c     propotional to v1 in the decoupling limit
      implicit double precision (a-h,o-z)
      double complex gm2,gm3
      double complex zu,zd
      double complex ls,ks,ds,es,us,ws
      double complex tmp,tmp1,tmp2
      double complex ckm
      double complex mu
      double complex yl,yu,yd
      logical init_alpha_susy
      common/hpar/hm1,hm2,hs,mu
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      common/vev/v1,v2
      common/ckm/ckm(3,3)
      common/yukawa/yl(3),yu(3),yd(3)
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/debug_4q/ih,ic,in,ig
      if (init_alpha_susy) call init_alpha_s_susy
      sig_dd_hol = (0.d0,0.d0)
c     Neutralino+gluino contribution
      do m=1,6
         do n=1,6
            tmp = (0.d0,0.d0)
            do l=1,3
               do k=1,3
                 tmp = tmp + zd(l,m)*dconjg(ds(l,k)*zd(k+3,n))
               end do
            end do
            if (in.eq.1) then
               sig_dd_hol = sig_dd_hol + tmp*dconjg(zd(i,m))*zd(j+3,n)
     $              * e2/ct2/9.d0*gm3*cp0(abs(gm3),sdm(m),sdm(n))
            end if
            if (ig.eq.1) then
               sig_dd_hol = sig_dd_hol - tmp*dconjg(zd(i,m))*zd(j+3,n)
     $              * 8/3.d0*g3d*g3d*gm1*cp0(gm1,sdm(m),sdm(n))
            end if
         end do
      end do
c     Chargino contribution
      if (ic.eq.1) then
         do m=1,6
            do n=1,6
               tmp  = (0.d0,0.d0)
               tmp1 = (0.d0,0.d0)
               tmp2 = (0.d0,0.d0)
               do l=1,3
                  tmp  = tmp + yu(l)*dconjg(ckm(l,i)*zu(l+3,m))
                  tmp1 = tmp1 + ckm(l,j)*zu(l,n)
                  do k=1,3
                     tmp2 = tmp2 + zu(k+3,m)*ws(l,k)*dconjg(zu(l,n))
                  end do
                  tmp2 = tmp2 + zu(l+3,m)*dconjg(mu*yu(l)*zu(l,n))
               end do
               sig_dd_hol = sig_dd_hol + yd(j)*mu*tmp*tmp1*tmp2
     $              * cp0(abs(mu),sum(m),sum(n))
            end do
         end do
      end if
      sig_dd_hol = - v1/sq2*sig_dd_hol/16/pi/pi
      return
      end

      double complex function sig_uu_ny(i,j)
c     Chirally not-suppressed diagonal u-quark self energy, terms *not
c     propotional* to Yukawa_u, in the decoupling limit
      implicit double precision (a-h,o-z)
      double complex gm2,gm3
      double complex zu,zd
      double complex ls,ks,ds,es,us,ws
      double complex tmp
      logical init_alpha_susy
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      common/vev/v1,v2
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/debug_4q/ih,ic,in,ig
      if (init_alpha_susy) call init_alpha_s_susy
      sig_uu_ny = (0.d0,0.d0)
c     Neutralino+gluino contribution
      do m=1,6
         do n=1,6
            tmp = (0.d0,0.d0)
            do l=1,3
               do k=1,3
                  tmp = tmp + dconjg(zu(l,m)*(v2*us(l,k) + v1*ws(l,k)))
     $                 * zu(k+3,n)
               end do
            end do
            if (in.eq.1) then
               sig_uu_ny = sig_uu_ny + tmp*zu(i,m)*dconjg(zu(j+3,n))
     $              * e2/ct2/9*dconjg(gm3)*cp0(abs(gm3),sum(m),sum(n))
            end if
            if (ig.eq.1) then
               sig_uu_ny = sig_uu_ny + tmp*zu(i,m)*dconjg(zu(j+3,n))
     $              * 4/3.d0*g3u*g3u*gm1*cp0(gm1,sum(m),sum(n))
            end if
         end do
      end do
      sig_uu_ny = - sq2*sig_uu_ny/16/pi/pi
      return
      end

      double complex function sig_uu_hol(i,j)
c     Chirally enhanced u-quark self energy, "holomorphic" terms, i.e.
c     propotional to v2 in the decoupling limit
      implicit double precision (a-h,o-z)
      double complex gm2,gm3
      double complex zu,zd
      double complex ls,ks,ds,es,us,ws
      double complex tmp,tmp1,tmp2
      double complex ckm
      double complex yl,yu,yd
      double complex mu
      logical init_alpha_susy
      common/hpar/hm1,hm2,hs,mu
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      common/yukawa/yl(3),yu(3),yd(3)
      common/ckm/ckm(3,3)
      common/vev/v1,v2
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/debug_4q/ih,ic,in,ig
      if (init_alpha_susy) call init_alpha_s_susy
c     Neutralino+gluino contribution only
      sig_uu_hol = (0.d0,0.d0)
      do m=1,6
         do n=1,6
            tmp = (0.d0,0.d0)
            do k=1,3
               do l=1,3
                  tmp = tmp + dconjg(us(l,k)*zu(l,m))*zu(k+3,n)
               end do
            end do
            if (in.eq.1) then
               sig_uu_hol = sig_uu_hol + tmp*zu(i,m)*dconjg(zu(j+3,n))
     $              * 2/9.d0*e2/ct2*gm3*cp0(abs(gm3),sum(m),sum(n))
            end if
            if (ig.eq.1) then
               sig_uu_hol = sig_uu_hol + tmp*zu(i,m)*dconjg(zu(j+3,n))
     $              * 8/3.d0*g3u*g3u*gm1*cp0(gm1,sum(m),sum(n))
            end if
         end do
      end do
      if (ic.eq.1) then
c     Chargino contribution
         do m=1,6
            do n=1,6
               tmp  = (0.d0,0.d0)
               tmp1 = (0.d0,0.d0)
               tmp2 = (0.d0,0.d0)
               do k=1,3
                  tmp1 = tmp1 + zd(k,m)*dconjg(ckm(j,k))
                  tmp2 = tmp2 + ckm(i,k)*dconjg(zd(k+3,n))*yd(k)
                  tmp  = tmp  + dconjg(mu*yd(k)*zd(k,m))*zd(k+3,n)
                  do l=1,3
                     tmp = tmp - zd(k+3,n)*es(l,k)*dconjg(zd(l,m))
                  end do
               end do
               sig_uu_hol = sig_uu_hol - mu*yu(j)*tmp1*tmp2*tmp
     $              * cp0(abs(mu),sdm(m),sdm(n))
            end do
         end do
      end if
      sig_uu_hol = - v2/sq2*sig_uu_hol/16/pi/pi
      return
      end

      double complex function sig_uu_nhol(i,j)
c     Chirally enhanced u-quark self energy, "nholomorphic" terms, i.e.
c     propotional to v2 in the decoupling limit
      implicit double precision (a-h,o-z)
      double complex gm2,gm3
      double complex zu,zd
      double complex ls,ks,ds,es,us,ws
      double complex tmp,tmp1,tmp2
      double complex ckm
      logical init_alpha_susy
      double complex mu
      double complex yl,yu,yd
      common/hpar/hm1,hm2,hs,mu
      common/yukawa/yl(3),yu(3),yd(3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      common/ckm/ckm(3,3)
      common/vev/v1,v2
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/debug_4q/ih,ic,in,ig
      if (init_alpha_susy) call init_alpha_s_susy
c     Neutralino+gluino contribution
      sig_uu_nhol = (0.d0,0.d0)
      do m=1,6
         if (in.eq.1) then
            sig_uu_nhol = sig_uu_nhol - 2/3.d0*e2/ct2*gm3*mu*yu(i)
     $           * zu(i+3,m)*dconjg(zu(j+3,m))
     $           * cp0(abs(mu),abs(gm3),sum(m))
     $           + e2/2.d0*zu(i,m)*dconjg(zu(j,m))*mu*yu(j)
     $           * (gm3/ct2/3*cp0(abs(mu),abs(gm3),sum(m))
     $           - gm2/st2*cp0(abs(mu),abs(gm2),sum(m)))
         end if
         do n=1,6
            tmp = (0.d0,0.d0)
            do k=1,3
               tmp = tmp + mu*yu(k)*dconjg(zu(k,m))*zu(k+3,n)
               do l=1,3
                  tmp = tmp + dconjg(ws(l,k)*zu(l,m))*zu(k+3,n)
               end do
            end do
         if (in.eq.1) then
            sig_uu_nhol = sig_uu_nhol + tmp*zu(i,m)*dconjg(zu(j+3,n))
     $           * 2/9.d0*e2/ct2*gm3*cp0(abs(gm3),sum(m),sum(n))
         end if
         if (ig.eq.1) then
            sig_uu_nhol = sig_uu_nhol + tmp*zu(i,m)*dconjg(zu(j+3,n))
     $           * 8/3.d0*g3u*g3u*gm1*cp0(gm1,sum(m),sum(n))
         end if
      end do
      end do
      if (ic.eq.1) then
c     Chargino contribution
         do m=1,6
            sig_uu_nhol = sig_uu_nhol - ic*e2/st2*mu*gm2*zu(i,m)
     $           * yu(j)*dconjg(zu(j,m))*cp0(abs(mu),abs(gm2),sum(m))
            do n=1,6
               tmp  = (0.d0,0.d0)
               tmp1 = (0.d0,0.d0)
               tmp2 = (0.d0,0.d0)
               do k=1,3
                  tmp1 = tmp1 + zd(k,n)*dconjg(ckm(j,k))
                  tmp2 = tmp2 + ckm(i,k)*dconjg(zd(k+3,m))*yd(k)
                  do l=1,3
                     tmp = tmp + zd(k+3,m)*ds(l,k)*dconjg(zd(l,n))
                  end do
               end do
               sig_uu_nhol = sig_uu_nhol - mu*yu(j)*tmp1*tmp2*tmp
     $              * cp0(abs(mu),sdm(m),sdm(n))
            end do
         end do
      end if
      sig_uu_nhol = - v1/sq2*sig_uu_nhol/16/pi/pi
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Complete up and down quark self energies, calculated at p^2=0    c  
c     without the resummation of enhanced chiral corrections!          c
c     Non-decoupling inifinite terms truncated                         c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Left d-quark self-energy (proportional to P_L)                    c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dsl_sig0_h(i,j)
c     Up quark + charged Higgs in loop
      implicit double precision (a-h,o-z)
      double complex b0,ckm
      double complex yl,yu,yd
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/km_mat/ckm(3,3)
      common/yukawa/yl(3),yu(3),yd(3)
      common/fmass/em(3),um(3),dm(3)
      dsl_sig0_h = (0.d0,0.d0)
      do l=1,3
         do k=1,2
            dsl_sig0_h = dsl_sig0_h + zh(2,k)*zh(1,k)*um(l)*yu(l)*yd(j)
     $           * dconjg(ckm(i,l))*ckm(j,l)*b0(0.d0,um(l),cm(k))
         end do
      end do
      return
      end
      
      double complex function dsl_sig0_n(i,j)
c     Neutralino and down squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd,zn
      double complex vl_ddn0,vr_ddn0
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      dsl_sig0_n = (0.d0,0.d0)
      do k=1,6
         do l=1,4
            dsl_sig0_n = dsl_sig0_n - vl_ddn0(i,k,l)
     $           * dconjg(vr_ddn0(j,k,l))
     $           * fnm(l)*b0(0.d0,fnm(l),sdm(k))
         end do
      end do
      return
      end

      double complex function dsl_sig0_c(i,j)
c     Chargino and up squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd,zpos,zneg
      double complex vl_duc0,vr_duc0
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      dsl_sig0_c = (0.d0,0.d0)
      do k=1,6
         do l=1,2
            dsl_sig0_c = dsl_sig0_c - vl_duc0(i,k,l)
     $           * dconjg(vr_duc0(j,k,l))
     $           * fcm(l)*b0(0.d0,fcm(l),sum(k))
        end do
      end do
      return
      end

      double complex function dsl_sig0_g(i,j)
c     Gluino and down squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd
      double complex gm2,gm3
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      dsl_sig0_g = (0.d0,0.d0)
      al = 8*g3d*g3d/3.d0
      do  k=1,6
         dsl_sig0_g = dsl_sig0_g + al*zd(i,k)*dconjg(zd(j+3,k))*gm1 
     $        * b0(0.d0,gm1,sdm(k))
      end do
      return
      end
      
      double complex function dsl_sig(i,j)
      implicit double precision (a-h,o-z)
c      double complex dsl_sig0_h
      double complex dsl_sig0_n,dsl_sig0_c,dsl_sig0_g
      logical cdiag,ndiag,tmp
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      common/nc_suppress/eps_d,eps_u,acc
      external init_d_self,init_nc_diag
      dsl_sig = (0.d0,0.d0)
c      if (ih.eq.1) dsl_sig = dsl_sig + dsl_sig0_h(i,j) 
      if ((in.eq.1).or.(ic.eq.1)) then
         eps = eps_d
         eps_d = acc
         tmp = cdiag().or.ndiag()
      end if
      if (in.eq.1) dsl_sig = dsl_sig + dsl_sig0_n(i,j)
      if (ic.eq.1) dsl_sig = dsl_sig + dsl_sig0_c(i,j)
      if ((in.eq.1).or.(ic.eq.1)) then
         eps_d = eps
         tmp = cdiag().or.ndiag()
      end if
      if (ig.eq.1) dsl_sig = dsl_sig + dsl_sig0_g(i,j) 
      dsl_sig = dsl_sig/16/pi/pi
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Right d-quark self-energy (proportional to P_R)                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dsr_sig0_h(i,j)
c     Up quark + charged Higgs in loop
      implicit double precision (a-h,o-z)
      double complex b0,ckm
      double complex yl,yu,yd
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/km_mat/ckm(3,3)
      common/yukawa/yl(3),yu(3),yd(3)
      common/fmass/em(3),um(3),dm(3)
      dsr_sig0_h = (0.d0,0.d0)
      do l=1,3
         do k=1,2
            dsr_sig0_h = dsr_sig0_h + zh(2,k)*zh(1,k)*um(l)*yu(l)*yd(i)
     $           * dconjg(ckm(i,l))*ckm(j,l)*b0(0.d0,um(l),cm(k))
         end do
      end do
      return
      end

      double complex function dsr_sig0_n(i,j)
c     Neutralino and down squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd,zn
      double complex vl_ddn0,vr_ddn0
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      dsr_sig0_n = (0.d0,0.d0)
      do k=1,6
         do l=1,4
            dsr_sig0_n = dsr_sig0_n - vr_ddn0(i,k,l)
     $           * dconjg(vl_ddn0(j,k,l))
     $           * fnm(l)*b0(0.d0,fnm(l),sdm(k))
        end do
      end do
      return
      end

      double complex function dsr_sig0_c(i,j)
c     Chargino and up squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd,zpos,zneg
      double complex vl_duc0,vr_duc0
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      dsr_sig0_c = (0.d0,0.d0)
      do k=1,6
         do l=1,2
            dsr_sig0_c = dsr_sig0_c - vr_duc0(i,k,l)
     $           * dconjg(vl_duc0(j,k,l))
     $           * fcm(l)*b0(0.d0,fcm(l),sum(k))
        end do
      end do
      return
      end

      double complex function dsr_sig0_g(i,j)
c     Gluino and down squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd
      double complex gm2,gm3
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      dsr_sig0_g = (0.d0,0.d0)
      al = 8*g3d*g3d/3.d0
      do  k=1,6
         dsr_sig0_g = dsr_sig0_g + al*zd(i+3,k)*dconjg(zd(j,k))
     $        * gm1*b0(0.d0,gm1,sdm(k))
      end do
      return
      end

      double complex function dsr_sig(i,j)
      implicit double precision (a-h,o-z)
c      double complex dsr_sig0_h
      double complex dsr_sig0_n,dsr_sig0_c,dsr_sig0_g
      logical cdiag,ndiag,tmp
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      common/nc_suppress/eps_d,eps_u,acc
      external init_d_self,init_nc_diag
      dsr_sig = (0.d0,0.d0)
c      if (ih.eq.1) dsr_sig = dsr_sig + dsr_sig0_h(i,j) 
      if ((in.eq.1).or.(ic.eq.1)) then
         eps = eps_d
         eps_d = acc
         tmp = cdiag().or.ndiag()
      end if
      if (in.eq.1) dsr_sig = dsr_sig + dsr_sig0_n(i,j)
      if (ic.eq.1) dsr_sig = dsr_sig + dsr_sig0_c(i,j)
      if ((in.eq.1).or.(ic.eq.1)) then
         eps_d = eps
         tmp = cdiag().or.ndiag()
      end if
      if (ig.eq.1) dsr_sig = dsr_sig + dsr_sig0_g(i,j) 
      dsr_sig = dsr_sig/16/pi/pi
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Left u-quark self-energy (proportional to P_L)                    c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function usl_sig0_h(i,j)
c     Down quark + charged Higgs in loop
      implicit double precision (a-h,o-z)
      double complex b0,ckm
      double complex yl,yu,yd
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/km_mat/ckm(3,3)
      common/yukawa/yl(3),yu(3),yd(3)
      usl_sig0_h = (0.d0,0.d0)
      do l=1,3
         do k=1,2
            usl_sig0_h = usl_sig0_h + zh(2,k)*zh(1,k)*dm(l)*yd(l)*yu(i)
     $           * ckm(l,i)*dconjg(ckm(l,j))*b0(0.d0,dm(l),cm(k))
         end do
      end do
      return
      end

      double complex function usl_sig0_n(i,j)
c     Neutralino and down squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd,zn
      double complex vl_uun0,vr_uun0
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      usl_sig0_n = (0.d0,0.d0)
      do k=1,6
         do l=1,4
            usl_sig0_n = usl_sig0_n -  vr_uun0(i,k,l)
     $           * dconjg(vl_uun0(j,k,l))
     $           * fnm(l)*b0(0.d0,fnm(l),sum(k))
         end do
      end do
      return
      end

      double complex function usl_sig0_c(i,j)
c     Chargino and down squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd,zpos,zneg
      double complex vl_udc0,vr_udc0
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      usl_sig0_c = (0.d0,0.d0)
      do k=1,6
         do l=1,2
            usl_sig0_c = usl_sig0_c - vr_udc0(i,k,l)
     $           * dconjg(vl_udc0(j,k,l))
     $           * fcm(l)*b0(0.d0,fcm(l),sdm(k))
         end do
      end do
      return
      end

      double complex function usl_sig0_g(i,j)
c     gluino and up squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd
      double complex gm2,gm3
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      al = 8/3.d0*g3u*g3u
      usl_sig0_g = (0.d0,0.d0)
      do k=1,6
         usl_sig0_g = usl_sig0_g + al*zu(j,k)*dconjg(zu(i+3,k))*gm1
     $        * b0(0.d0,gm1,sum(k))
      end do
      return
      end

      double complex function usl_sig(i,j)
      implicit double precision (a-h,o-z)
c      double complex usl_sig0_h
      double complex usl_sig0_n,usl_sig0_c,usl_sig0_g
      logical cdiag,ndiag,tmp
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      common/nc_suppress/eps_d,eps_u,acc
      external init_d_self,init_nc_diag
      usl_sig = (0.d0,0.d0)
c      if (ih.eq.1) usl_sig = usl_sig + usl_sig0_h(i,j) 
      if ((in.eq.1).or.(ic.eq.1)) then
         eps = eps_u
         eps_u = acc
         tmp = cdiag().or.ndiag()
      end if
      if (in.eq.1) usl_sig = usl_sig + usl_sig0_n(i,j)
      if (ic.eq.1) usl_sig = usl_sig + usl_sig0_c(i,j)
      if ((in.eq.1).or.(ic.eq.1)) then
         eps_u = eps
         tmp = cdiag().or.ndiag()
      end if
      if (ig.eq.1) usl_sig = usl_sig + usl_sig0_g(i,j) 
      usl_sig = usl_sig/16/pi/pi
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Right u-quark self-energy (proportional to P_R)                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function usr_sig0_h(i,j)
c     Down quark + charged Higgs in loop
      implicit double precision (a-h,o-z)
      double complex b0,ckm
      double complex yl,yu,yd
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/km_mat/ckm(3,3)
      common/yukawa/yl(3),yu(3),yd(3)
      usr_sig0_h = (0.d0,0.d0)
      do l=1,3
         do k=1,2
            usr_sig0_h = usr_sig0_h + zh(2,k)*zh(1,k)*dm(l)*yd(l)*yu(j)
     $           * ckm(l,i)*dconjg(ckm(l,j))*b0(0.d0,dm(l),cm(k))
         end do
      end do
      return
      end

      double complex function usr_sig0_n(i,j)
c     Neutralino and down squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd,zn
      double complex vl_uun0,vr_uun0
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      usr_sig0_n = (0.d0,0.d0)
      do k=1,6
         do l=1,4
            usr_sig0_n = usr_sig0_n - vl_uun0(i,k,l)
     $           * dconjg(vr_uun0(j,k,l))
     $           * fnm(l)*b0(0.d0,fnm(l),sum(k))
         end do
      end do
      return
      end

      double complex function usr_sig0_c(i,j)
c     Chargino and down squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd,zpos,zneg
      double complex vl_udc0,vr_udc0
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      usr_sig0_c = (0.d0,0.d0)
      do k=1,6
         do l=1,2
            usr_sig0_c = usr_sig0_c - vl_udc0(i,k,l)
     $           * dconjg(vr_udc0(j,k,l))
     $           * fcm(l)*b0(0.d0,fcm(l),sdm(k))
         end do
      end do
      return
      end

      double complex function usr_sig0_g(i,j)
c     gluino and up squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd
      double complex gm2,gm3
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      al = 8/3.d0*g3u*g3u
      usr_sig0_g = (0.d0,0.d0)
      do k=1,6
         usr_sig0_g = usr_sig0_g + al*zu(j+3,k)*dconjg(zu(i,k))*gm1
     $        * b0(0.d0,gm1,sum(k))
      end do
      return
      end

      double complex function usr_sig(i,j)
      implicit double precision (a-h,o-z)
c     double complex usr_sig0_h
      double complex usr_sig0_n,usr_sig0_c,usr_sig0_g
      logical cdiag,ndiag,tmp
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      common/nc_suppress/eps_d,eps_u,acc
      external init_d_self,init_nc_diag
      usr_sig = (0.d0,0.d0)
c     if (ih.eq.1) usr_sig = usr_sig + usr_sig0_h(i,j) 
      if ((in.eq.1).or.(ic.eq.1)) then
         eps = eps_u
         eps_u = acc
         tmp = cdiag().or.ndiag()
      end if
      if (in.eq.1) usr_sig = usr_sig + usr_sig0_n(i,j)
      if (ic.eq.1) usr_sig = usr_sig + usr_sig0_c(i,j)
      if ((in.eq.1).or.(ic.eq.1)) then
         eps_u = eps
         tmp = cdiag().or.ndiag()
      end if
      if (ig.eq.1) usr_sig = usr_sig + usr_sig0_g(i,j) 
      usr_sig = usr_sig/16/pi/pi
      return
      end

