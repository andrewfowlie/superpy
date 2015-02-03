c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM}
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: L_SELF0_DLIM.F
c     Released: 20: 2:2012(J.R.)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains initialization routines for the various       c
c     decompositions of the lepton self energies, calculated           c
c     at p^2 = 0 and in the decoupling limit only                      c
c     Conventions follow the arXiv:1103.4272 [hep-ph], i.e. comparing  c
c     to l_self0.f one has sig_ll(i,j) ~ esl_sig(i,j)^*                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function eps_ll(i)
c     Diagonal "epsilon_d" terms from d-quark self energies in the
c     decoupling limit
      implicit double precision (a-h,o-z)
      double complex gm2,gm3,zv,zl,mu
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/gmass/gm1,gm2,gm3
      common/hpar/hm1,hm2,hs,mu
      common/debug_4q/ih,ic,in,ig
      eps_ll = (0.d0,0.d0)
c     Neutralino contribution, terms proportional to Yukawa 
      if (in.eq.1) then
         do m=1,6
            eps_ll = eps_ll - e2/st2/2.d0*gm2*mu*abs(zl(i,m))**2
     $           * cp0(abs(mu),abs(gm2),slm(m))
     $           + e2/ct2*gm3*mu/2.d0*(abs(zl(i,m))**2
     $           - 2*abs(zl(i+3,m))**2)*cp0(abs(mu),abs(gm3),slm(m))
            do n=1,6
               eps_ll = eps_ll + e2/ct2*abs(zl(i,m)*zl(i+3,n))**2*mu*gm3
     $              * cp0(abs(gm3),slm(m),slm(n))
            end do
         end do
      end if
c     Chargino contribution (sign and yukawa)
      if (ic.eq.1) then
         do m=1,3
            eps_ll = eps_ll - e2/st2*mu*gm2*abs(zv(i,m))**2
     $           * cp0(abs(mu),abs(gm2),vm(m))
         end do
      end if
      eps_ll = - eps_ll/16/pi/pi
      return
      end

      double complex function sig_ll_yt(i)
c     e and mu self energy, part propotional to Ytau, in the decoupling limit
      implicit double precision (a-h,o-z)
      double complex mu
      double complex gm2,gm3
      double complex zv,zl
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/gmass/gm1,gm2,gm3
      common/vev/v1,v2
      common/hpar/hm1,hm2,hs,mu
      common/yukawa/yl(3),yu(3),yd(3)
      common/debug_4q/ih,ic,in,ig
      sig_ll_yt = (0.d0,0.d0)
      if (in.ne.1) return
c     Neutralino contribution
      do m=1,6
         do n=1,6
            sig_ll_yt = sig_ll_yt - e2/ct2*yl(3)*v2/sq2*mu*gm3
     $           * dconjg(zl(6,n)*zl(i,m))*zl(3,m)*zl(i+3,n)
     $           * cp0(abs(gm3),slm(m),slm(n))
         end do
      end do
      sig_ll_yt = - sig_ll_yt/16/pi/pi
      return     
      end

      double complex function sig_ll_y(i,j)
c     Chirally enhanced diagonal lepton self energy, terms propotional
c     to Yukawa, in the decoupling limit
      implicit double precision (a-h,o-z)
      double complex tmp
      double complex mu
      double complex gm2,gm3
      double complex zv,zl
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/gmass/gm1,gm2,gm3
      common/vev/v1,v2
      common/hpar/hm1,hm2,hs,mu
      common/yukawa/yl(3),yu(3),yd(3)
      common/debug_4q/ih,ic,in,ig
      sig_ll_y = (0.d0,0.d0)
c     Neutralino contribution
      if (in.eq.1) then
         do m=1,6
            sig_ll_y = sig_ll_y + e2/st2*v2/sq2/2.d0*yl(j)
     $           * gm2*mu*zl(j,m)*dconjg(zl(i,m))
     $           * cp0(abs(mu),abs(gm2),slm(m))
     $           - e2/ct2*v2/sq2/2.d0*gm3*mu
     $           * cp0(abs(mu),abs(gm3),slm(m))
     $           * (dconjg(zl(i,m))*zl(j,m)*yl(j)
     $           - 2*dconjg(zl(i+3,m))*zl(j+3,m)*yl(i))
            do n=1,6
               tmp = (0.d0,0.d0)
               do l=1,3
                  tmp = tmp + zl(l,m)*dconjg(zl(l+3,n))*yl(l)
               end do
               sig_ll_y = sig_ll_y - mu*v2/sq2*dconjg(zl(i,m))*zl(j+3,n)
     $              *e2/ct2*gm3*cp0(abs(gm3),slm(m),slm(n))*tmp
            end do
         end do
      end if
c     Chargino contribution
      if (ic.eq.1) then
         do m=1,3
            sig_ll_y = sig_ll_y + e2/st2*v2/sq2*yl(j)*gm2*mu*zv(i,m)
     $           * dconjg(zv(j,m))*cp0(abs(mu),abs(gm2),vm(m))
         end do
      end if
      sig_ll_y = - sig_ll_y/16/pi/pi
      return
      end

      double complex function sig_ll_ny(i,j)
c     Chirally enhanced diagonal leptons self energy, terms *not
c     propotional* to Yukawa, in the decoupling limit
      implicit double precision (a-h,o-z)
      double complex gm2,gm3
      double complex zv,zl
      double complex ls,ks,ds,es,us,ws
      double complex tmp
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/gmass/gm1,gm2,gm3
      common/vev/v1,v2
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/debug_4q/ih,ic,in,ig
c     Neutralino contribution
      sig_ll_ny = (0.d0,0.d0)
      if (in.ne.1) return
      do m=1,6
         do n=1,6
            tmp = (0.d0,0.d0)
            do l=1,3
               do k=1,3
                  tmp = tmp + zl(l,m)*dconjg(zl(k+3,n)*(v1*ls(l,k)
     $                 - v2*ks(l,k)))
               end do
            end do
            sig_ll_ny = sig_ll_ny - tmp/sq2*dconjg(zl(i,m))*zl(j+3,n)
     $           *e2/ct2*gm3*cp0(abs(gm3),slm(m),slm(n))
         end do
      end do
      sig_ll_ny = - sig_ll_ny/16/pi/pi
      return 
      end

      double complex function sig_ll_nhol(i,j)
c     Chirally enhanced lepton self energy, "non-holomorphic" terms,
c     i.e.  propotional to v2 in the decoupling limit
      implicit double precision (a-h,o-z)
      double complex tmp
      double complex mu
      double complex gm2,gm3
      double complex zv,zl
      double complex ls,ks,ds,es,us,ws
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/gmass/gm1,gm2,gm3
      common/vev/v1,v2
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/hpar/hm1,hm2,hs,mu
      common/yukawa/yl(3),yu(3),yd(3)
      common/debug_4q/ih,ic,in,ig
      sig_ll_nhol = (0.d0,0.d0)
      if (in.eq.1) then
c     Neutralino contribution, terms without Yukawa couplings
         do m=1,6
            do n=1,6
               tmp = (0.d0,0.d0)
               do l=1,3
                  do k=1,3
                     tmp = tmp - zl(l,m)*dconjg(zl(k+3,n)*ks(l,k))
                  end do
               end do
               sig_ll_nhol = sig_ll_nhol - tmp*dconjg(zl(i,m))*zl(j+3,n)
     $              *e2/ct2*gm3*cp0(abs(gm3),slm(m),slm(n))
            end do
         end do
c     Neutralino contribution, terms proportional to Yukawa couplings
         do m=1,6
            sig_ll_nhol = sig_ll_nhol + e2/st2/2.d0*yl(j)*gm2*mu
     $           * zl(j,m)*dconjg(zl(i,m))*cp0(abs(mu),abs(gm2),slm(m))
     $           - e2/ct2*gm3*mu/2.d0*cp0(abs(mu),abs(gm3),slm(m))
     $           * (dconjg(zl(i,m))*zl(j,m)*yl(j)
     $           - 2*dconjg(zl(i+3,m))*zl(j+3,m)*yl(i))
            do n=1,6
               tmp = (0.d0,0.d0)
               do l=1,3
                  tmp = tmp + zl(l,m)*dconjg(zl(l+3,n))*yl(l)
               end do
               sig_ll_nhol = sig_ll_nhol - tmp*dconjg(zl(i,m))*zl(j+3,n)
     $              * mu*e2/ct2*gm3*cp0(abs(gm3),slm(m),slm(n))
            end do
         end do
      end if
c     Chargino contribution
      if (ic.eq.1) then
         do m=1,3
            sig_ll_nhol = sig_ll_nhol + e2/st2*yl(j)*gm2*mu*zv(i,m)
     $           * dconjg(zv(j,m))*cp0(abs(mu),abs(gm2),vm(m))
         end do
      end if
      sig_ll_nhol = - v2/sq2*sig_ll_nhol/16/pi/pi
      return
      end

      double complex function sig_ll_hol(i,j)
c     Chirally enhanced d-quark self energy, "holomorphic" terms, i.e.
c     propotional to v1 in the decoupling limit
      implicit double precision (a-h,o-z)
      double complex gm2,gm3
      double complex zv,zl
      double complex ls,ks,ds,es,us,ws
      double complex tmp
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/gmass/gm1,gm2,gm3
      common/vev/v1,v2
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/debug_4q/ih,ic,in,ig
      sig_ll_hol = (0.d0,0.d0)
      if (in.ne.1) return
c     Neutralino contribution
      do m=1,6
         do n=1,6
            tmp = (0.d0,0.d0)
            do l=1,3
               do k=1,3
                 tmp = tmp + zl(l,m)*dconjg(ls(l,k)*zl(k+3,n))
               end do
            end do
            sig_ll_hol = sig_ll_hol - tmp*dconjg(zl(i,m))*zl(j+3,n)
     $           *e2/ct2*gm3*cp0(abs(gm3),slm(m),slm(n))
         end do
      end do
      sig_ll_hol = - v1/sq2*sig_ll_hol/16/pi/pi
      return 
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Complete scalar left and right lepton self energies, calculated  c
c     at p^2=0 without the resummation of enhanced chiral corrections! c
c     Non-decoupling inifinite terms truncated                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Left lepton self-energy (proportional to P_L)                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function esl_sig_n(i,j)
c     Neutralino and down squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zv,zl,zn
      double complex vl_lln0,vr_lln0
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/neut/fnm(4),zn(4,4)
      esl_sig_n = (0.d0,0.d0)
      do k=1,6
         do l=1,4
            esl_sig_n = esl_sig_n - vl_lln0(i,k,l)
     $           * dconjg(vr_lln0(j,k,l))
     $           * fnm(l)*b0(0.d0,fnm(l),slm(k))
         end do
      end do
      return
      end

      double complex function esl_sig_c(i,j)
c     Chargino and up squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zv,zl,zpos,zneg
      double complex vl_lsnc0,vr_lsnc0
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      esl_sig_c = (0.d0,0.d0)
      do k=1,3
         do l=1,2
            esl_sig_c = esl_sig_c - vl_lsnc0(i,k,l)
     $           * dconjg(vr_lsnc0(j,k,l))*fcm(l)*b0(0.d0,fcm(l),vm(k))
        end do
      end do
      return
      end

      double complex function esl_sig(i,j)
      implicit double precision (a-h,o-z)
      double complex esl_sig_n,esl_sig_c
      logical cdiag,ndiag,tmp
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      common/nc_suppress/eps_d,eps_u,acc
      external init_nc_diag
      esl_sig = (0.d0,0.d0)
      if ((in.eq.1).or.(ic.eq.1)) then
         eps = eps_d
         eps_d = acc
         tmp = cdiag().or.ndiag()
      end if
      if (in.eq.1) esl_sig = esl_sig + esl_sig_n(i,j)
      if (ic.eq.1) esl_sig = esl_sig + esl_sig_c(i,j)
      if ((in.eq.1).or.(ic.eq.1)) then
         eps_d = eps
         tmp = cdiag().or.ndiag()
      end if
      esl_sig = esl_sig/16/pi/pi
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Right lepton self-energy (proportional to P_R)                    c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function esr_sig_n(i,j)
c     Neutralino and down squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zv,zl,zn
      double complex vl_lln0,vr_lln0
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/neut/fnm(4),zn(4,4)
      esr_sig_n = (0.d0,0.d0)
      do k=1,6
         do l=1,4
            esr_sig_n = esr_sig_n - vr_lln0(i,k,l)
     $           * dconjg(vl_lln0(j,k,l))
     $           * fnm(l)*b0(0.d0,fnm(l),slm(k))
        end do
      end do
      return
      end

      double complex function esr_sig_c(i,j)
c     Chargino and up squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zv,zl,zpos,zneg
      double complex vl_lsnc0,vr_lsnc0
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      esr_sig_c = (0.d0,0.d0)
      do k=1,3
         do l=1,2
            esr_sig_c = esr_sig_c - vr_lsnc0(i,k,l)
     $           * dconjg(vl_lsnc0(j,k,l))*fcm(l)*b0(0.d0,fcm(l),vm(k))
        end do
      end do
      return
      end

      double complex function esr_sig(i,j)
      implicit double precision (a-h,o-z)
      double complex esr_sig_n,esr_sig_c
      logical cdiag,ndiag,tmp
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      common/nc_suppress/eps_d,eps_u,acc
      external init_nc_diag
      esr_sig = (0.d0,0.d0)
      if ((in.eq.1).or.(ic.eq.1)) then
         eps = eps_d
         eps_d = acc
         tmp = cdiag().or.ndiag()
      end if
      if (in.eq.1) esr_sig = esr_sig + esr_sig_n(i,j)
      if (ic.eq.1) esr_sig = esr_sig + esr_sig_c(i,j)
      if ((in.eq.1).or.(ic.eq.1)) then
         eps_d = eps
         tmp = cdiag().or.ndiag()
      end if
      esr_sig = esr_sig/16/pi/pi
      return
      end

