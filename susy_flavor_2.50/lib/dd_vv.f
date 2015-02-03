c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor
      
c     FILENAME: DD_VV.F
c     Released: 25:03:1995 (J.R.)
c     Revised: 18.10.2013 (J.R.)
c     Bare -> Physical CKM replaced in gauge and Higgs diagrams
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expression for the coefficient of the    c
c     Hamiltonian for the B->pi+v+vbar and the K->pi+v+vbar       c
c     processes (at the quark level).                             c
c     General form of the Hamiltonian is:                         c
c     iH = A^V_L H^V_L + A^V_R H^V_R                              c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^V_L                                            c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dd_vv_l_sm(i,j,k,l)
c     Gauge contributions
      implicit double precision (a-h,o-z)
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      dd_vv_l_sm = (0.d0,0.d0)
      if (k.ne.l) return
      xt =  umu(3)*umu(3)/wm2
c     Penguins + Boxes
      dd_vv_l_sm = - e2*e2/st2/st2/wm2*ckm_phys(3,i)
     $     * dconjg(ckm_phys(3,j))*(g_dv(xt) - 4*f_dv(xt))
      return
      end

      double precision function f_dv(x)
      implicit double precision (a-h,o-z)
      f_dv = x/(1 - x)/4*(1 + log(x)/(1 - x))
      return 
      end

      double precision function g_dv(x)
      implicit double precision (a-h,o-z)
      g_dv = x/(1 - x)/8*(6 - x + (2 + 3*x)*log(x)/(1 - x))
      return 
      end

      double complex function dd_vv_l_h(i,j,k,l)
c     Higgs contributions
      implicit double precision (a-h,o-z)
      double complex yh_eff_l
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      dd_vv_l_h = (0.d0,0.d0)
      if (k.ne.l) return
      do m=1,3
c     Penguins
         dd_vv_l_h = dd_vv_l_h + yh_eff_l(i,m,1)*dconjg(yh_eff_l(j,m,1))
     $        * umu(m)**2*cp0(cm(1),umu(m),umu(m))
      end do
      dd_vv_l_h = e2/st2/wm2/4*dd_vv_l_h
      return
      end

      double complex function dd_vv_l_c(i,j,k,l)
c     chargino contributions
      implicit double precision (a-h,o-z)
      double complex vl_duc,v_nlc,zpos,zneg,zu,zd,zv,zl,ztmp
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dd_vv_l_c = (0.d0,0.d0)
c     Boxes
      do m=1,2
         do n=1,2
            do mm=1,6
               do nn=1,6
                  dd_vv_l_c = dd_vv_l_c + vl_duc(i,mm,m)*v_nlc(l,nn,m)
     $                 * dconjg(vl_duc(j,mm,n)*v_nlc(k,nn,n))/2
     $                 *fcm(m)*fcm(n)*dp0(fcm(m),fcm(n),sum(mm),slm(nn))
               end do
            end do
         end do
      end do
      if (k.ne.l) return
c     Penguins
      do m=1,2
        do n=1,2
          do mm=1,6
            dd_vv_l_c = dd_vv_l_c + e2/8/st2/wm2
     1          * vl_duc(i,mm,m)*dconjg(vl_duc(j,mm,n))
     2          *(zpos(1,n)*dconjg(zpos(1,m))*cp1(sum(mm),fcm(m),fcm(n))
     3          - 2*zneg(1,m)*dconjg(zneg(1,n))*fcm(m)*fcm(n)
     4          * cp0(sum(mm),fcm(m),fcm(n)))
          end do
        end do
      end do
      do mm=1,6
        do nn=1,6
          ztmp = (0.d0,0.d0)
          do n=1,3
            ztmp = ztmp + zu(n,mm)*dconjg(zu(n,nn))
          end do
          do m=1,2
            dd_vv_l_c = dd_vv_l_c - e2/8/st2/wm2*ztmp*vl_duc(i,mm,m)
     1               *dconjg(vl_duc(j,nn,m))*cp1(fcm(m),sum(mm),sum(nn))
          end do
        end do
      end do
      return
      end

      double complex function dd_vv_l_n(i,j,k,l)
c     neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vl_ddn,vl_nnz,zn,zu,zd,zv,zl
      double complex ztmp,tmp1,tmp2
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/neut/fnm(4),zn(4,4)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dd_vv_l_n = (0.d0,0.d0)
c     Boxes
      do m=1,4
         tmp1 = zn(1,m)*st - zn(2,m)*ct
         do n=1,4
            tmp2 = zn(1,n)*st - zn(2,n)*ct
            do mm=1,6
               do nn=1,3
                  dd_vv_l_n = dd_vv_l_n + e2/4/sct2*zv(k,nn)
     $                 * vl_ddn(i,mm,m)*dconjg(zv(l,nn)*vl_ddn(j,mm,n))
     $                 * ( - dconjg(tmp1)*tmp2/2
     $                 * dp1(fnm(m),fnm(n),sdm(mm),vm(nn))
     $                 + tmp1*dconjg(tmp2)*fnm(m)*fnm(n)
     $                 * dp0(fnm(m),fnm(n),sdm(mm),vm(nn)))
               end do
            end do
         end do
      end do
      if (k.ne.l) return
c     Penguins
      do m=1,4
         do n=1,4
            do mm=1,6
               dd_vv_l_n = dd_vv_l_n - e2/4/st2/wm2
     $              * vl_ddn(i,mm,m)*dconjg(vl_ddn(j,mm,n))
     $              * (dconjg(vl_nnz(n,m))*cp1(sdm(mm),fnm(m),fnm(n))
     $              + 2*vl_nnz(n,m)*fnm(m)*fnm(n)
     $              * cp0(sdm(mm),fnm(m),fnm(n)))
            end do
         end do
      end do
      do mm=1,6
         do nn=1,6
            ztmp = (0.d0,0.d0)
            do n=1,3
               ztmp = ztmp + zd(n+3,nn)*dconjg(zd(n+3,mm))
            end do
            do m=1,4
               dd_vv_l_n = dd_vv_l_n - e2/8/st2/wm2*ztmp*vl_ddn(i,mm,m)
     $              *dconjg(vl_ddn(j,nn,m))*cp1(fnm(m),sdm(mm),sdm(nn))
            end do
         end do
      end do
      return
      end

      double complex function dd_vv_l_gl(i,j,k,l)
c     gluino contributions
      implicit double precision (a-h,o-z)
      double complex zu0,zd0,gm2,gm3,ztmp
      double complex zd
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      dd_vv_l_gl = (0.d0,0.d0)
      if (k.ne.l) return
      do mm=1,6
         do nn=1,6
            ztmp = (0.d0,0.d0)
            do m=1,3
               ztmp = ztmp + zd(m+3,nn)*dconjg(zd(m+3,mm))
            end do
            dd_vv_l_gl = dd_vv_l_gl - e2*g3d*g3d/3/st2/wm2*ztmp
     $           * zd(i,mm)*dconjg(zd(j,nn))*cp1(gm1,sdm(mm),sdm(nn))
         end do
      end do
      return
      end

      double complex function dd_vv_l(i,j,k,l)
c     Full A^V_L formfactor
      implicit double precision (a-h,o-z)
      double complex dd_vv_l_sm,dd_vv_l_h,dd_vv_l_c,dd_vv_l_n,dd_vv_l_gl
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      dd_vv_l = (0.d0,0.d0)
      if (ih.eq.1) then
         dd_vv_l = dd_vv_l + dd_vv_l_sm(i,j,k,l)
         dd_vv_l = dd_vv_l + dd_vv_l_h(i,j,k,l)
      end if
      if (in.eq.1) dd_vv_l = dd_vv_l + dd_vv_l_n(i,j,k,l)
      if (ic.eq.1) dd_vv_l = dd_vv_l + dd_vv_l_c(i,j,k,l)
      if (ig.eq.1) dd_vv_l = dd_vv_l + dd_vv_l_gl(i,j,k,l)
      dd_vv_l = dd_vv_l/16/pi/pi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^V_R                                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dd_vv_r_sm(i,j,k,l)
c     Gauge contributions
      implicit double precision (a-h,o-z)
      dd_vv_r_sm = (0.d0,0.d0)
      return
      end

      double complex function dd_vv_r_h(i,j,k,l)
c     Higgs contributions
      implicit double precision (a-h,o-z)
      double complex yh_eff_r
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      dd_vv_r_h = (0.d0,0.d0)
      if (k.ne.l) return
      do m=1,3
         do n=1,2
c     Penguins
            dd_vv_r_h = dd_vv_r_h + yh_eff_r(i,m,n)
     $           * dconjg(yh_eff_r(j,m,n))
     $           * umu(m)**2*cp0(cm(n),umu(m),umu(m))
         end do
      end do
      dd_vv_r_h = - e2/st2/wm2/4*dd_vv_r_h
      return
      end

      double complex function dd_vv_r_c(i,j,k,l)
c     chargino contributions
      implicit double precision (a-h,o-z)
      double complex vr_duc,v_nlc,zpos,zneg,zu,zd,zv,zl,ztmp
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dd_vv_r_c = (0.d0,0.d0)
c     Boxes
      do m=1,2
         do n=1,2
            do mm=1,6
               do nn=1,6
                  dd_vv_r_c = dd_vv_r_c - vr_duc(i,mm,m)*v_nlc(l,nn,m)
     $                 * dconjg(vr_duc(j,mm,n)*v_nlc(k,nn,n))/4
     $                 * dp1(fcm(m),fcm(n),sum(mm),slm(nn))
               end do
            end do
         end do
      end do
      if (k.ne.l) return
c     Penguins
      do m=1,2
         do n=1,2
            do mm=1,6
               dd_vv_r_c = dd_vv_r_c + e2/8/st2/wm2
     $              * vr_duc(i,mm,m)*dconjg(vr_duc(j,mm,n))
     $              * (zneg(1,m)*dconjg(zneg(1,n))
     $              * cp1(sum(mm),fcm(m),fcm(n))
     $              -  2*zpos(1,n)*dconjg(zpos(1,m))*fcm(m)*fcm(n)
     $              *  cp0(sum(mm),fcm(m),fcm(n)))
          end do
        end do
      end do
      do mm=1,6
         do nn=1,6
            ztmp = (0.d0,0.d0)
            do n=1,3
               ztmp = ztmp - zu(n+3,mm)*dconjg(zu(n+3,nn))
            end do
            do m=1,2
               dd_vv_r_c = dd_vv_r_c - e2/8/st2/wm2*ztmp*vr_duc(i,mm,m)
     $              * dconjg(vr_duc(j,nn,m))*cp1(fcm(m),sum(mm),sum(nn))
            end do
         end do
      end do
      return
      end

      double complex function dd_vv_r_n(i,j,k,l)
c     neutralino contributions
      implicit double precision (a-h,o-z)
      double complex vr_ddn,vl_nnz,zn,zu,zd,zv,zl
      double complex ztmp,tmp1,tmp2
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/neut/fnm(4),zn(4,4)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dd_vv_r_n = (0.d0,0.d0)
c     Boxes
      do m=1,4
         tmp1 = zn(1,m)*st - zn(2,m)*ct
         do n=1,4
            tmp2 = zn(1,n)*st - zn(2,n)*ct
            do mm=1,6
               do nn=1,3
                  dd_vv_r_n = dd_vv_r_n - e2/4/sct2*zv(k,nn)
     $                 * vr_ddn(i,mm,m)*dconjg(zv(l,nn)*vr_ddn(j,mm,n))
     $                 * (tmp1*dconjg(tmp2)/2
     $                 * dp1(fnm(m),fnm(n),sdm(mm),vm(nn))
     $                 - dconjg(tmp1)*tmp2*fnm(m)*fnm(n)
     $                 * dp0(fnm(m),fnm(n),sdm(mm),vm(nn)))
               end do
            end do
         end do
      end do
      if (k.ne.l) return
c     Penguins
      do m=1,4
         do n=1,4
            do mm=1,6
               dd_vv_r_n = dd_vv_r_n + e2/4/st2/wm2
     $              * vr_ddn(i,mm,m)*dconjg(vr_ddn(j,mm,n))
     $              * (vl_nnz(n,m)*cp1(sdm(mm),fnm(m),fnm(n))
     $              + 2*dconjg(vl_nnz(n,m))*fnm(m)*fnm(n)
     $              * cp0(sdm(mm),fnm(m),fnm(n)))
            end do
         end do
      end do
      do mm=1,6
         do nn=1,6
            ztmp = (0.d0,0.d0)
            do n=1,3
               ztmp = ztmp - zd(n,nn)*dconjg(zd(n,mm))
            end do
            do m=1,4
               dd_vv_r_n = dd_vv_r_n - e2/8/st2/wm2*ztmp*vr_ddn(i,mm,m)
     $              * dconjg(vr_ddn(j,nn,m))*cp1(fnm(m),sdm(mm),sdm(nn))
            end do
         end do
      end do
      return
      end

      double complex function dd_vv_r_gl(i,j,k,l)
c     gluino contributions
      implicit double precision (a-h,o-z)
      double complex zu0,zd0,ztmp,gm2,gm3
      double complex zd
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      dd_vv_r_gl = (0.d0,0.d0)
      if (k.ne.l) return
      do mm=1,6
         do nn=1,6
            ztmp = (0.d0,0.d0)
            do m=1,3
               ztmp = ztmp + zd(m,nn)*dconjg(zd(m,mm))
            end do
            dd_vv_r_gl = dd_vv_r_gl + e2*g3d*g3d/3/st2/wm2*ztmp
     $           * zd(i+3,mm)*dconjg(zd(j+3,nn))
     $           * cp1(gm1,sdm(mm),sdm(nn))
         end do
      end do
      return
      end
      
      double complex function dd_vv_r(i,j,k,l)
c     Full A^V_R formfactor
      implicit double precision (a-h,o-z)
      double complex dd_vv_r_sm,dd_vv_r_h,dd_vv_r_c,dd_vv_r_n,dd_vv_r_gl
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ig
      dd_vv_r = (0.d0,0.d0)
      if (ih.eq.1) then
         dd_vv_r = dd_vv_r + dd_vv_r_sm(i,j,k,l)
         dd_vv_r = dd_vv_r + dd_vv_r_h(i,j,k,l)
      end if
      if (in.eq.1) dd_vv_r = dd_vv_r + dd_vv_r_n(i,j,k,l)
      if (ic.eq.1) dd_vv_r = dd_vv_r + dd_vv_r_c(i,j,k,l)
      if (ig.eq.1) dd_vv_r = dd_vv_r + dd_vv_r_gl(i,j,k,l)
      dd_vv_r = dd_vv_r/16/pi/pi
      return
      end



