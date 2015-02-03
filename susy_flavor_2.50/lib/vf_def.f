c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM} 
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor
 
c     FILENAME: VF_DEF.F
c     Released: 5: 1:1993(J.R.)
c     Revised: 28: 1:1993 (P.Ch.)
c     Vertices vl(r)_lsnc (lepton-sneutrino-chargino) added
c     Revised: 15: 2:2010 (J.R)
c     Resummed quark vertices added
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expressions for fermion vertices              c
c     Convention is such that the Feynman rule for the vertex is       c
c     (-i)*(expression calculated below)*(projector P_l or P_R)        c
c     Compare with the paper: J.Rosiek@Phys.Rev.D41(1990)p.3464;       c
c     erratum, hep-ph/9511250                                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
cccccccccccccccccccccccccccccccccc
c     SUSY vertices              c
cccccccccccccccccccccccccccccccccc

      double complex function v_nnn(i,j,k)
c     Neutrino-sneutrino-neutralino vertex
      implicit double precision (a-h,o-z)
      double complex zv,zl,zn
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/neut/fnm(4),zn(4,4)
      v_nnn = -e/sct/sq2*dconjg(zv(i,j))*(zn(1,k)*st - zn(2,k)*ct)
      return
      end
 
      double complex function vl_lln0(i,j,k)
c     Lepton-slepton-neutralino left vertex
      implicit double precision (a-h,o-z)
      double complex zv,zl,zn
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/neut/fnm(4),zn(4,4)
      common/yukawa/yl(3),yu(3),yd(3)
      vl_lln0 = -e/sct/sq2*zl(i,j)*(zn(1,k)*st + zn(2,k)*ct) 
     $     - dconjg(yl(i))*zl(i+3,j)*zn(3,k)
      return
      end
 
      double complex function vr_lln0(i,j,k)
c     Lepton-slepton-neutralino right vertex
      implicit double precision (a-h,o-z)
      double complex zv,zl,zn
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/neut/fnm(4),zn(4,4)
      common/yukawa/yl(3),yu(3),yd(3)
      vr_lln0 = sq2*e/ct*zl(i+3,j)*dconjg(zn(1,k))
     $     - yl(i)*zl(i,j)*dconjg(zn(3,k))
      return
      end
 
      double complex function vl_lsnc0(i,j,k)
c     Lepton-sneutrino-chargino left vertex
      implicit double precision (a-h,o-z)
      double complex zv,zl,zpos,zneg
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      vl_lsnc0 = e/st*dconjg(zv(i,j))*zpos(1,k)
      return
      end
 
      double complex function vr_lsnc0(i,j,k)
c     Lepton-sneutrino-chargino right vertex
      implicit double precision (a-h,o-z)
      double complex zv,zl,zpos,zneg
      double complex yl,yu,yd
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/yukawa/yl(3),yu(3),yd(3)
      vr_lsnc0 = yl(i)*dconjg(zv(i,j)*zneg(2,k))
      return
      end
 
      double complex function v_nlc(i,j,k)
c     Neutrino-slepton-chargino left vertex
      implicit double precision (a-h,o-z)
      double complex zv,zl,zpos,zneg
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/yukawa/yl(3),yu(3),yd(3)
      v_nlc = e/st*zl(i,j)*zneg(1,k) + dconjg(yl(i))*zl(i+3,j)*zneg(2,k)
      return
      end
 
      double complex function vl_uun0(i,j,k)
c     Up quark-up squark-neutralino left vertex
      implicit double precision (a-h,o-z)
      double complex zu,zd,zn
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      common/yukawa/yl(3),yu(3),yd(3)
      vl_uun0 = - e/sct/sq2*dconjg(zu(i,j))*(zn(1,k)*st/3 + zn(2,k)*ct)
     $     - dconjg(yu(i)*zu(i+3,j))*zn(4,k)
      return
      end

      double complex function vr_uun0(i,j,k)
c     Up quark-up squark-neutralino right vertex
      implicit double precision (a-h,o-z)
      double complex zu,zd,zn
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      common/yukawa/yl(3),yu(3),yd(3)
      vr_uun0 = 2*sq2*e/3.d0/ct*dconjg(zu(i+3,j)*zn(1,k))
     $     - yu(i)*dconjg(zu(i,j)*zn(4,k))
      return
      end
 
      double complex function vl_ddn0(i,j,k)
c     Down quark-down squark-neutralino left vertex
      implicit double precision (a-h,o-z)
      double complex zu,zd,zn
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      common/yukawa/yl(3),yu(3),yd(3)
      vl_ddn0 =  e/sct/sq2*zd(i,j)*(zn(1,k)*st/3 - zn(2,k)*ct)
     1       - dconjg(yd(i))*zd(i+3,j)*zn(3,k)
      return
      end
 
      double complex function vr_ddn0(i,j,k)
c     Down quark-down squark-neutralino right vertex
      implicit double precision (a-h,o-z)
      double complex zu,zd,zn
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      common/yukawa/yl(3),yu(3),yd(3)
      vr_ddn0 = sq2*e/3/ct*zd(i+3,j)*dconjg(zn(1,k))
     $     - yd(i)*zd(i,j)*dconjg(zn(3,k))
      return
      end
 
      double complex function vl_duc0(i,j,k)
c     Down quark-up squark-chargino left vertex
      implicit double precision (a-h,o-z)
      double complex zu,zd,zpos,zneg
      double complex ckm
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/yukawa/yl(3),yu(3),yd(3)
      common/km_mat/ckm(3,3)
      vl_duc0 = (0.d0,0.d0)
      do l=1,3
        vl_duc0 = vl_duc0 
     $        + dconjg(ckm(i,l))*(e/st*dconjg(zu(l,j))*zpos(1,k)
     $        - dconjg(yu(l)*zu(l+3,j))*zpos(2,k))
      end do
      return
      end
 
      double complex function vr_duc0(i,j,k)
c     Down quark-up squark-chargino right vertex
      implicit double precision (a-h,o-z)
      double complex zu,zd,zpos,zneg
      double complex ckm
      double complex yl,yu,yd
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/yukawa/yl(3),yu(3),yd(3)
      common/km_mat/ckm(3,3)
      vr_duc0 = (0.d0,0.d0)
      do l=1,3
         vr_duc0 = vr_duc0 + yd(i)*dconjg(zu(l,j)*zneg(2,k)*ckm(i,l))
      end do
      return
      end
 
      double complex function vl_udc0(i,j,k)
c     Up quark-down squark-chargino left vertex
      implicit double precision (a-h,o-z)
      double complex zu,zd,zpos,zneg
      double complex ckm
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/yukawa/yl(3),yu(3),yd(3)
      common/km_mat/ckm(3,3)
      vl_udc0 = (0.d0,0.d0)
      do l=1,3
         vl_udc0 = vl_udc0 - (e/st*zd(l,j)*zneg(1,k) 
     $        - dconjg(yd(l))*zd(l+3,j)*zneg(2,k))*ckm(l,i)
      end do
      return
      end
 
      double complex function vr_udc0(i,j,k)
c     Up quark-down squark-chargino right vertex
      implicit double precision (a-h,o-z)
      double complex zu,zd,zpos,zneg
      double complex ckm
      double complex yl,yu,yd
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/yukawa/yl(3),yu(3),yd(3)
      common/km_mat/ckm(3,3)
      vr_udc0 = (0.d0,0.d0)
      do l=1,3
         vr_udc0 = vr_udc0 + yu(i)*zd(l,j)*dconjg(zpos(2,k))*ckm(l,i)
      end do
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Effective vertices (chirally enhanced effects resummed)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function vl_uun(i,j,k)
c     Up quark-up squark-neutralino left vertex, resummed
      implicit double precision (a-h,o-z)
      double complex vl_uun0
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      common/resum_level/il,nmax,errin,errout
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      if (il.eq.0) then
         vl_uun = vl_uun0(i,j,k)
      else
         vl_uun = (0.d0,0.d0)
         do l=1,3
            vl_uun = vl_uun + uul(l,i)*vl_uun0(l,j,k)
         end do
      end if
      return
      end

      double complex function vr_uun(i,j,k)
c     UP quark-up squark-neutralino right vertex, resummed
      implicit double precision (a-h,o-z)
      double complex vr_uun0
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      common/resum_level/il,nmax,errin,errout
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      if (il.eq.0) then
         vr_uun = vr_uun0(i,j,k)
      else
         vr_uun = (0.d0,0.d0)
         do l=1,3
            vr_uun = vr_uun + uur(l,i)*vr_uun0(l,j,k)
         end do
      end if
      return
      end
 
      double complex function vl_ddn(i,j,k)
c     Down quark-down squark-neutralino left vertex, resummed
      implicit double precision (a-h,o-z)
      double complex vl_ddn0
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      common/resum_level/il,nmax,errin,errout
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      if (il.eq.0) then
         vl_ddn = vl_ddn0(i,j,k)
      else
         vl_ddn = (0.d0,0.d0)
         do l=1,3
            vl_ddn = vl_ddn + udl(l,i)*vl_ddn0(l,j,k)
         end do
      end if
      return
      end

      double complex function vr_ddn(i,j,k)
c     Down quark-down squark-neutralino right vertex, resummed
      implicit double precision (a-h,o-z)
      double complex vr_ddn0
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      common/resum_level/il,nmax,errin,errout
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      if (il.eq.0) then
         vr_ddn = vr_ddn0(i,j,k)
      else
         vr_ddn = (0.d0,0.d0)
         do l=1,3
            vr_ddn = vr_ddn + udr(l,i)*vr_ddn0(l,j,k)
         end do
      end if
      return
      end

      double complex function vl_duc(i,j,k)
c     Down quark-up squark-chargino left vertex, resummed
      implicit double precision (a-h,o-z)
      double complex vl_duc0
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      common/resum_level/il,nmax,errin,errout
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      if (il.eq.0) then
         vl_duc = vl_duc0(i,j,k)
      else
         vl_duc = (0.d0,0.d0)
         do l=1,3
            vl_duc = vl_duc + udl(l,i)*vl_duc0(l,j,k)
         end do
      end if
      return
      end
 
      double complex function vr_duc(i,j,k)
c     Down quark-up squark-chargino right vertex, resummed
      implicit double precision (a-h,o-z)
      double complex vr_duc0
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      common/resum_level/il,nmax,errin,errout
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      if (il.eq.0) then
         vr_duc = vr_duc0(i,j,k)
      else
         vr_duc = (0.d0,0.d0)
         do l=1,3
            vr_duc = vr_duc + udr(l,i)*vr_duc0(l,j,k)
         end do
      end if
      return
      end
 
      double complex function vl_udc(i,j,k)
c     Up quark-down squark-chargino left vertex, resummed
      implicit double precision (a-h,o-z)
      double complex vl_udc0
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      common/resum_level/il,nmax,errin,errout
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      if (il.eq.0) then
         vl_udc = vl_udc0(i,j,k)
      else
         vl_udc = (0.d0,0.d0)
         do l=1,3
            vl_udc = vl_udc + uul(l,i)*vl_udc0(l,j,k)
         end do
      end if
      return
      end
 
      double complex function vr_udc(i,j,k)
c     Up quark-down squark-chargino right vertex, resummed
      implicit double precision (a-h,o-z)
      double complex vr_udc0
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      common/resum_level/il,nmax,errin,errout
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      if (il.eq.0) then
         vr_udc = vr_udc0(i,j,k)
      else
         vr_udc = (0.d0,0.d0)
         do l=1,3
            vr_udc = vr_udc + uur(l,i)*vr_udc0(l,j,k)
         end do
      end if
      return
      end

      double complex function zd(i,j)
c     Effective down squark mixing matrix, resummed
      implicit double precision (a-h,o-z)
      double complex zd0,zu0,zd_eff,zu_eff
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/sqmass_eff/sum_eff(6),sdm_eff(6),zu_eff(6,6),zd_eff(6,6)
      common/resum_level/il,nmax,errin,errout
      if (il.eq.0) then
         zd = zd0(i,j)
      else
         zd = zd_eff(i,j)
      end if
      return
      end
 
      double complex function zu(i,j)
c     Effective up squark mixing matrix, resummed
      implicit double precision (a-h,o-z)
      double complex zd0,zu0,zd_eff,zu_eff
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/sqmass_eff/sum_eff(6),sdm_eff(6),zu_eff(6,6),zd_eff(6,6)
      common/resum_level/il,nmax,errin,errout
      if (il.eq.0) then
         zu = zu0(i,j)
      else
         zu = zu_eff(i,j)
      end if
      return
      end

      double complex function vl_lln(i,j,k)
c     Lepton-slepton-neutralino left vertex, resummed
      implicit double precision (a-h,o-z)
      double complex vl_lln0
      double complex ull,ulr
      common/resum_level/il,nmax,errin,errout
      common/lepton_switch/ull(3,3),ulr(3,3)
      if (il.eq.0) then
         vl_lln = vl_lln0(i,j,k)
      else
         vl_lln = (0.d0,0.d0)
         do l=1,3
            vl_lln = vl_lln + ull(l,i)*vl_lln0(l,j,k)
         end do
      end if
      return
      end

      double complex function vr_lln(i,j,k)
c     Lepton-slepton-neutralino right vertex, resummed
      implicit double precision (a-h,o-z)
      double complex vr_lln0
      double complex ull,ulr
      common/resum_level/il,nmax,errin,errout
      common/lepton_switch/ull(3,3),ulr(3,3)
      if (il.eq.0) then
         vr_lln = vr_lln0(i,j,k)
      else
         vr_lln = (0.d0,0.d0)
         do l=1,3
            vr_lln = vr_lln + ulr(l,i)*vr_lln0(l,j,k)
         end do
      end if
      return
      end

      double complex function vl_lsnc(i,j,k)
c     Lepton-sneutrino-chargino left vertex, resummed
      implicit double precision (a-h,o-z)
      double complex vl_lsnc0
      double complex ull,ulr
      common/resum_level/il,nmax,errin,errout
      common/lepton_switch/ull(3,3),ulr(3,3)
      if (il.eq.0) then
         vl_lsnc = vl_lsnc0(i,j,k)
      else
         vl_lsnc = (0.d0,0.d0)
         do l=1,3
            vl_lsnc = vl_lsnc + ull(l,i)*vl_lsnc0(l,j,k)
         end do
      end if
      return
      end

      double complex function vr_lsnc(i,j,k)
c     Lepton-sneutrino-chargino right vertex, resummed
      implicit double precision (a-h,o-z)
      double complex vr_lsnc0
      double complex ull,ulr
      common/resum_level/il,nmax,errin,errout
      common/lepton_switch/ull(3,3),ulr(3,3)
      if (il.eq.0) then
         vr_lsnc = vr_lsnc0(i,j,k)
      else
         vr_lsnc = (0.d0,0.d0)
         do l=1,3
            vr_lsnc = vr_lsnc + ulr(l,i)*vr_lsnc0(l,j,k)
         end do
      end if
      return
      end
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Effective charged Higgs vertices
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function yh_eff_r(i,j,k)
c     right charged Higgs-quark effective Yukawa 
      implicit double precision (a-h,o-z)
      double complex yhl,yhr,ysu,ypu,ysd,ypd
      logical init_yukawa_eff
      common/yukawa_eff/yhl(3,3,2),yhr(3,3,2),ysu(3,3,2),ypu(3,3,2),
     $     ysd(3,3,2),ypd(3,3,2),init_yukawa_eff
      if (init_yukawa_eff) call yukawa_eff_init
      yh_eff_r =  yhr(i,j,k)
      return
      end

      double complex function yh_eff_l(i,j,k)
c     left charged Higgs-quark effective Yukawa 
      implicit double precision (a-h,o-z)
      double complex yhl,yhr,ysu,ypu,ysd,ypd
      logical init_yukawa_eff
      common/yukawa_eff/yhl(3,3,2),yhr(3,3,2),ysu(3,3,2),ypu(3,3,2),
     $     ysd(3,3,2),ypd(3,3,2),init_yukawa_eff
      if (init_yukawa_eff) call yukawa_eff_init
      yh_eff_l =  yhl(i,j,k)
      return
      end

      double complex function yhl_eff_r(i,j,k)
c     right charged Higgs-lepton effective Yukawa 
      implicit double precision (a-h,o-z)
      double complex yhlr,ysl,ypl
      logical init_yukawa_l
      common/yukawa_lept/yhlr(3,3,2),ysl(3,3,2),ypl(3,3,2),init_yukawa_l
      if (init_yukawa_l) call yukawa_eff_init
      yhl_eff_r =  yhlr(i,j,k)
      return
      end

