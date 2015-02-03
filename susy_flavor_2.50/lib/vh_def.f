c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM} 
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor
 
c     FILENAME: VH_DEF.F
c     Revised: 25: 4:1996 (J.R.)
c     Vertices UDH,UUHH,DDHH added
c     Revised:  9: 6:2013 (J.R.)
c     Improved compatibility with notation of hep-ph/9511250
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains definitions of more complicated mixing        c
c     matrices and expressions for Higgs boson vertices                c
c     Compare with the paper: J.Rosiek@Phys.Rev.D41(1990)p.3464;       c
c     erratum, hep-ph/9511250                                          c
c     Common factors like e/2/sct etc. factorized from expressions     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      double precision function ah(i,j)
c     AH mixing matrix definition
      implicit double precision (a-h,o-z)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      ah = zh(1,i)*zh(1,j) - zh(2,i)*zh(2,j)
      return
      end
 
      double precision function ar(i,j)
c     AR mixing matrix definition
      implicit double precision (a-h,o-z)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      ar = zr(1,i)*zr(1,j) - zr(2,i)*zr(2,j)
      return
      end
 
      double precision function am(i,j)
c     AM mixing matrix definition
      implicit double precision (a-h,o-z)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      am = zr(1,i)*zh(1,j) - zr(2,i)*zh(2,j)
      return
      end
 
      double precision function ap(i,j)
c     AP mixing matrix definition
      implicit double precision (a-h,o-z)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      ap = zr(1,i)*zh(2,j) + zr(2,i)*zh(1,j)
      return
      end
 
      double precision function br(i)
c     BR mixing matrix definition
      implicit double precision (a-h,o-z)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/vev/v1,v2
      br = v1*zr(1,i) - v2*zr(2,i)
      return
      end
 
      double precision function cr(i)
c     CR mixing matrix definition - v(i)*zr(i,j)
      implicit double precision (a-h,o-z)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/vev/v1,v2
      cr = v1*zr(1,i) + v2*zr(2,i)
      return
      end
 
      double complex function vl_ccs(i,j,k)
c     Scalar - chargino - chargino left vertex
c     incoming chi_i^+, outgoing chi_j^+
      implicit double precision (a-h,o-z)
      double complex zpos,zneg
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      vl_ccs = - e/st/sq2*(zr(1,k)*zpos(1,i)*zneg(2,j) 
     $     + zr(2,k)*zpos(2,i)*zneg(1,j))
      return
      end
 
      double complex function vr_ccs(i,j,k)
c     Scalar - chargino - chargino right vertex
c     incoming chi_i^+, outgoing chi_j^+
      implicit double precision (a-h,o-z)
      double complex zpos,zneg
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      vr_ccs = - e/st/sq2*dconjg(zr(1,k)*zpos(1,j)*zneg(2,i) 
     $     + zr(2,k)*zpos(2,j)*zneg(1,i))
      return
      end
 
      double complex function vl_nns(i,j,k)
c     Scalar - neutralino -neutralino left vertex
      implicit double precision (a-h,o-z)
      double complex zn
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/neut/fnm(4),zn(4,4)
      vl_nns = e/2/sct*((zr(1,k)*zn(3,i) - zr(2,k)*zn(4,i))
     $     *(zn(1,j)*st - zn(2,j)*ct)
     $     + (zr(1,k)*zn(3,j) - zr(2,k)*zn(4,j))
     $     *(zn(1,i)*st - zn(2,i)*ct))
      return
      end

      double complex function vr_nns(i,j,k)
c     Scalar - neutralino -neutralino left vertex
      implicit double precision (a-h,o-z)
      double complex zn
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/neut/fnm(4),zn(4,4)
      vr_nns = e/2/sct*dconjg((zr(1,k)*zn(3,i) - zr(2,k)*zn(4,i))
     $     *(zn(1,j)*st - zn(2,j)*ct)
     $     + (zr(1,k)*zn(3,j) - zr(2,k)*zn(4,j))
     $     *(zn(1,i)*st - zn(2,i)*ct))
      return
      end

      double complex function v_dds(i,j,k)
c     Scalar - 2 d-squark vertex
      implicit double precision (a-h,o-z)
      double complex vdds,vi_dds
      logical init_dds
      common/rsd/vdds(6,6,2),init_dds
      if (init_dds) then
         do l=1,6
            do m=1,6
               do n=1,2
                  vdds(l,m,n) = vi_dds(l,m,n)
               end do
            end do
         end do
         init_dds = .false.
      end if
      v_dds = vdds(i,j,k)
      return
      end
 
      double complex function vi_dds(i,j,k)
c     Scalar - 2 d-squark vertex initialization
      implicit double precision (a-h,o-z)
      double complex zu,zd
      double complex h,ls,ks,ds,es,us,ws
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/hpar/hm1,hm2,hs,h
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/vev/v1,v2
      common/yukawa/yl(3),yu(3),yd(3)
      common/delta/del(6,6)
      vi_dds = e2/6/ct2*br(k)*del(i,j)
      do l=1,3
         vi_dds = vi_dds
     $        + e2/12/sct2*br(k)*(3 - 4*st2)*dconjg(zd(l,i))*zd(l,j)
         vi_dds = vi_dds - v1*abs(yd(l))**2*zr(1,k)
     $        *(dconjg(zd(l,i))*zd(l,j) + dconjg(zd(l+3,i))*zd(l+3,j))
         vi_dds = vi_dds - zr(2,k)/sq2
     $        *(dconjg(yd(l)*h*zd(l,i))*zd(l+3,j) 
     $        + yd(l)*h*zd(l,j)*dconjg(zd(l+3,i)))
         do m=1,3
            vi_dds = vi_dds - zr(1,k)/sq2
     $           *(dconjg(ds(l,m)*zd(m+3,i))*zd(l,j)
     $           + ds(l,m)*zd(m+3,j)*dconjg(zd(l,i)))
            vi_dds = vi_dds + zr(2,k)/sq2
     $           *(dconjg(es(l,m)*zd(m+3,i))*zd(l,j)
     $           + es(l,m)*zd(m+3,j)*dconjg(zd(l,i)))
         end do
      end do
      return
      end

      double complex function v_uus(i,j,k)
c     Scalar - 2 u-squark vertex
      implicit double precision (a-h,o-z)
      double complex vuus,vi_uus
      logical init_uus
      common/rsu/vuus(6,6,2),init_uus
      if (init_uus) then
         do l=1,6
            do m=1,6
               do n=1,2
                  vuus(l,m,n) = vi_uus(l,m,n)
               end do
            end do
         end do
         init_uus = .false.
      end if
      v_uus = vuus(i,j,k)
      return
      end
 
      double complex function vi_uus(i,j,k)
c     Scalar - 2 u-squark vertex initialization
      implicit double precision (a-h,o-z)
      double complex zu,zd
      double complex h,ls,ks,ds,es,us,ws
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/hpar/hm1,hm2,hs,h
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/vev/v1,v2
      common/yukawa/yl(3),yu(3),yd(3)
      common/delta/del(6,6)
      vi_uus = - e2/3/ct2*br(k)*del(i,j)
      do l=1,3
         vi_uus = vi_uus
     $        - e2/12/sct2*br(k)*(3 - 8*st2)*dconjg(zu(l,i))*zu(l,j)
         vi_uus = vi_uus - v2*abs(yu(l))**2*zr(2,k)
     $        * (dconjg(zu(l,i))*zu(l,j) + dconjg(zu(l+3,i))*zu(l+3,j))
         vi_uus = vi_uus + zr(1,k)/sq2
     $        * (dconjg(yu(l)*h*zu(l+3,i))*zu(l,j) 
     $        + yu(l)*h*zu(l+3,j)*dconjg(zu(l,i)))
         do m=1,3
            vi_uus = vi_uus + zr(2,k)/sq2
     $           *(dconjg(us(l,m)*zu(l,i))*zu(m+3,j)
     $           + us(l,m)*zu(l,j)*dconjg(zu(m+3,i)))
            vi_uus = vi_uus + zr(1,k)/sq2
     $           *(dconjg(ws(l,m)*zu(l,i))*zu(m+3,j)
     $           + ws(l,m)*zu(l,j)*dconjg(zu(m+3,i)))
         end do
      end do
      return
      end
 
      double complex function v_lls(i,j,k)
c     Scalar - 2 slepton vertex
      implicit double precision (a-h,o-z)
      double complex vlls,vi_lls
      logical init_lls
      common/rsl/vlls(6,6,2),init_lls
      if (init_lls) then
         do  l=1,6
            do  m=1,6
               do  n=1,2
                  vlls(l,m,n) = vi_lls(l,m,n)
               end do
            end do
         end do
         init_lls = .false.
      end if
      v_lls = vlls(i,j,k)
      return
      end
 
      double complex function vi_lls(i,j,k)
c     Scalar - 2 slepton vertex initialization
      implicit double precision (a-h,o-z)
      double complex zv,zl,h,ls,ks,ds,es,us,ws
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/hpar/hm1,hm2,hs,h
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/yukawa/yl(3),yu(3),yd(3)
      common/vev/v1,v2
      common/delta/del(6,6)
      vi_lls = e2/2/ct2*br(k)*del(i,j)
      do l=1,3
         vi_lls = vi_lls 
     $        + e2/4/sct2*br(k)*(1 - 4*st2)*dconjg(zl(l,i))*zl(l,j)
         vi_lls = vi_lls - v1*abs(yl(l))**2*zr(1,k)
     $        *(dconjg(zl(l,i))*zl(l,j) + dconjg(zl(l+3,i))*zl(l+3,j))
         vi_lls = vi_lls - zr(2,k)/sq2
     $        *(dconjg(yl(l)*h*zl(l,i))*zl(l+3,j) 
     $        + yl(l)*h*zl(l,j)*dconjg(zl(l+3,i)))
         do m=1,3
            vi_lls = vi_lls - zr(1,k)/sq2
     $           *(dconjg(ls(l,m)*zl(m+3,i))*zl(l,j)
     $           + ls(l,m)*zl(m+3,j)*dconjg(zl(l,i)))
            vi_lls = vi_lls + zr(2,k)/sq2
     $           *(dconjg(ks(l,m)*zl(m+3,i))*zl(l,j)
     $           + ks(l,m)*zl(m+3,j)*dconjg(zl(l,i)))
         end do
      end do
      return
      end

      double precision function v_sss(i,j,k)
c     3 - scalar vertex
      implicit double precision (a-h,o-z)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/vev/v1,v2
      v_sss = 3*v1*zr(1,i)*zr(1,j)*zr(1,k)
     $     + 3*v2*zr(2,i)*zr(2,j)*zr(2,k)
     $     - v1*(zr(2,i)*zr(2,j)*zr(1,k) + zr(2,i)*zr(1,j)*zr(2,k)
     $     + zr(1,i)*zr(2,j)*zr(2,k))
     $     - v2*(zr(1,i)*zr(1,j)*zr(2,k) + zr(1,i)*zr(2,j)*zr(1,k)
     $     + zr(2,i)*zr(1,j)*zr(1,k))
      return
      end
 
      double precision function v_hhs(i,j,k)
c     Scalar - 2 charged Higgs vertex
      implicit double precision (a-h,o-z)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/delta/d(6,6)
      v_hhs = e/2/st/ct2*ah(j,k)*br(i)
     $     + wm*(ap(i,j)*d(1,k) + ap(i,k)*d(1,j))
c     Leading second order correction added:
c     2       + shh_vert(i,j,k)
c     mt^4 term:  
c     - (um(3)*um(3)/pi/v2)**2*st/e/v2*zr(2,i)*zh(2,j)*zh(2,k)
c     *log(um(3)*um(3)/sum(1)/sum(6))
      return
      end
 
      double precision function v_ssss(i,j,k,l)
c     4 - scalar vertex
      implicit double precision (a-h,o-z)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      v_ssss = 3*zr(1,i)*zr(1,j)*zr(1,k)*zr(1,l)
     $       + 3*zr(2,i)*zr(2,j)*zr(2,k)*zr(2,l)
     $       - zr(1,i)*zr(1,j)*zr(2,k)*zr(2,l)
     $       - zr(1,i)*zr(2,j)*zr(1,k)*zr(2,l)
     $       - zr(1,i)*zr(2,j)*zr(2,k)*zr(1,l)
     $       - zr(2,i)*zr(1,j)*zr(1,k)*zr(2,l)
     $       - zr(2,i)*zr(1,j)*zr(2,k)*zr(1,l)
     $       - zr(2,i)*zr(2,j)*zr(1,k)*zr(1,l)
      return
      end
 
      double precision function v_hhss(i,j,k,l)
c     2 charged Higgses + 2 scalars vertex
      implicit double precision (a-h,o-z)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      v_hhss = ar(i,j)*ah(k,l)/ct2 + ap(i,k)*ap(j,l) + ap(i,l)*ap(j,k)
      return
      end
 
      double complex function v_llss(i,j,k,l)
c     2 scalar - 2 slepton vertex
      implicit double precision (a-h,o-z)
      double complex zv,zl
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/yukawa/yl(3),yu(3),yd(3)
      common/delta/del(6,6)
      v_llss = e2/2/ct2*ar(i,j)*del(k,l)
      do m=1,3
         v_llss = v_llss + e2/4/sct2*ar(i,j)*(1 - 4*st2)
     $        *dconjg(zl(m,k))*zl(m,l)
         v_llss = v_llss - yl(m)**2*zr(1,i)*zr(1,j)
     $        *(dconjg(zl(m,k))*zl(m,l) + dconjg(zl(m+3,k))*zl(m+3,l))
      end do
      return
      end
 
      double complex function v_uuss(i,j,k,l)
c     2 scalar - 2 U-squark vertex
      implicit double precision (a-h,o-z)
      double complex zu,zd
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/yukawa/yl(3),yu(3),yd(3)
      common/delta/del(6,6)
      v_uuss = - e2/3/ct2*ar(i,j)*del(k,l)
      do m=1,3
         v_uuss = v_uuss - e2/12/sct2*ar(i,j)*(3 - 8*st2)
     $        *dconjg(zu(m,k))*zu(m,l)
         v_uuss = v_uuss - yu(m)**2*zr(2,i)*zr(2,j)
     $        *(dconjg(zu(m,k))*zu(m,l) + dconjg(zu(m+3,k))*zu(m+3,l))
      end do
      return
      end
 
      double complex function v_ddss(i,j,k,l)
c     2 scalar - 2 D-squark vertex
      implicit double precision (a-h,o-z)
      double complex zu,zd
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/yukawa/yl(3),yu(3),yd(3)
      common/delta/del(6,6)
      v_ddss = e2/6/ct2*ar(i,j)*del(k,l)
      do m=1,3
         v_ddss = v_ddss + e2/12/sct2*ar(i,j)*(3 - 4*st2)
     $        *dconjg(zd(m,k))*zd(m,l)
         v_ddss = v_ddss - yd(m)**2*zr(1,i)*zr(1,j)
     $        *(dconjg(zd(m,k))*zd(m,l) + dconjg(zd(m+3,k))*zd(m+3,l))
      end do
      return
      end
 
      double complex function v_llp(i,j,k)
c     Pseudoscalar - 2 slepton vertex
      implicit double precision (a-h,o-z)
      double complex vllp,vi_llp
      logical init_llp
      common/psl/vllp(6,6,2),init_llp
      if (init_llp) then
         do l=1,6
            do m=1,6
               do n=1,2
                  vllp(l,m,n) = vi_llp(l,m,n)
               end do
            end do
         end do
         init_llp = .false.
      end if
      v_llp = vllp(i,j,k)
      return
      end
 
      double complex function v_uup(i,j,k)
c     Pseudoscalar - 2 u-squark vertex
      implicit double precision (a-h,o-z)
      double complex vuup,vi_uup
      logical init_uup
      common/psu/vuup(6,6,2),init_uup
      if (init_uup) then
         do l=1,6
            do m=1,6
               do n=1,2
                  vuup(l,m,n) = vi_uup(l,m,n)
               end do
            end do
         end do
         init_uup = .false.
      end if
      v_uup = vuup(i,j,k)
      return
      end
      
      double complex function v_ddp(i,j,k)
c     Pseudoscalar - 2 d-squark vertex
      implicit double precision (a-h,o-z)
      double complex vddp,vi_ddp
      logical init_ddp
      common/psd/vddp(6,6,2),init_ddp
      if (init_ddp) then
         do l=1,6
            do m=1,6
               do n=1,2
                  vddp(l,m,n) = vi_ddp(l,m,n)
               end do
            end do
         end do
         init_ddp = .false.
      end if
      v_ddp = vddp(i,j,k)
      return
      end
 
      double complex function vi_llp(i,j,k)
c     Pseudoscalar - 2 slepton vertex
      implicit double precision (a-h,o-z)
      double complex zv,zl,h,ls,ks,ds,es,us,ws
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/hpar/hm1,hm2,hs,h
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/yukawa/yl(3),yu(3),yd(3)
      vi_llp = (0.d0,0.d0)
      do l=1,3
         vi_llp = vi_llp - yl(l)/sq2*zh(2,k)
     $        *(dconjg(h*zl(l,i))*zl(l+3,j) 
     $        - h*zl(l,j)*dconjg(zl(l+3,i)))
         do m=1,3
            vi_llp = vi_llp - zh(1,k)*(dconjg(ls(l,m)*zl(m+3,i))*zl(l,j)
     $           - ls(l,m)*zl(m+3,j)*dconjg(zl(l,i)))/sq2
            vi_llp = vi_llp - zh(2,k)*(dconjg(ks(l,m)*zl(m+3,i))*zl(l,j)
     $           - ks(l,m)*zl(m+3,j)*dconjg(zl(l,i)))/sq2
         end do
      end do
      return
      end
 
      double complex function vi_ddp(i,j,k)
c     Pseudoscalar - 2 D-squark vertex
      implicit double precision (a-h,o-z)
      double complex zu,zd
      double complex h,ls,ks,ds,es,us,ws
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/hpar/hm1,hm2,hs,h
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/yukawa/yl(3),yu(3),yd(3)
      vi_ddp = (0.d0,0.d0)
      do l=1,3
         vi_ddp = vi_ddp - yd(l)/sq2*zh(2,k)
     $        *(dconjg(h*zd(l,i))*zd(l+3,j) 
     $        - h*zd(l,j)*dconjg(zd(l+3,i)))
         do m=1,3
            vi_ddp = vi_ddp - zh(1,k)*(dconjg(ds(l,m)*zd(m+3,i))*zd(l,j)
     $           - ds(l,m)*zd(m+3,j)*dconjg(zd(l,i)))/sq2
            vi_ddp = vi_ddp - zh(2,k)*(dconjg(es(l,m)*zd(m+3,i))*zd(l,j)
     $           - es(l,m)*zd(m+3,j)*dconjg(zd(l,i)))/sq2
         end do
      end do
      return
      end
 
      double complex function vi_uup(i,j,k)
c     Pseudoscalar - 2 U-squark vertex
      implicit double precision (a-h,o-z)
      double complex zu,zd
      double complex h,ls,ks,ds,es,us,ws
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/hpar/hm1,hm2,hs,h
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/yukawa/yl(3),yu(3),yd(3)
      vi_uup = (0.d0,0.d0)
      do l=1,3
         vi_uup = vi_uup + yu(l)/sq2*zh(1,k)
     $        *(dconjg(h*zu(l+3,i))*zu(l,j) 
     $        - h*zu(l+3,j)*dconjg(zu(l,i)))
         do m=1,3
            vi_uup = vi_uup + zh(2,k)*(dconjg(us(l,m)*zu(l,i))*zu(m+3,j)
     $           - us(l,m)*zu(l,j)*dconjg(zu(m+3,i)))/sq2
            vi_uup = vi_uup - zh(1,k)*(dconjg(ws(l,m)*zu(l,i))*zu(m+3,j)
     $           - ws(l,m)*zu(l,j)*dconjg(zu(m+3,i)))/sq2
         end do
      end do
      return
      end
 
      double complex function v_nnp(i,j,k)
c     Pseudoscalar - 2 neutralino vertex
      implicit double precision (a-h,o-z)
      double complex zn
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/neut/fnm(4),zn(4,4)
      v_nnp = ((zh(1,k)*zn(3,j)-zh(2,k)*zn(4,j))*(zn(1,i)*st-zn(2,i)*ct)
     $     + (zh(1,k)*zn(3,i)-zh(2,k)*zn(4,i))*(zn(1,j)*st-zn(2,j)*ct))
      return
      end
 
      double complex function v_ccp(i,j,k)
c     Pseudoscalar - 2 chargino vertex
      implicit double precision (a-h,o-z)
      double complex zpos,zneg
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      v_ccp = zh(1,k)*zpos(1,i)*zneg(2,j) + zh(2,k)*zpos(2,i)*zneg(1,j)
      return
      end
 
      double precision function v_pppp(i,j,k,l)
c     4 - pseudoscalar vertex
      implicit double precision (a-h,o-z)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      v_pppp = 3*zh(1,i)*zh(1,j)*zh(1,k)*zh(1,l)
     $     + 3*zh(2,i)*zh(2,j)*zh(2,k)*zh(2,l)
     $     - zh(1,i)*zh(1,j)*zh(2,k)*zh(2,l)
     $     - zh(1,i)*zh(2,j)*zh(1,k)*zh(2,l)
     $     - zh(1,i)*zh(2,j)*zh(2,k)*zh(1,l)
     $     - zh(2,i)*zh(1,j)*zh(1,k)*zh(2,l)
     $     - zh(2,i)*zh(1,j)*zh(2,k)*zh(1,l)
     $     - zh(2,i)*zh(2,j)*zh(1,k)*zh(1,l)
      return
      end
 
      double precision function v_hhpp(i,j,k,l)
c     2 charged Higgses + 2 pseudoscalars vertex
      implicit double precision (a-h,o-z)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/eps/ep(2,2)
      v_hhpp = ah(i,j)*ah(k,l)/ct2 + ep(i,k)*ep(j,l) + ep(i,l)*ep(j,k)
      return
      end
 
      double complex function v_llpp(i,j,k,l)
c     2 pseudoscalars - 2 sleptons vertex
      implicit double precision (a-h,o-z)
      double complex zv,zl
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/yukawa/yl(3),yu(3),yd(3)
      common/delta/del(6,6)
      v_llpp = e2/2/ct2*ah(i,j)*del(k,l)
      do m=1,3
         v_llpp = v_llpp + e2/4/sct2*ah(i,j)*(1-4*st2)
     $        *dconjg(zl(m,k))*zl(m,l)
         v_llpp = v_llpp - yl(m)**2*zh(1,i)*zh(1,j)
     $        *(dconjg(zl(m,k))*zl(m,l) + dconjg(zl(m+3,k))*zl(m+3,l))
      end do
      return
      end
 
      double complex function v_uupp(i,j,k,l)
c     2 pseudoscalars - 2 U-squarks vertex
      implicit double precision (a-h,o-z)
      double complex zu,zd
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/yukawa/yl(3),yu(3),yd(3)
      common/delta/del(6,6)
      v_uupp = - e2/3/ct2*ah(i,j)*del(k,l)
      do m=1,3
         v_uupp = v_uupp - e2/12/sct2*ah(i,j)*(3-8*st2)
     $        *dconjg(zu(m,k))*zu(m,l)
         v_uupp = v_uupp - yu(m)**2*zh(2,i)*zh(2,j)
     $        *(dconjg(zu(m,k))*zu(m,l) + dconjg(zu(m+3,k))*zu(m+3,l))
      end do
      return
      end
 
      double complex function v_ddpp(i,j,k,l)
c     2 pseudoscalars - 2 D-squarks vertex
      implicit double precision (a-h,o-z)
      double complex zu,zd
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/yukawa/yl(3),yu(3),yd(3)
      common/delta/del(6,6)
      v_ddpp = e2/6/ct2*ah(i,j)*del(k,l)
      do m=1,3
         v_ddpp = v_ddpp + e2/12/sct2*ah(i,j)*(3-4*st2)
     $        *dconjg(zd(m,k))*zd(m,l)
         v_ddpp = v_ddpp - yd(m)**2*zh(1,i)*zh(1,j)
     $        *(dconjg(zd(m,k))*zd(m,l)+dconjg(zd(m+3,k))*zd(m+3,l))
      end do
      return
      end

      double complex function v_udh(i,j,k)
c     Charged Higgs - u-squark - d-squark vertex
      implicit double precision (a-h,o-z)
      double complex vudh,vi_udh
      logical init_udh
      common/hud/vudh(6,6,2),init_udh
      if (init_udh) then
         do l=1,6
            do m=1,6
               do n=1,2
                  vudh(l,m,n) = vi_udh(l,m,n)
               end do
            end do
         end do
         init_udh = .false.
      end if
      v_udh = vudh(i,j,k)
      return
      end
 
      double complex function vi_udh(i,j,k)
c     Charged Higgs - u-squark - d-squark vertex initialization
      implicit double precision (a-h,o-z)
      double complex zd,zu
      double complex h,ls,ks,ds,es,us,ws,ckm
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/hpar/hm1,hm2,hs,h
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/yukawa/yl(3),yu(3),yd(3)
      common/km_mat/ckm(3,3)
      common/vev/v1,v2
      common/delta/del(6,6)
      vi_udh = (0.d0,0.d0)
      do l=1,3
         do m=1,3
            vi_udh = vi_udh + sq2*dconjg(ckm(l,m))*(( - e2/4/st2*cr(k) 
     $           + v1*yd(l)**2/2*zh(1,k) + v2*yu(m)**2/2*zh(2,k))
     $           * dconjg(zd(l,j)*zu(m,i))
     $           - wm*st/e*del(1,k)*yu(m)*yd(l)
     $           * dconjg(zd(l+3,j)*zu(m+3,i))
     $           + yu(m)/sq2*zh(1,k)*dconjg(h*zd(l,j)*zu(m+3,i))
     $           - yd(l)/sq2*zh(2,k)*h*dconjg(zd(l+3,j)*zu(m,i)))
            do n=1,3
               vi_udh = vi_udh + (ws(n,m)*zh(1,k) - us(n,m)*zh(2,k))
     $              * dconjg(zd(l,j)*zu(m+3,i)*ckm(l,n))
     $              + dconjg((ds(n,l)*zh(1,k) + es(n,l)*zh(2,k))
     $              * zd(l+3,j)*zu(m,i)*ckm(n,m))
            end do
         end do
      end do
      return
      end
 
      double complex function v_ddhh(i,j,k,l)
c     2 charged Higgs - 2 d-squark vertex
      implicit double precision (a-h,o-z)
      double complex vddhh,vi_ddhh
      integer o,p
      logical init_ddhh
      common/hhdd/vddhh(6,6,2,2),init_ddhh
      if (init_ddhh) then
         do m=1,6
            do n=1,6
               do o=1,2
                  do p=1,2
                     vddhh(m,n,o,p) = vi_ddhh(m,n,o,p)
                  end do
               end do
            end do
         end do
         init_ddhh = .false.
      end if
      v_ddhh = vddhh(i,j,k,l)
      return
      end
 
      double complex function vi_ddhh(i,j,k,l)
c     2 charged Higgs - 2 D-squark vertex initialization
      implicit double precision (a-h,o-z)
      integer o
      double complex zu,zd,ckm
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/yukawa/yl(3),yu(3),yd(3)
      common/delta/del(6,6)
      common/km_mat/ckm(3,3)
      vi_ddhh = e2/6/ct2*ah(k,l)*del(i,j)
      do m=1,3
         vi_ddhh = vi_ddhh - e2/12/sct2*ah(k,l)*(3 - 2*st2)
     $        *dconjg(zd(m,i))*zd(m,j)
     $        - yd(m)**2*zh(1,k)*zh(1,l)*dconjg(zd(m+3,i))*zd(m+3,j)
         do n=1,3
            do o=1,3  
               vi_ddhh = vi_ddhh - yu(n)**2*zh(2,k)*zh(2,l)
     $              *ckm(m,n)*zd(m,j)*dconjg(ckm(o,n)*zd(o,i))
            end do
         end do
      end do
      return
      end
 
      double complex function v_uuhh(i,j,k,l)
c     2 charged Higgs - 2 u-squark vertex
      implicit double precision (a-h,o-z)
      double complex vuuhh,vi_uuhh
      integer o,p
      logical init_uuhh
      common/hhuu/vuuhh(6,6,2,2),init_uuhh
      if (init_uuhh) then
         do m=1,6
            do n=1,6
               do o=1,2
                  do p=1,2
                     vuuhh(m,n,o,p) = vi_uuhh(m,n,o,p)
                  end do
               end do
            end do
         end do
         init_uuhh = .false.
      end if
      v_uuhh = vuuhh(i,j,k,l)
      return
      end
 
      double complex function vi_uuhh(i,j,k,l)
c     2 charged Higgs - 2 u-squark vertex initialization
      implicit double precision (a-h,o-z)
      integer o
      double complex zu,zd,ckm
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/yukawa/yl(3),yu(3),yd(3)
      common/delta/del(6,6)
      common/km_mat/ckm(3,3)
      vi_uuhh = - e2/3/ct2*ah(k,l)*del(i,j)
      do m=1,3
         vi_uuhh = vi_uuhh + e2/12/sct2*ah(k,l)*(3 + 2*st2)
     $        *dconjg(zu(m,i))*zu(m,j)
     $        - yu(m)**2*zh(2,k)*zh(2,l)*dconjg(zu(m+3,i))*zu(m+3,j)
         do n=1,3
            do o=1,3  
               vi_uuhh = vi_uuhh - yd(m)**2*zh(1,k)*zh(1,l)
     $              *ckm(m,n)*zu(n,j)*dconjg(ckm(m,o)*zu(o,i))
            end do
         end do
      end do
      return
      end
 
      double complex function v_udhs(i,j,k,l)
c     Charged Higgs - scalar - u-squark - d-squark vertex
      implicit double precision (a-h,o-z)
      double complex vudhs,vi_udhs
      integer o,p
      logical init_udhs
      common/hsud/vudhs(6,6,2,2),init_udhs
      if (init_udhs) then
         do m=1,6
            do n=1,6
               do o=1,2
                  do p=1,2
                     vudhs(m,n,o,p) = vi_udhs(m,n,o,p)
                  end do
               end do
            end do
         end do
         init_udhs = .false.
      end if
      v_udhs = vudhs(i,j,k,l)
      return
      end
 
      double complex function vi_udhs(i,j,k,l)
c     2 charged Higgs - 2 u-squark vertex initialization
      implicit double precision (a-h,o-z)
      double complex zu,zd,ckm
      double complex yl,yu,yd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/yukawa/yl(3),yu(3),yd(3)
      common/km_mat/ckm(3,3)
      vi_udhs = (0.d0,0.d0)
      do m=1,3
         do n=1,3
            vi_udhs = vi_udhs 
     $           + (zu(n,i)*zd(m,j)*(yu(n)**2*zh(2,k)*zr(2,l)
     $           + yd(m)**2*zh(1,k)*zr(1,l) 
     $           - e2/2/st2*(zh(1,k)*zr(1,l) + zh(2,k)*zr(2,l)))
     $           - yu(n)*yd(m)*ap(l,k)*zu(n+3,i)*zd(m+3,j))/sq2*ckm(m,n)
         end do
      end do
      return
      end
 
      double complex function vl_nch(i,j,k)
c     Neutralino-chargino-charged Higgs left vertex
      implicit double precision (a-h,o-z)
      double complex zpos,zneg,zn
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/neut/fnm(4),zn(4,4)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      vl_nch = (zneg(2,j)*(zn(1,i)*st + zn(2,i)*ct)/sq2
     $     - zneg(1,j)*zn(3,i)*ct)*zh(1,k)
      return
      end

      double complex function vr_nch(i,j,k)
c     Neutralino-chargino-charged Higgs right vertex
      implicit double precision (a-h,o-z)
      double complex zpos,zneg,zn
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/neut/fnm(4),zn(4,4)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      vr_nch = - dconjg(zpos(2,j)*(zn(1,i)*st + zn(2,i)*ct)/sq2
     $     + zpos(1,j)*zn(4,i)*ct)*zh(2,k)
      return
      end
 

