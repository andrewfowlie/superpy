c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM} 
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor
  
c     FILENAME: VG_DEF.F
c     Last revised: 28: 9:1992
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expressions for gauge boson vertices          c
c     Compare with the paper: J.Rosiek@Phys.Rev.D41(1990)p.3464;       c
c     erratum, hep-ph/9511250                                          c
c     Common factors like e/2/sct etc. factorized from expressions.    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      double complex function v_llz(i,j)
c     Z0 - 2 slepton vertex
      implicit double precision (a-h,o-z)
      double complex zv,zl
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/delta/del(6,6)
      v_llz = - 2*st2*del(i,j)
      do k=1,3
         v_llz = v_llz + dconjg(zl(k,i))*zl(k,j)
      end do
      v_llz = e/2/sct*v_llz
      return
      end
 
      double complex function v_ddz(i,j)
c     Z0 - 2 d-squark vertex
      implicit double precision (a-h,o-z)
      double complex zd,zu
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/delta/del(6,6)
      v_ddz = - 2.d0/3*st2*del(i,j)
      do 10 k=1,3
10      v_ddz = v_ddz + dconjg(zd(k,i))*zd(k,j)
      return
      end
 
      double complex function v_uuz(i,j)
c     Z0 - 2 u-squark vertex
      implicit double precision (a-h,o-z)
      double complex zd,zu
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/delta/del(6,6)
      v_uuz = - 4.d0/3*st2*del(i,j)
      do 10 k=1,3
10      v_uuz = v_uuz + dconjg(zu(k,i))*zu(k,j)
      return
      end
 
      double complex function vl_nnz(i,j)
c     Z0 - 2 neutralino left vertex
      implicit double precision (a-h,o-z)
      double complex zn
      common/neut/fnm(4),zn(4,4)
      vl_nnz = (dconjg(zn(4,i))*zn(4,j) - dconjg(zn(3,i))*zn(3,j))/2
      return
      end
 
      double complex function vr_nnz(i,j)
c     Z0 - 2 neutralino right vertex
c     vr_nnz(i,j) = - dconjg(vl_nnz(i,j))
      implicit double precision (a-h,o-z)
      double complex zn
      common/neut/fnm(4),zn(4,4)
      vr_nnz = (zn(3,i)*dconjg(zn(3,j)) - zn(4,i)*dconjg(zn(4,j)))/2
      return
      end
 
      double complex function vl_ccz(i,j)
c     Z0 - 2 chargino left vertex
      implicit double precision (a-h,o-z)
      double complex zpos,zneg
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/delta/del(6,6)
      vl_ccz = dconjg(zpos(1,i))*zpos(1,j) + del(i,j)*(ct2 - st2)
      return
      end
 
      double complex function vr_ccz(i,j)
c     Z0 - 2 chargino right vertex
      implicit double precision (a-h,o-z)
      double complex zpos,zneg
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/delta/del(6,6)
      vr_ccz = zneg(1,i)*dconjg(zneg(1,j)) + del(i,j)*(ct2 - st2)
      return
      end
 
      double complex function v_llzz(i,j)
c     2 Z0 - 2 slepton vertex
      implicit double precision (a-h,o-z)
      double complex zv,zl
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/delta/del(6,6)
      v_llzz = st2/ct2*del(i,j)
      do 10 k=1,3
10      v_llzz = v_llzz + (1 - 4*st2)/sct2/4*dconjg(zl(k,i))*zl(k,j)
      return
      end
 
      double complex function v_ddzz(i,j)
c     2 Z0 - 2 d-squark vertex
      implicit double precision (a-h,o-z)
      double complex zd,zu
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/delta/del(6,6)
      v_ddzz = st2/ct2/9*del(i,j)
      do 10 k=1,3
10      v_ddzz = v_ddzz + (3 - 4*st2)/sct2/12*dconjg(zd(k,i))*zd(k,j)
      return
      end
 
      double complex function v_uuzz(i,j)
c     2 Z0 - 2 u-squark vertex
      implicit double precision (a-h,o-z)
      double complex zd,zu
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/delta/del(6,6)
      v_uuzz = 4.d0*st2/ct2/9*del(i,j)
      do 10 k=1,3
10      v_uuzz = v_uuzz + (3 - 8*st2)/sct2/12*dconjg(zu(k,i))*zu(k,j)
      return
      end
 
      double complex function vl_cnw(i,j)
c     W - chargino - neutralino left vertex
      implicit double precision (a-h,o-z)
      double complex zpos,zneg,zn
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/neut/fnm(4),zn(4,4)
      vl_cnw = zn(2,i)*dconjg(zpos(1,j)) - zn(4,i)*dconjg(zpos(2,j))/sq2
      return
      end
 
      double complex function vr_cnw(i,j)
c     W - chargino - neutralino right vertex
      implicit double precision (a-h,o-z)
      double complex zpos,zneg,zn
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/neut/fnm(4),zn(4,4)
      vr_cnw = dconjg(zn(2,i))*zneg(1,j) + dconjg(zn(3,i))*zneg(2,j)/sq2
      return
      end




