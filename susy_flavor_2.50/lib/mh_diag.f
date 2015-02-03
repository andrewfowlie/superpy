c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: MH_DIAG.F
c     Revised: 25: 4:1996(J.R.)
c     Initialization variables for charged Higgs-squark vertices added
c     Subroutine set_2hdm added: see comments inside.
c     Revised: 25: 3:1997(J.R.)
c     Proper handling of complex chargino and neutralino parameters

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains set of routines computing masses and mixing   c
c     angles of physical particles in the MSSM.                        c
c     Compare with the paper: J.Rosiek@Phys.Rev.D41(1990)p.3464;       c
c     erratum hep-ph/9511250                                           c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      logical function sldiag()
c     Slepton masses and mixing angles
c     sldiag=true if one or more masses negative
      implicit double precision (a-h,o-z)
      logical init_lls,init_llp
      double complex zv,zl,h
      double complex lms,rms,ums,dms,qms
      double complex sl_mat,sv_mat
      double complex vlls,vllp
      double precision mv(3,3),imv(3,3),zv1(3,3),zv2(3,3)
      double precision msl(6,6),imsl(6,6),zl1(6,6),zl2(6,6),work(12)
      double complex ls,ks,ds,es,us,ws
      double complex yl,yu,yd
      logical zzs_stat,zps_stat
      common/sl_matrix/sl_mat(6,6),sv_mat(3,3)
      common/zzs_stat/zzs_cra,zzs_crb,zzs_cr,ss,nh,zzs_stat
      common/zps_stat/zps_cr,ssa,nha,zps_stat
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/hpar/hm1,hm2,hs,h
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/vev/v1,v2
      common/yukawa/yl(3),yu(3),yd(3)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/msoft/lms(3,3),rms(3,3),ums(3,3),dms(3,3),qms(3,3)
      common/rsl/vlls(6,6,2),init_lls
      common/psl/vllp(6,6,2),init_llp
c     Scalar - slepton vertices should be recalculated after this procedure
      init_lls = .true.
      init_llp = .true.
c     Cross sections should be recalculated after this procedure
      zzs_stat = .false.
      zps_stat = .false.
c     sldiag variable initialization
      sldiag = .false.
c     Sneutrino masses and mixing angles
      vmin = v1*v1 - v2*v2
      do i=1,3
        do j=1,3
          sv_mat(i,j) = lms(i,j)
        end do
        sv_mat(i,i) = sv_mat(i,i) + e2/8/sct2*vmin
      end do
      do i=1,3
        do j=1,3
          imv(i,j) = dimag(sv_mat(i,j))
          mv(i,j) = dble(sv_mat(i,j))
        end do
      end do
      call eisch1(3,3,mv,imv,vm,zv1,zv2,ierr,work)
      do i=1,3
        if (vm(i).le.0.d0) then
          sldiag = .true.
          vm(i) = - sqrt( - vm(i))
        else
          vm(i) = sqrt(vm(i))
        end if
        do j=1,3
          zv(i,j) = dcmplx(zv1(i,j),zv2(i,j))
        end do
      end do
      call correct_unitarity(zv,3)
c     Slepton masses and mixing angles
      do i=1,3
        do j=1,3
          sl_mat(i,j)       = lms(j,i)
          sl_mat(i+3,j+3)   = rms(i,j)
          sl_mat(i,j+3)     = - (v2*ks(i,j) - v1*ls(i,j))/sq2
          if (i.eq.j) then
            sl_mat(i,j)     = sl_mat(i,j)
     $          + e2/8*(1 - 2*ct2)/sct2*vmin + abs(v1*yl(i))**2/2
            sl_mat(i+3,j+3) = sl_mat(i+3,j+3) - e2/4/ct2*vmin
     $           + abs(v1*yl(i))**2/2
            sl_mat(i,j+3)   = sl_mat(i,j+3) + v2/sq2*dconjg(h*yl(i))
          end if
          sl_mat(j+3,i)     = dconjg(sl_mat(i,j+3))
        end do
      end do

c      call cr_mat_print(sv_mat,3,6)
c      call ci_mat_print(sv_mat,3,6)
c      call cr_mat_print(sl_mat,6,6)
c      call ci_mat_print(sl_mat,6,6)
c      stop

      do i=1,6
        do j=1,6
          msl(i,j)  = dble(sl_mat(i,j))
          imsl(i,j) = dimag(sl_mat(i,j))
        end do
      end do
      call eisch1(6,6,msl,imsl,slm,zl1,zl2,ierr,work)
      do i=1,6
        if (slm(i).le.0.d0) then
          sldiag = .true.
          slm(i) = - sqrt( - slm(i))
        else
          slm(i) = sqrt(slm(i))
        end if
        do j=1,6
          zl(i,j) = dcmplx(zl1(i,j),zl2(i,j))
        end do
      end do
      call correct_unitarity(zl,6)
      return
      end

      logical function sqdiag()
c     Squark masses and mixing angles
c     sqdiag=true if one or more masses negative
      implicit double precision (a-h,o-z)
      logical init_uus,init_dds,init_uup,init_ddp
      logical init_ddhh,init_uuhh,init_udh,init_udhs
      double complex vuus,vdds,vuup,vddp
      double complex vddhh,vuuhh,vudh,vudhs
      double complex zu,zd,h,ckm
      double complex lms,rms,ums,dms,qms
      double complex ls,ks,ds,es,us,ws
      double complex sd_mat,su_mat
      double complex yl,yu,yd
      double precision msu(6,6),imsu(6,6),zu1(6,6),zu2(6,6)
      double precision msd(6,6),imsd(6,6),zd1(6,6),zd2(6,6),work(12)
      logical zzs_stat,zps_stat
      common/sq_matrix/sd_mat(6,6),su_mat(6,6)
      common/zzs_stat/zzs_cra,zzs_crb,zzs_cr,ss,nh,zzs_stat
      common/zps_stat/zps_cr,ssa,nha,zps_stat
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/hpar/hm1,hm2,hs,h
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/vev/v1,v2
      common/yukawa/yl(3),yu(3),yd(3)
      common/msoft/lms(3,3),rms(3,3),ums(3,3),dms(3,3),qms(3,3)
      common/km_mat/ckm(3,3)
      common/rsu/vuus(6,6,2),init_uus
      common/rsd/vdds(6,6,2),init_dds
      common/psu/vuup(6,6,2),init_uup
      common/psd/vddp(6,6,2),init_ddp
      common/hud/vudh(6,6,2),init_udh
      common/hhdd/vddhh(6,6,2,2),init_ddhh
      common/hhuu/vuuhh(6,6,2,2),init_uuhh
      common/hsud/vudhs(6,6,2,2),init_udhs
c     Scalar - squark vertices should be recalculated after this procedure
      init_uus = .true.
      init_dds = .true.
      init_uup = .true.
      init_ddp = .true.
      init_udh = .true.
      init_ddhh = .true.
      init_uuhh = .true.
      init_udhs = .true.
c      Cross sections should be recalculated after this procedure
      zzs_stat = .false.
      zps_stat = .false.
c      sqdiag variable initialization
      sqdiag = .false.
      vmin = v1*v1 - v2*v2
      do i=1,3
        do j=1,3
c      D-squark mass matrix initialization
          sd_mat(i,j)       = qms(j,i)
          sd_mat(i+3,j+3)   = dms(i,j)
          sd_mat(i,j+3)     = - (v2*es(i,j) - v1*ds(i,j))/sq2
          if (i.eq.j) then
            sd_mat(i,j)     = sd_mat(i,j) - e2*vmin*(1 + 2*ct2)/24/sct2
     $          + abs(yd(i)*v1)**2/2
            sd_mat(i+3,j+3) = sd_mat(i+3,j+3) - e2*vmin/12/ct2
     $          + abs(yd(i)*v1)**2/2
            sd_mat(i,j+3)   = sd_mat(i,j+3) + v2/sq2*dconjg(h*yd(i))
          end if
          sd_mat(j+3,i)     = dconjg(sd_mat(i,j+3))
c      U-squark mass matrix initialization
          su_mat(i,j) = (0.d0,0.d0)
          do k=1,3
            do l=1,3
              su_mat(i,j)   = su_mat(i,j)
     $            + dconjg(qms(k,l)*ckm(l,i))*ckm(k,j)
            end do
          end do
          su_mat(i+3,j+3)   = dconjg(ums(i,j))
          su_mat(i,j+3)     = - dconjg(v2*us(i,j) + v1*ws(i,j))/sq2
          if (i.eq.j) then
            su_mat(i,j)     = su_mat(i,j) - e2*vmin*(1 - 4*ct2)/24/sct2
     $            + abs(yu(i)*v2)**2/2
            su_mat(i+3,j+3) = su_mat(i+3,j+3) + e2/6/ct2*vmin
     $           + abs(yu(i)*v2)**2/2
            su_mat(i,j+3)   = su_mat(i,j+3) - v1/sq2*h*yu(i)
          end if
          su_mat(j+3,i)     = dconjg(su_mat(i,j+3))
        end do
      end do

c      call cr_mat_print(sd_mat,6,6)
c      call ci_mat_print(sd_mat,6,6)
c      stop

      do i=1,6
        do j=1,6
          msd(i,j)  = dble(sd_mat(i,j))
          imsd(i,j) = dimag(sd_mat(i,j))
          msu(i,j)  = dble(su_mat(i,j))
          imsu(i,j) = dimag(su_mat(i,j))
        end do
      end do
      call eisch1(6,6,msd,imsd,sdm,zd1,zd2,ierr,work)
      call eisch1(6,6,msu,imsu,sum,zu1,zu2,ierr,work)
      do i=1,6
        if (sdm(i).le.0.d0) then
          sqdiag = .true.
          sdm(i) = - sqrt( - sdm(i))
        else
          sdm(i) = sqrt(sdm(i))
        end if
        if (sum(i).le.0.d0) then
          sqdiag = .true.
          sum(i) = - sqrt( - sum(i))
        else
          sum(i) = sqrt(sum(i))
        end if
        do j=1,6
          zu(i,j) = dcmplx(zu1(i,j),zu2(i,j))
          zd(i,j) = dcmplx(zd1(i,j),zd2(i,j))
        end do
      end do
      call correct_unitarity(zd,6)
      call correct_unitarity(zu,6)
      return
      end

      logical function cdiag()
c     Chargino masses and mixing angles
c     cdiag=true if one or two masses less then 1 MeV
      implicit double precision (a-h,o-z)
      double complex zpos,zneg,x(2,2)
      double complex ztmp(2,2),sgn(2)
      double complex h,gm2,gm3
      double precision z1(2,2),z2(2,2),z3(2,2),z4(2,2)
      double precision z5(2,2),z6(2,2),z7(2,2),z8(2,2),work(4)
      double precision diag(2,2),swap(2,2)
      logical zzs_stat,zps_stat
      common/zzs_stat/zzs_cra,zzs_crb,zzs_cr,ss,nh,zzs_stat
      common/zps_stat/zps_cr,ssa,nha,zps_stat
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/vev/v1,v2
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/hpar/hm1,hm2,hs,h
      common/gmass/gm1,gm2,gm3
      common/nc_suppress/eps_d,eps_u,acc
      external init_nc_diag
      data swap/0.d0,1.d0,1.d0,0.d0/
      cdiag = .false.
c     Cross sections should be recalculated after this procedure
      zzs_stat = .false.
      zps_stat = .false.
c     Chargino mass matrix initialization
      x(1,1) = gm2
      x(1,2) = e*v2/sq2/st*eps_u
      x(2,1) = e*v1/sq2/st*eps_d
      x(2,2) = h

c      call cr_mat_print(x,2,6)
c      call ci_mat_print(x,2,6)
c      stop

c     Build X*Xherm and X*hermX
      do i=1,2
        do j=1,2
          z1(i,j) = 0.d0
          z2(i,j) = 0.d0
          z3(i,j) = 0.d0
          z4(i,j) = 0.d0
          do k=1,2
            z1(i,j) = z1(i,j) + dble(dconjg(x(k,i))*x(k,j))
            z2(i,j) = z2(i,j) + dimag(dconjg(x(k,i))*x(k,j))
            z3(i,j) = z3(i,j) + dble(dconjg(x(i,k))*x(j,k))
            z4(i,j) = z4(i,j) + dimag(dconjg(x(i,k))*x(j,k))
          end do
        end do
      end do
      call eisch1(2,2,z1,z2,fcm,z5,z6,ierr,work)
      call eisch1(2,2,z3,z4,fcm,z7,z8,ierr,work)
c      Build mixing matrices up to possible rows swap and mass sign chang
      do i=1,2
        do j=1,2
          zpos(i,j) = dcmplx(z5(i,j),z6(i,j))
          zneg(i,j) = dcmplx(z7(i,j),z8(i,j))
        end do
      end do
c      Check is row swapping is necessary. If yes, perform it.
      do i=1,2
        do j=1,2
          ztmp(i,j) = (0.d0,0.d0)
          do k=1,2
            do l=1,2
              ztmp(i,j) = ztmp(i,j) + zneg(k,i)*zpos(l,j)*x(k,l)
            end do
          end do
          diag(i,j) = abs(ztmp(i,j))
        end do
      end do
      if (max(diag(1,1),diag(2,2)).lt.max(diag(1,2),diag(2,1))) then
        do i=1,2
          do j=1,2
            ztmp(i,j) = (0.d0,0.d0)
            do k=1,2
              ztmp(i,j) = ztmp(i,j) + zpos(i,k)*swap(k,j)
            end do
          end do
        end do
        do i=1,2
          do j=1,2
            zpos(i,j) = ztmp(i,j)
          end do
        end do
      end if
c      Find complex chargino masses and rephasing vector
      do i=1,2
        sgn(i) = (0.d0,0.d0)
        do k=1,2
          do l=1,2
            sgn(i) = sgn(i) + zneg(k,i)*zpos(l,i)*x(k,l)
          end do
        end do
        if (abs(sgn(i)).eq.0) then
          sgn(i) = (1.d0,0.d0)
        else
          sgn(i) = sgn(i)/abs(sgn(i))
        end if
      end do
      do i=1,2
        do j=1,2
          zpos(i,j) = dconjg(sgn(j))*zpos(i,j)
        end do
      end do
      call correct_unitarity(zpos,2)
      call correct_unitarity(zneg,2)
c      Finally, calculate chargino masses. Check if not too small.
      do i=1,2
        fcm(i) = 0.d0
        do k=1,2
          do l=1,2
            fcm(i) = fcm(i) + dble(zneg(k,i)*zpos(l,i)*x(k,l))
          end do
        end do
        if (fcm(i).le.zm/2) cdiag = .true.
      end do
      return
      end

      logical function ndiag()
c     Neutralino masses and mixing angles
c     ndiag=true if one or more neutralinos is lighter than 1 MeV
      implicit double precision (a-h,o-z)
      double complex zn,h,y(4,4),sgn(4),gm2,gm3
      double precision z1(4,4),z2(4,4),z3(4,4),z4(4,4),work(8)
      logical zzs_stat,zps_stat
      common/zzs_stat/zzs_cra,zzs_crb,zzs_cr,ss,nh,zzs_stat
      common/zps_stat/zps_cr,ssa,nha,zps_stat
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/vev/v1,v2
      common/neut/fnm(4),zn(4,4)
      common/gmass/gm1,gm2,gm3
      common/hpar/hm1,hm2,hs,h
      common/nc_suppress/eps_d,eps_u,acc
      external init_nc_diag
c     Cross sections should be recalculated after this procedure
      zzs_stat = .false.
      zps_stat = .false.
c     ndiag variable initialization
      ndiag = .false.
c     Neutralino mass matrix initialization
      g3 = e/ct/2
      g2 = e/st/2
      y(1,1) = gm3
      y(1,2) = 0.d0
      y(1,3) = - v1*g3*eps_d
      y(1,4) = v2*g3*eps_u
      y(2,2) = gm2
      y(2,3) = v1*g2*eps_d
      y(2,4) = - v2*g2*eps_u
      y(3,3) = 0.d0
      y(3,4) = - h
      y(4,4) = 0.d0
      do i=2,4
         do j=1,i-1
            y(i,j) = y(j,i)
         end do
      end do
c     Build Y^{\dagger}*Y
      do i=1,4
         do j=1,4
            z1(i,j) = 0.d0
            z2(i,j) = 0.d0
            do k=1,4
               z1(i,j) = z1(i,j) + dble(y(k,j)*dconjg(y(k,i)))
               z2(i,j) = z2(i,j) + dimag(y(k,j)*dconjg(y(k,i)))
            end do
         end do
      end do
      call eisch1(4,4,z1,z2,fnm,z3,z4,ierr,work)
c     Build phase matrix
      do i=1,4
         do j=1,4
            zn(i,j) = dcmplx(z3(i,j),z4(i,j))
         end do
      end do
c     Find complex neutralino masses and rephasing vector
      do i=1,4
         sgn(i) = (0.d0,0.d0)
         do k=1,4
            do l=1,4
               sgn(i) = sgn(i) + zn(k,i)*y(k,l)*zn(l,i)
            end do
         end do
         if (abs(sgn(i)).eq.0) then
            sgn(i) = (1.d0,0.d0)
         else
            sgn(i) = sqrt(sgn(i)/abs(sgn(i)))
         end if
      end do
c     Find diagonalization matrix
      do i=1,4
         do j=1,4
            zn(i,j) = zn(i,j)*dconjg(sgn(j))
         end do
      end do
      call correct_unitarity(zn,4)
c     Find neutralino masses. Check if not too low
      do i=1,4
         fnm(i) = 0.d0
         do k=1,4
            do l=1,4
               fnm(i) = fnm(i) + dble(zn(k,i)*y(k,l)*zn(l,i))
            end do
         end do
         if (fnm(i).le.zm/2) ndiag = .true.
      end do
      return
      end

      logical function hdiag()
c     Higgs masses and mixing matrices
c     hdiag=true if one or more masses less equal to zero
      implicit double precision (a-h,o-z)
      logical init_lls,init_uus,init_dds
      logical init_llp,init_uup,init_ddp
      logical init_ddhh,init_uuhh,init_udh,init_udhs
      double complex h,vuus,vdds,vlls,vuup,vddp,vllp
      double complex vddhh,vuuhh,vudh,vudhs
      logical zzs_stat,zps_stat
      common/zzs_stat/zzs_cra,zzs_crb,zzs_cr,ss,nh,zzs_stat
      common/zps_stat/zps_cr,ssa,nha,zps_stat
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/hmass_mssm/cmm(2),rmm(2),pmm(2),zzr(2,2),zzh(2,2)
      common/thdm/alpha0,beta0,am0,hm10,hm20,hmc0,istat0
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/vev/v1,v2
      common/hangle/ca,sa,cb,sb
      common/hpar/hm1,hm2,hs,h
      common/rsl/vlls(6,6,2),init_lls
      common/rsu/vuus(6,6,2),init_uus
      common/rsd/vdds(6,6,2),init_dds
      common/psl/vllp(6,6,2),init_llp
      common/psu/vuup(6,6,2),init_uup
      common/psd/vddp(6,6,2),init_ddp
      common/hud/vudh(6,6,2),init_udh
      common/hhdd/vddhh(6,6,2,2),init_ddhh
      common/hhuu/vuuhh(6,6,2,2),init_uuhh
      common/hsud/vudhs(6,6,2,2),init_udhs
      external init_phys,init_const,init_control
c     Scalar - sfermion vertices should be recalculated after this procedure
      init_lls = .true.
      init_uus = .true.
      init_dds = .true.
      init_llp = .true.
      init_uup = .true.
      init_ddp = .true.
      init_udh = .true.
      init_ddhh = .true.
      init_uuhh = .true.
      init_udhs = .true.
c     Cross sections should be recalculated after this procedure
      zzs_stat = .false.
      zps_stat = .false.
c     hdiag variable initialization
      hdiag=.false.
c     Build pseudoscalar and charged Higgs mixing matrix
      zh(1,1) = v2
      zh(2,2) = v2
      zh(1,2) = -v1
      zh(2,1) = v1
      do i=1,2
        do j=1,2
          zh(i,j)  = zh(i,j)/sqrt(v1*v1 + v2*v2)
        end do
      end do
c     Store sin(beta), cos(beta) in common/hangle/
      sb = zh(1,1)
      cb = zh(2,1)
c     Calculate pseudoscalar and charged Higgs masses
      pm(1) = hm1 + hm2 + 2*abs(h*h)
      cm(1) = pm(1) + wm*wm
c     Goldstone masses (already not squared)
      pm(2) = zm
      cm(2) = wm
c     Calculate pseudoscalar and charged Higgs masses
c     Check if all masses positive
      if (pm(1).le.0.d0) then
        hdiag = .true.
        pm(1) = - sqrt( - pm(1))
c     Next if nested because if cm(1)<0 then always pm(1)<0 and hdiag=true
        if (cm(1).le.0.d0) then
          cm(1) = - sqrt( - cm(1))
        else
          cm(1) = sqrt(cm(1))
        end if
      else
        pm(1) = sqrt(pm(1))
        cm(1) = sqrt(cm(1))
      end if
c     Build scalar mixing matrix
      a = hm1 + abs(h*h) + e2/8/sct2*(3*v1*v1 - v2*v2)
      b = hm2 + abs(h*h) + e2/8/sct2*(3*v2*v2 - v1*v1)
      c = hs - e2/4/sct2*v1*v2
      x = atan(2*c/(a - b))/2
      if (x.le.0.d0) then
        ca = cos(x)
        sa = sin(x)
      else
        ca = sin(x)
        sa = - cos(x)
      end if
      zr(1,1) = ca
      zr(2,2) = ca
      zr(1,2) = - sa
      zr(2,1) = sa
c     Calculate squared scalar masses.
      rm(1) = a*ca*ca + b*sa*sa + 2*c*sa*ca
      rm(2) = a*sa*sa + b*ca*ca - 2*c*sa*ca
c     Check if scalar masses positive
      do i=1,2
         if (rm(i).le.0) then
            hdiag = .true.
            rm(i) = - sqrt( - rm(i))
         else
            rm(i) = sqrt(rm(i))
         end if
      end do
c     Initialize common/hmass_mssm/
      istat0 = 0
      do i=1,2
         cmm(i) = cm(i)
         pmm(i) = pm(i)
         rmm(i) = rm(i)
         do j=1,2
            zzr(i,j) = zr(i,j)
            zzh(i,j) = zh(i,j)
         end do
      end do
      return
      end

      subroutine correct_unitarity(z,n)
c     improves numerical accuracy of complex matrix unitarity
      implicit double precision (a-h,o-z)
      parameter (nmax=6)
      double complex z(n,n),eps(nmax,nmax),v(nmax,nmax)
c     check array sizes
      if (n.gt.nmax) stop
     $     'increase array sizes in routine correct_unitarity'
c     build matrix eps = 1 - (Z^+ Z - 1)/2 = 3/2 - (Z^+ Z)/2
c     also, store z in auxiliary matrix v
      do i=1,n
         do j=1,n
            eps(i,j) = (0.d0,0.d0)
            v(i,j) = z(i,j)
            do k=1,n
               eps(i,j) = eps(i,j) - dconjg(z(k,i))*z(k,j)/2
            end do
         end do
         eps(i,i) = 1.5d0 + eps(i,i)
      end do
c     correct initial matrix Z -> Z*eps
      do i=1,n
         do j=1,n
            z(i,j) = (0.d0,0.d0)
            do k=1,n
               z(i,j) = z(i,j) + v(i,k)*eps(k,j)
            end do
         end do
      end do
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     The following set of routines computes approximate 2-loop masses c
c     and mixing angles of Higgs particles in the MSSM                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function alpha_eff(nh)
      implicit double precision (a-h,o-z)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/hmass_EPA/pm,hm1,hm2,sa,ca,sb,cb
c     Effective alpha (argument nh unused).
      alpha_eff = atan(sa/ca)
      if (alpha_eff.gt.0) alpha_eff = alpha_eff - pi
      return
      end

      subroutine mhcorr_app2(ierr)
c     approximate version of 2-loop corrections to mh/mH based on 
c     hep-ph/9903404
      implicit double precision (a-h,o-z)
      double complex sd_mat,su_mat,mu
      common/sq_matrix/sd_mat(6,6),su_mat(6,6)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/fmass/em(3),um(3),dm(3)
      common/qmass_pole/ump(3),dmp(3)     
      common/hmass/cm(2),rm(2),ppm(2),zr(2,2),zh(2,2)
      common/hmass_EPA/pm,hm1,hm2,sa,ca,sb,cb
      common/hpar/hmpar1,hmpar2,hmpar12,mu
      common/fermi/g_fermi
      pm = ppm(1)               ! M_A
      sb = zh(1,1)              ! sin(beta)
      cb = zh(2,1)              ! cos(beta)
      tm = ump(3)**2            ! pole mt^2
      bm = dmp(3)**2            ! pole mb^2
c     stop mass eigenstates^2 (no flavor mixing included!)
      stm1 = dble(su_mat(3,3) + su_mat(6,6) 
     $     - sqrt((su_mat(3,3) - su_mat(6,6))**2 
     $     + 4*abs(su_mat(3,6))**2))/2
      stm2 = dble(su_mat(3,3) + su_mat(6,6) 
     $     + sqrt((su_mat(3,3) - su_mat(6,6))**2 
     $     + 4*abs(su_mat(3,6))**2))/2
c     sbottom mass eigenstates^2 (no flavor mixing included!)
      sbm1 = dble(sd_mat(3,3) + sd_mat(6,6) 
     $     - sqrt((sd_mat(3,3) - sd_mat(6,6))**2 
     $     + 4*abs(sd_mat(3,6))**2))/2
      sbm2 = dble(sd_mat(3,3) + sd_mat(6,6) 
     $     + sqrt((sd_mat(3,3) - sd_mat(6,6))**2 
     $     + 4*abs(sd_mat(3,6))**2))/2

      if ((stm1.le.0).or.(stm2.le.0).or.(sbm1.le.0).or.(sbm2.le.0)) then
c     negative eigemass^2
        ierr = 2
        return
      end if

c     averaged stop and sbottom mass^2
      stm =  sqrt(stm1*stm2 + tm*(stm1 + stm2) + tm*tm)
      sbm =  sqrt(sbm1*sbm2 + bm*(sbm1 + sbm2) + bm*bm)
c     "M^t_LR^2/M_S^2" and "M^b_LR^2/M_S^2" of hep-ph/9903404
      dsm = abs(su_mat(3,6)/um(3))**2/stm
      dbm = abs(sd_mat(3,6)/dm(3))**2/sbm
c     on-shell s_W^2 
      sw2 = 1 - wm2/zm2
c     auxiliary variables
      pt = 1.d0 - 8*sw2/3 + 32.d0/9*sw2*sw2
      xt = zm2/tm               ! MZ^2/mt^2
      xb = bm/zm2               ! mb^2/MZ^2
      xs = tm/stm               ! mt^2/M_stop^2
      sl = log(xs)
      gf = g_fermi*sq2/pi/pi*zm2*zm2

c     corrections from top/stop sector to CP-even Higgs self energy,
c     factor g_fermi*sq2/pi/pi*zm**4 extracted
c     1-loop
      se11 = cb*cb*pt*sl/8.d0
      se12 = - cb/sb*(pt*sb*sb - 3/xt)*sl/8.d0
      se22 = (- 2*xt + 1.1d0*xt*xt + (12 - 6*xt*sb*sb + pt*xt*xt*sb**4)
     $     *sl + dsm*(- 12 + 4*xt + 6*xs) + dsm**2*(1 - 4*xs + 3*xs*xs)
     $     + dsm**3*xs*(3/5.d0 - 12/5.d0*xs + 2*xs*xs) + dsm**4*xs*xs
     $     *(3/7.d0 - 12/7.d0*xs + 1.5d0*xs*xs))/xt/xt/sb/sb/8
c     2-loop
      se22 = se22 + 3*alfas(sqrt(stm))/pi*(sl*sl - 2*sl - 2*sqrt(dsm) -
     $     dsm*sl + dsm*dsm/4)/xt/xt/sb/sb

c     corrections from other SUSY sectors (hep-ph/9307201)
c     no treshold (non-logarithmic) corrections from Ab!
      pt = 3*pt
      pf = 21 - 40*sw2 + 160/3.d0*sw2*sw2
      pg = - 44 + 106*sw2 - 62*sw2*sw2
      pg1 = 10 + 34*sw2 - 26*sw2*sw2
      ph = - 10 + 2*sw2 - 2*sw2*sw2
      ph1 = 8 - 22*sw2 + 10*sw2*sw2

      se11 = se11 + cb*cb*(36*xb*xb/cb**4 - 18*xb/cb/cb + pf + pg + ph)
     $     *log(sbm/zm2)/24.d0
      se12 = se12 + sb*cb*(9*xb/cb/cb - pf - pg1 - ph1)*log(sbm/zm2)
     $     /24.d0
      se22 = se22 + sb*sb*(pf + pg + ph)*log(sbm/zm2)/24.d0

c     CP-even Higgs mass matrix
      a =   zm2*cb*cb + pm*pm*sb*sb - gf*se11 
      b =   zm2*sb*sb + pm*pm*cb*cb - gf*se22 
      c = - (zm2 + pm*pm)*sb*cb - gf*se12

      x = atan(2*c/(a - b))/2
      if (x.gt.0.d0) then
        ix = int(2*x/pi)
        x = x - (ix + 1)*pi/2.d0
      end if

      nrot = 0
 10   sa = sin(x)
      ca = cos(x)
      hm1 = a*ca*ca + b*sa*sa + 2*c*ca*sa
      hm2 = a*sa*sa + b*ca*ca - 2*c*ca*sa
      hm1 = sign(sqrt(abs(hm1)),hm1)
      hm2 = sign(sqrt(abs(hm2)),hm2)
      if (hm1.lt.hm2) then
        if (nrot.eq.0) then
          nrot = 1
          x = x - pi/2
          goto 10
        else
          ierr = 3
          return
        end if
      end if

c     subleading (non-logarithmic sbottom sector and 2-loop Y_t)
c     corrections to mh (hep-ph/990340 eq.25) in large M_A limit
c     trilinear sbottom corrections
      dmhb = 1.5d0*gf*xb*xb*dbm*(1 - dbm/12) - gf*(cb*cb - sb*sb)*(3*xb
     $     *dbm + abs(sd_mat(3,6)*(sd_mat(3,6) + 2*dconjg(mu)*dm(3)*sb
     $     /cb))/sbm/zm2)/8
c     2-loop Yt corrections
      if ((abs(stm1-stm2)/(stm1+stm2)).le.1.d-4) then
         dsf = 0.d0
      else
         dsf = 2 - (stm2 + stm1)/(stm2 - stm1)*log(stm2/stm1)
      end if
      sin2t = abs(su_mat(3,6)**2/(su_mat(3,6)**2 + (su_mat(3,3) -
     $     su_mat(6,6))**2/4))  ! sin^2(2 theta_t)
      dst = (stm2 - stm1)/4/um(3)**2*sin2t
      xxt = dst*dst*dsf + 2*dst*log(stm2/stm1)
      tl = log(stm1*stm2/um(3)**4)/2
      dmhy = (3*g_fermi*um(3)**3/4.d0/pi/pi)**2*tl*(xxt + tl)
c     corrected physical m_h
      hm2sq = hm2*hm2 + dmhb + dmhy
      if (hm2sq.le.0.d0) then
         ierr = 4
         return
      end if
      hm2 = sqrt(hm2sq)
      return
      end

      subroutine set_2hdm(alpha,beta,am,hm1,hm2,hmc,istat)
c     Set/reset Higgs sector parameters for 2HDM/MSSM
c     alpha, beta:   mixing angles
c     am:            pseudoscalar mass
c     hm1,hm2:       H/h masses
c     hmc:           H^+ mass
c     istat defines the required action (istat value is stored in
c     common/2hdm_args/):
c     istat=1    Restores tree level MSSM mixing angles and masses
c                Subroutine mh_diag has to be called first!
c     istat=2    MSSM alpha_eff and 1-loop mh/mH used in EPA approximation
c                Subroutines corr_EPA or fcorr_EPA has to be called first!
c     istat=3    Sets completely free values of masses and mixing angles.
c                Values of parameters alpha,...,hmc important only in this case
      implicit double precision (a-h,o-z)
      logical init_lls,init_uus,init_dds
      logical init_llp,init_uup,init_ddp
      logical init_ddhh,init_uuhh,init_udh,init_udhs
      logical zzs_stat,zps_stat
      double complex vddhh,vuuhh,vudh,vudhs
      double complex vuus,vdds,vlls,vuup,vddp,vllp
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alp,wm,wm2,zm,zm2,pi,sq2
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/vev/v1,v2
      common/hangle/ca,sa,cb,sb
      common/hmass_mssm/cmm(2),rmm(2),pmm(2),zzr(2,2),zzh(2,2)
      common/thdm/alpha0,beta0,am0,hm10,hm20,hmc0,istat0
      common/hmass_EPA/pme,hm1e,hm2e,sae,cae,sbe,cbe
      common/zzs_stat/zzs_cra,zzs_crb,zzs_cr,ss,nh,zzs_stat
      common/zps_stat/zps_cr,ssa,nha,zps_stat
      common/rsl/vlls(6,6,2),init_lls
      common/rsu/vuus(6,6,2),init_uus
      common/rsd/vdds(6,6,2),init_dds
      common/psl/vllp(6,6,2),init_llp
      common/psu/vuup(6,6,2),init_uup
      common/psd/vddp(6,6,2),init_ddp
      common/hud/vudh(6,6,2),init_udh
      common/hhdd/vddhh(6,6,2,2),init_ddhh
      common/hhuu/vuuhh(6,6,2,2),init_uuhh
      common/hsud/vudhs(6,6,2,2),init_udhs
c     Scalar - sfermion vertices should be recalculated after this procedure
      init_lls = .true.
      init_uus = .true.
      init_dds = .true.
      init_llp = .true.
      init_uup = .true.
      init_ddp = .true.
      init_udh = .true.
      init_ddhh = .true.
      init_uuhh = .true.
      init_udhs = .true.
c     Store parameter values
      alpha0 = alpha
      beta0  = beta
      am0    = am
      hm10   = hm1
      hm20   = hm2
      hmc0   = hmc
      istat0 = istat
c     Cross sections should be recalculated after this procedure
      zzs_stat = .false.
      zps_stat = .false.
      if (istat.eq.1) then
         do i=1,2
            cm(i) = cmm(i)
            pm(i) = pmm(i)
            rm(i) = rmm(i)
            do j=1,2
               zr(i,j) = zzr(i,j)
               zh(i,j) = zzh(i,j)
            end do
         end do
         sa = zr(2,1)
         ca = zr(1,1)
         sb = zh(1,1)
         cb = zh(2,1)
      else if (istat.eq.2) then
         rm(1) = hm1e
         rm(2) = hm2e
         zr(1,1) = cae
         zr(1,2) = - sae
         zr(2,1) = sae
         zr(2,2) = cae
      else if (istat.eq.3) then
         rm(1) = hm1
         rm(2) = hm2
         cm(1) = hmc
         pm(1) = am
         sa = sin(alpha)
         ca = cos(alpha)
         sb = sin(beta)
         cb = cos(beta)
         v1 = 2*wm*st/e*cb
         v2 = 2*wm*st/e*sb
         zr(1,1) = ca
         zr(1,2) = - sa
         zr(2,1) = sa
         zr(2,2) = ca
         zh(1,1) = sb
         zh(1,2) = - cb
         zh(2,1) = cb
         zh(2,2) = sb
      else
         stop 'istat out of range in set_2hdm'
      end if
      call init_tree_yukawa
      return
      end

      block data init_nc_diag
      implicit double precision (a-h,o-z)
      common/nc_suppress/eps_d,eps_u,acc
      data eps_d,eps_u/2*1.d0/
      data acc/1.d-7/
      end

