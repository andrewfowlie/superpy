c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: YUK_REN.F
c     Released: 15:02:2011(J.R.)
c     Changelog:
c     20:09:2012 (J.R.) corrected bug in calculation of effective
c     gluino-fermion vertices (zu_eff, zd_eff)
c     20:09:2012 (J.R.) corrected calculation of effective neutral
c     Yukawa couplings for ilev=0 (no resummation) - previously taken at
c     tree level instead of unresummed 1-loop
c     11:07:2013 (J.R.) minor additions, common/yukawa_ren/ stores both
c     tree level and bare Yukawa couplings, set_yukawa switches them

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Renormalization of the Yukawa couplings and CKM matrix          c
c     in the decoupling limit v/M_SUSY << 1                           c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine reset_tree_level
c     reset fermion and sfermion sector to tree level couplings
      implicit double precision(a-h,o-z)
      logical sqdiag,sldiag
      call init_tree_yukawa     ! reset Yukawa 
      call set_ckm(0)           ! reinitialize CKM
c     reinitialize sfermion sector, shouldn't fail at this point!
c     incorrect SUSY parameters should be eliminated by prior checks
      if (sldiag().or.sqdiag()) stop 'Check SUSY parameters!'
      return
      end

      subroutine set_ckm(is)
c     initialize CKM to bare or physical values
c        is = 0      physical
c        is = 1      bare
      implicit double precision (a-h,o-z)
      double complex ckm,ckm_herm
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      common/ckm/ckm(3,3)
      common/km_mat/ckm_herm(3,3)
      do i=1,3
         do j=1,3
            if (is.eq.0) then
               ckm(i,j) = ckm_phys(i,j)
               ckm_herm(i,j) = dconjg(ckm_phys(j,i))
            else
               ckm(i,j) = ckm0(i,j)
               ckm_herm(i,j) = dconjg(ckm0(j,i))
            end if
         end do
      end do
      return
      end

      subroutine set_yukawa(is)
c     initialize Yukawa couplings to bare or physical values
c        is = 0      physical
c        is = 1      bare
      implicit double precision (a-h,o-z)
      double complex yl,yu,yd
      double complex yl0,yu0,yd0,yl_bare,yu_bare,yd_bare
      common/yukawa/yl(3),yu(3),yd(3)      
      common/yukawa_ren/yl0(3),yu0(3),yd0(3),
     $     yl_bare(3),yu_bare(3),yd_bare(3)
      do i=1,3
         if (is.eq.0) then
            yl(i) = yl0(i)
            yu(i) = yu0(i)
            yd(i) = yd0(i)
         else
            yl(i) = yl_bare(i)
            yu(i) = yu_bare(i)
            yd(i) = yd_bare(i)
         end if
      end do
      return
      end
      
      subroutine init_tree_yukawa()
c     Initialization of Yukawa coupling
      implicit double precision (a-h,o-z)
      double complex yl,yu,yd
      double complex yl0,yu0,yd0,yl_bare,yu_bare,yd_bare
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/fmass/em(3),um(3),dm(3)
      common/yukawa/yl(3),yu(3),yd(3)
      common/yukawa_ren/yl0(3),yu0(3),yd0(3),
     $     yl_bare(3),yu_bare(3),yd_bare(3)
      common/vev/v1,v2
      do i=1,3
         yu(i) = sq2*um(i)/v2
         yd(i) = - sq2*dm(i)/v1
         yl(i) = - sq2*em(i)/v1
      end do
      do i=1,3
         yu0(i) = yu(i)
         yd0(i) = yd(i)
         yl0(i) = yl(i)
      end do
      return
      end

      subroutine renormalize_yukawa_declim
c     Yukawa renormalization, decoupling limit
      implicit double precision (a-h,o-z)
      double complex yl,yu,yd
      double complex yl0,yu0,yd0,yl_bare,yu_bare,yd_bare
      double complex eps_dd,sig_dd_ny,sig_dd_yb,sig_uu_ny
      double complex eps_ll,sig_ll_ny,sig_ll_yt
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/vev/v1,v2
      common/yukawa/yl(3),yu(3),yd(3)      
      common/yukawa_ren/yl0(3),yu0(3),yd0(3),
     $     yl_bare(3),yu_bare(3),yd_bare(3)
      common/fmass/em(3),um(3),dm(3)
c     SUSY sector initialized at this point, but with CKM_phys instead
c     of CKM_0!!
c
c     calculate effective Yukawa couplings for quarks
c     Up-quark Yukawa first, cot beta contributions not important
      do i=3,1,-1 
         yu(i) =   sq2*(um(i) - sig_uu_ny(i,i))/v2
      end do
c     Down Yukawa, starting from heavier to lighter quarks
c     Bottom Yukawa first, does not depend significantly on (uncorrected
c     yet) Y_c and Y_d
      yd(3) = - sq2*(dm(3) - sig_dd_ny(3))/(v1 + v2*eps_dd(3))
c     Strange and down Yukawa with the use of bare yd(3)
      do i=1,2
         yd(i) = - sq2*(dm(i) - sig_dd_ny(i) - sig_dd_yb(i))/(v1 
     $        + v2*eps_dd(i))
      end do
c     Lepton Yukawa, starting from heavier to lighter
c     Tau Yukawa first, does not depend significantly on (uncorrected
c     yet) Y_mu and Y_e
      yl(3) = - sq2*(em(3) - sig_ll_ny(3,3))/(v1 + v2*eps_ll(3))
c     Mu and e Yukawa with the use of bare yl(3)
      do i=1,2
         yl(i) = - sq2*(em(i) - sig_ll_ny(i,i) - sig_ll_yt(i))/(v1 
     $        + v2*eps_ll(i))
      end do
      do i=1,3
         yu_bare(i) = yu(i)
         yd_bare(i) = yd(i)
         yl_bare(i) = yl(i)
      end do
      return
      end

      subroutine renormalize_ckm_declim
c     Corrections to tree level CKM in the decoupling limit
      implicit double precision (a-h,o-z)
      double complex sig_dd_nckm
      double complex usr_sig
      double complex tt11,tt12,tt13,tt23
      double complex eps_fc,eps,den
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      common/fmass/em(3),um(3),dm(3)
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
c     note: usr_sig is Andi Crivellin's sigma_uLR^*
      tt11 = ckm_phys(1,2)*dimag(usr_sig(1,2)/um(2)
     $     - sig_dd_nckm(1,2)/dm(2))     
     $     + dimag(usr_sig(1,2)*sig_dd_nckm(1,2))/um(2)/dm(2)
      tt12 = ckm_phys(1,2) + dconjg(usr_sig(1,2))/um(2) 
     $     - sig_dd_nckm(1,2)/dm(2)
      tt23 = ckm_phys(2,3) + dconjg(usr_sig(2,3))/um(3)
     $     - sig_dd_nckm(2,3)/dm(3)
      tt13 = ckm_phys(1,3) + dconjg(usr_sig(1,3))/um(3) 
     $     + dconjg(usr_sig(1,2))/um(2)*ckm_phys(2,3)
     $     - sig_dd_nckm(1,3)/dm(3) 
     $     - sig_dd_nckm(2,3)/dm(3)*ckm_phys(1,2)
     $     + (sig_dd_nckm(1,2)/dm(2) - dconjg(usr_sig(1,2))/um(2)) 
     $     * sig_dd_nckm(2,3)/dm(3)
c     Tree level CKM matrix. Store it in /ckm_switch/
      eps = eps_fc()
      den = 1 - eps
      ckm0(1,1) = 1 - abs(tt12)**2/2 + (0.d0,1.d0)*tt11
      ckm0(1,2) = tt12
      ckm0(1,3) = tt13/den
      ckm0(2,1) = - dconjg(tt12)
      ckm0(2,2) = 1 - abs(tt12)**2/2 - (0.d0,1.d0)*tt11
      ckm0(2,3) = tt23/den
      ckm0(3,1) = dconjg((tt12*tt23 - tt13)/den)
      ckm0(3,2) = - dconjg(tt23/den)
      ckm0(3,3) = (1.d0,0.d0)
      call correct_unitarity(ckm0,3)
      return
      end

      subroutine u_rotation(sig,fm,u,iloop)
c     fermion rotation matrices at loop level 0,1,2
      implicit double precision(a-h,o-z)
      double complex sig,u(3,3)
      double precision fm(3)
      external sig
      do i=1,3
         do j=1,3
            u(i,j) = (0.d0,0.d0)
         end do
         u(i,i) = (1.d0,0.d0) 
      end do
      if (iloop.eq.0) return    ! no loop corrections
      do i=1,3
         do j=1,3
            if (i.ne.j) then
               u(i,j) = (fm(j)*dconjg(sig(i,j)) + fm(i)*sig(j,i))
     $              / (fm(j)*fm(j) - fm(i)*fm(i))
            end if
         end do
      end do
      if (iloop.eq.1) then
         call correct_unitarity(u,3)
         return                 ! 1-loop corrections only
      end if
      u(1,1) = u(1,1) - (abs(sig(1,2)/fm(2))**2 
     $     + abs(sig(1,3)/fm(3))**2)/2 
      u(2,2) = u(2,2) 
     $     - (abs(sig(1,2)/fm(2))**2 + abs(sig(2,3)/fm(3))**2)/2 
      u(3,3) = u(3,3) - (abs(sig(1,3)/fm(3))**2 
     $     + abs(sig(2,3)/fm(3))**2)/2 
      u(1,2) = u(1,2) - dconjg(sig(1,3)*sig(3,2))/fm(2)/fm(3)
      u(1,3) = u(1,3) + dconjg(sig(1,2))*sig(3,2)/fm(3)/fm(3)
      u(2,1) = u(2,1) + sig(1,3)*sig(3,2)/fm(2)/fm(3)
      u(2,3) = u(2,3) + dconjg(sig(2,1))*sig(3,1)/fm(3)/fm(3)
      u(3,1) = u(3,1) + sig(1,2)*sig(2,3)/fm(2)/fm(3)
      u(3,2) = u(3,2) - dconjg(sig(1,2))*sig(1,3)/fm(2)/fm(3)
      call correct_unitarity(u,3)
      return                    ! leading reducible 2-loop terms included
      end

      subroutine evaluate_u_rotations
      implicit double precision(a-h,o-z)
      double complex usl_sig,usr_sig,dsl_sig,dsr_sig,esl_sig,esr_sig
      double complex ckm_phys,ckm0,udl,udr,uul,uur,ull,ulr
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      common/lepton_switch/ull(3,3),ulr(3,3)
      common/fmass/em(3),um(3),dm(3)
      external usl_sig,usr_sig,dsl_sig,dsr_sig,esl_sig,esr_sig
      call u_rotation(dsl_sig,dm,udl,2)
      call u_rotation(dsr_sig,dm,udr,2)
      call u_rotation(usr_sig,um,uul,2)
      call u_rotation(usl_sig,um,uur,2)
      call u_rotation(esl_sig,em,ull,1)
      call u_rotation(esr_sig,em,ulr,1)
      return
      end
     
      subroutine declim_ren(ierr)
c     renomalization sequence for the Yukawa couplings and the CKM
c     matrix in the decoupling limit v1,v2 << M_SUSY
      implicit double precision(a-h,o-z)
      logical sqdiag,sldiag
      double complex yhl,yhr,ysu,ypu,ysd,ypd
      double complex yhlr,ysl,ypl
      logical init_yukawa_eff,init_yukawa_l
      common/yukawa_eff/yhl(3,3,2),yhr(3,3,2),ysu(3,3,2),ypu(3,3,2),
     $     ysd(3,3,2),ypd(3,3,2),init_yukawa_eff
      common/yukawa_lept/yhlr(3,3,2),ysl(3,3,2),ypl(3,3,2),init_yukawa_l
      common/debug_4q/ih,ic,in,ig
      common/declim_ih/ihr
      ierr = 0
c     reset status of effective yukawa calculation
      init_yukawa_eff = .true.
      init_yukawa_l = .true.
c     Only chirally enhanced contributions required here, no Higgs diagrams
      ihr = ih
      ih = 0
c     renormalization sequence
      call renormalize_yukawa_declim
      call renormalize_ckm_declim
      call set_ckm(1)           ! initialize ckm to bare values
      if (sldiag().or.sqdiag()) then
         call  reset_tree_level ! reset all to tree level
         ierr = 2               ! set error code
         return                 ! and exit
      end if
      call evaluate_u_rotations ! wave function rotations
c     reset Higgs diagrams status
      ih = ihr
      return
      end

      subroutine exact_iterative_renormalization(errin,errout,nmax,ierr)
c     iterative Yukawa and CKM renormalization beyond the decoupling limit
      implicit double precision (a-h,o-z)
      double complex yl,yu,yd,yu_new(3),yd_new(3),yl_new(3)
      double complex yl0,yu0,yd0,yl_bare,yu_bare,yd_bare
      double complex dsl_sig,usl_sig,esl_sig
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      double complex ckm
      logical sldiag,sqdiag
      double complex yhl,yhr,ysu,ypu,ysd,ypd
      double complex yhlr,ysl,ypl
      logical init_yukawa_eff,init_yukawa_l
      common/yukawa_eff/yhl(3,3,2),yhr(3,3,2),ysu(3,3,2),ypu(3,3,2),
     $     ysd(3,3,2),ypd(3,3,2),init_yukawa_eff
      common/yukawa_lept/yhlr(3,3,2),ysl(3,3,2),ypl(3,3,2),init_yukawa_l
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/vev/v1,v2
      common/yukawa/yl(3),yu(3),yd(3)
      common/yukawa_ren/yl0(3),yu0(3),yd0(3),
     $     yl_bare(3),yu_bare(3),yd_bare(3)
      common/fmass/em(3),um(3),dm(3)
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      common/ckm/ckm(3,3)
      ierr = 0
c     reset status of effective yukawa calculation
      init_yukawa_eff = .true.
      init_yukawa_l = .true.
      do n=1,nmax
         errout = 0.d0
c     Iterate up and down Yukawas, from 3rd to 1st generation
         do i=3,1,-1
            yu_new(i) =   sq2*(um(i) - usl_sig(i,i))/v2
            errout = max(errout,abs(yu(i)/yu_new(i) - 1))
            yu(i) = yu_new(i)
            yd_new(i) = - sq2*(dm(i) - dconjg(dsl_sig(i,i)))/v1
            errout = max(errout,abs(yd(i)/yd_new(i) - 1))
            yd(i) = yd_new(i)
            yl_new(i) = - sq2*(em(i) - dconjg(esl_sig(i,i)))/v1
            errout = max(errout,abs(yl(i)/yl_new(i) - 1))
            yl(i) = yl_new(i)
c     store "bare" Yukawas in separate common block
            yu_bare(i) = yu(i)
            yd_bare(i) = yd(i)
            yl_bare(i) = yl(i)
         end do
c     Iterate bare CKM matrix         
         call evaluate_u_rotations
c     CKM^0 = U_uL V U_dL^+
         do i=1,3
            do j=1,3
               ckm0(i,j) = (0.d0,0.d0)
               do k=1,3
                  do l=1,3
                     ckm0(i,j) = ckm0(i,j) + uul(i,k)*ckm_phys(k,l)
     $                    * dconjg(udl(j,l))
                  end do
               end do
            end do
         end do
         errout = max(errout,rel_diff_mat_norm(ckm,ckm0,3))
         call set_ckm(1)
         if (sldiag().or.sqdiag()) then
            ierr = 2
            return
         end if
         if (errout.le.errin) return
      end do
      ierr = 1
      return
      end

      subroutine set_resummation_level(ilev,ierr)
c     fixes the treatment of enhanced chiral correction resummation
c       ilev = 0     no resummation, SUSY corrections structly 1-loop
c       ilev = 1     resummation using the decoupling limit
c       ilev = 2     exact iterative solution, may not always converge
      implicit double precision(a-h,o-z)
      common/resum_level/il,nmax,errin,errout
      external resummation_data
      if (ilev.eq.0) then
         call reset_tree_level
         ierr = 0
         il = 0                 ! no resummation
      else if (ilev.eq.1) then
         call declim_ren(ierr)  ! calculate bare Yukawa and CKM 
         if (ierr.eq.0) then
            il = 1              ! renormalization successful
         else
            il = 0              ! calculations failed, negative sfermion mass^2
         end if
      else if (ilev.eq.2) then
         call exact_iterative_renormalization(errin,errout,nmax,ierr)
         if (ierr.eq.0) then
            il = 2              ! renormalization successful
         else                   
c     iterative renormalization failed, reset all to tree level
            call reset_tree_level
c     and try if at least decoupling limit approach works
            call declim_ren(ierr) ! calculate initial bare Yukawa and CKM 
            if (ierr.eq.0) then
               il = 1           ! decoupling limit successful
               ierr = -1        ! but ilev=2 not reached, negative error code
            else
               il = 0           ! calculations failed again
            end if
         end if
      else
         stop 'Incorrect resummation level requested!' ! wrong parameter
      end if
      call yukawa_eff_init
      call sq_mix_eff(il)
      return
      end

      double precision function rel_diff_mat_norm(a,b,n)
c     relative norm of difference of two complex matrices
      implicit double precision(a-h,o-z)
      double complex a(n,n),b(n,n)
      diff = 0.d0
      do i=1,n
         do j=1,n
            diff = diff + 4*abs((a(i,j) - b(i,j))/(a(i,j) + b(i,j)))**2
         end do
      end do
      rel_diff_mat_norm = sqrt(diff)
      return
      end

      subroutine yukawa_eff_init
c     initialize effective Yukawa matrices
      implicit double precision (a-h,o-z)
      double complex yhl,yhr,ysu,ypu,ysd,ypd
      double complex yhlr,ysl,ypl
      double complex ed(3,3),eu(3,3),el(3,3)
      double complex tmp1(3,3),tmp2(3,3),tmp3(3,3)
      double complex sig_dd_nhol,sig_uu_nhol,sig_ll_nhol
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      double complex ull,ulr
      logical init_yukawa_eff,init_yukawa_l
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      common/yukawa_eff/yhl(3,3,2),yhr(3,3,2),ysu(3,3,2),ypu(3,3,2),
     $     ysd(3,3,2),ypd(3,3,2),init_yukawa_eff
      common/lepton_switch/ull(3,3),ulr(3,3)
      common/yukawa_lept/yhlr(3,3,2),ysl(3,3,2),ypl(3,3,2),init_yukawa_l
      common/vev/v1,v2
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/resum_level/il,nmax,errin,errout
      external resummation_data
      do k=1,2
         do i=1,3
            do j=1,3
               if (i.eq.j) then
                  ysu(i,j,k) = - zr(2,k)/v2*um(i)
                  ypu(i,j,k) = - zh(2,k)/v2*um(i)
                  ysd(i,j,k) = - zr(1,k)/v1*dm(i)
                  ypd(i,j,k) = - zh(1,k)/v1*dm(i)
                  ysl(i,j,k) = - zr(1,k)/v1*em(i)
                  ypl(i,j,k) = - zh(1,k)/v1*em(i)
                  yhlr(i,j,k) = sq2*em(i)*zh(1,k)/v1
               else
                  ysu(i,j,k) = (0.d0,0.d0)
                  ypu(i,j,k) = (0.d0,0.d0)
                  ysd(i,j,k) = (0.d0,0.d0)
                  ypd(i,j,k) = (0.d0,0.d0)
                  ysl(i,j,k) = (0.d0,0.d0)
                  ypl(i,j,k) = (0.d0,0.d0)
                  yhlr(i,j,k) = (0.d0,0.d0)
               end if
               yhl(i,j,k) = sq2*um(j)*ckm_phys(j,i)*zh(2,k)/v2
               yhr(i,j,k) = sq2*dm(i)*ckm_phys(j,i)*zh(1,k)/v1
           end do
         end do
      end do
      do i=1,3
         do j=1,3
            ed(i,j) = sig_dd_nhol(i,j)
            eu(i,j) = sig_uu_nhol(i,j)
            el(i,j) = sig_ll_nhol(i,j)
         end do
      end do
      if (il.ne.0) then        
c     Rotate self-energy corrections to include higher order resummation
         do i=1,3
            do j=1,3
               tmp1(i,j) = (0.d0,0.d0)
               tmp2(i,j) = (0.d0,0.d0)
               tmp3(i,j) = (0.d0,0.d0)
               do l=1,3
                  do m=1,3
                     tmp1(i,j) = tmp1(i,j) 
     $                    + dconjg(udl(l,i))*ed(l,m)*udr(m,j)
                     tmp2(i,j) = tmp2(i,j) 
     $                    + dconjg(uul(l,i))*eu(l,m)*uur(m,j)
                     tmp3(i,j) = tmp3(i,j) 
     $                    + dconjg(ull(l,i))*el(l,m)*ulr(m,j)
                  end do
               end do
            end do
         end do
         do i=1,3
            do j=1,3
               ed(i,j) = tmp1(i,j)
               eu(i,j) = tmp2(i,j)
               el(i,j) = tmp3(i,j)
            end do
         end do
c     Correction to charged Higgs Yukawa are also rotated by the
c     physical CKM
         do i=1,3
            do j=1,3
               tmp1(i,j) = (0.d0,0.d0)
               tmp2(i,j) = (0.d0,0.d0)
               do n=1,3
                  tmp1(i,j) = tmp1(i,j) + dconjg(eu(n,j))*ckm_phys(n,i)
                  tmp2(i,j) = tmp2(i,j) + ed(n,i)*ckm_phys(j,n)
               end do
            end do
         end do
      end if
      do i=1,3
         do j=1,3
            do k=1,2
               ysu(i,j,k) = ysu(i,j,k) 
     $              + (zr(2,k)/v2 - zr(1,k)/v1)*eu(i,j)
               ypu(i,j,k) = ypu(i,j,k) 
     $              + (zh(2,k)/v2 + zh(1,k)/v1)*eu(i,j)
               ysd(i,j,k) = ysd(i,j,k) 
     $              - (zr(2,k)/v2 - zr(1,k)/v1)*ed(i,j)
               ypd(i,j,k) = ypd(i,j,k) 
     $              + (zh(2,k)/v2 + zh(1,k)/v1)*ed(i,j)
               ysl(i,j,k) = ysl(i,j,k) 
     $              - (zr(2,k)/v2 - zr(1,k)/v1)*el(i,j)
               ypl(i,j,k) = ypl(i,j,k) 
     $              + (zh(2,k)/v2 + zh(1,k)/v1)*el(i,j)
c     For charged Yukawa correct only physical Higgs couplings and only
c     for il>0
               if ((k.eq.1).and.(il.ne.0)) then
                  yhl(i,j,k) = yhl(i,j,k) 
     $                 - sq2*zh(2,k)*v2/v1/v1*tmp1(i,j)
                  yhr(i,j,k) = yhr(i,j,k) - sq2*zh(1,k)/v1*tmp2(i,j) 
                  yhlr(i,j,k) = yhlr(i,j,k) - sq2*zh(1,k)/v1*el(j,i) 
               end if
           end do
         end do
      end do
      init_yukawa_eff = .false.
      init_yukawa_l = .false.
      return
      end

      subroutine sq_mix_eff(il)
c     Resummed Z_U, Z_D matrices
      implicit double precision(a-h,o-z)
      double complex zu,zd
      double complex zu_eff,zd_eff
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/sqmass_eff/sum_eff(6),sdm_eff(6),zu_eff(6,6),zd_eff(6,6)
      do i=1,6
         sum_eff(i) = sum(i)
         sdm_eff(i) = sdm(i)
         do j=1,6
            zu_eff(i,j) = zu(i,j)
            zd_eff(i,j) = zd(i,j)
         end do
      end do
      if (il.eq.0) return
      do i=1,3
         do j=1,6
            zu_eff(i,j) = (0.d0,0.d0)
            zd_eff(i,j) = (0.d0,0.d0)
            do l=1,3
               zd_eff(i,j) = zd_eff(i,j) + udl(l,i)*zd(l,j)
               zd_eff(i+3,j) = zd_eff(i+3,j) + udr(l,i)*zd(l+3,j)
               zu_eff(i,j) = zu_eff(i,j) + dconjg(uul(l,i))*zu(l,j)
               zu_eff(i+3,j) = zu_eff(i+3,j)
     $              + dconjg(uur(l,i))*zu(l+3,j)
            end do
         end do
      end do
      call correct_unitarity(zd_eff,6)
      call correct_unitarity(zu_eff,6)
      return
      end

      subroutine chiral_corr_size(corr_l,corr_d,corr_u,corr_ckm)
      implicit double precision(a-h,o-z)
      double complex ckm,ckm_phys,ckm0,udl,udr,uul,uur
      double complex yl,yu,yd
      dimension corr_l(3),corr_d(3),corr_u(3),corr_ckm(3,3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/yukawa/yl(3),yu(3),yd(3)
      common/ckm/ckm(3,3)
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      common/vev/v1,v2
      common/fmass/em(3),um(3),dm(3)
      do i=1,3
         corr_l(i) = abs(1 + v1*yl(i)/sq2/em(i))
         corr_d(i) = abs(1 + v1*yd(i)/sq2/dm(i))
         corr_u(i) = abs(1 - v2*yu(i)/sq2/um(i))
         do j=1,3
            corr_ckm(i,j)=abs(1 - ckm(i,j)/ckm_phys(i,j))
         end do
      end do
      return
      end

      block data resummation_data
      implicit double precision (a-h,o-z)
      double complex yhl,yhr,ysu,ypu,ysd,ypd
      double complex yhlr,ysl,ypl
      logical init_yukawa_eff,init_yukawa_l
      common/yukawa_eff/yhl(3,3,2),yhr(3,3,2),ysu(3,3,2),ypu(3,3,2),
     $     ysd(3,3,2),ypd(3,3,2),init_yukawa_eff
      common/yukawa_lept/yhlr(3,3,2),ysl(3,3,2),ypl(3,3,2),init_yukawa_l
      common/resum_level/il,nmax,errin,errout
      data init_yukawa_eff,init_yukawa_l/2*.true./
      data il/0/                ! default resummation level (no resummation)
      data nmax/100/            ! default iterations in CKM+Y calculations 
      data errin/1.d-3/         ! default  accuracy in CKM+Y calculations 
      end

