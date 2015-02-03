c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: SFLAV_IO.F
c     Released: 20:02:2010(J.R.)
c     Changelog:
c     20:09:2012 (J.R.) corrected lacking hermitization of sfermion
c     mass matrices for input_type = 2
c     12:07:2013(J.R.) generalized and relaxed input format, soft marix
c     blocks are now optional, sfermion masses can be given in BLOCK
c     EXTPAR.
c     Also, hadronic parameters can be now initialized in optional BLOCK
c     SFLAV_HADRON.

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     SLHA2 based input routine for SUSY_FLAVOR                       c
c     Test output routines for MSSM parameters and SUSY spectrum      c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      logical function find_block(ifl,name)
      implicit double precision (a-h,o-z)
      character*18 name,line
      len = len_trim(name)
      find_block=.false.
      rewind(ifl)
 10   read(ifl,'(a)',END=100) line
      if (line(7:6+len).ne.name(1:len)) goto 10
      find_block = .true.
 100  return
      end

      logical function read_var(ifl,vec,n)
c     read vec(i) from the ifl file, skipping comments
c     return true if next BLOCK statement or EOF reached
      implicit double precision (a-h,o-z)
      dimension vec(n)
      character*1 line
      read_var = .true.
 10   read(ifl,'(a)',END=100) line
      if (line(1:1).eq.'#') goto 10
      if ((line(1:1).eq.'B').or.(line(1:1).eq.'b')) return
      backspace(ifl)
      read(ifl,*)i,tmp
      if (i.ne.0) vec(i) = tmp
      read_var = .false.
 100  return
      end

      logical function read_mat_el(ifl,arr,icmpl)
c     read real (icmpl=0) or imaginary (icmpl=1) part if vec(i,j) from
c     the ifl file, skipping comments return true if next BLOCK
c     statement or EOF reached
      implicit double precision (a-h,o-z)
      double complex arr(3,3)
      character*1 line
      read_mat_el = .true.
 10   read(ifl,'(a)',END=100) line
      if (line(1:1).eq.'#') goto 10
      if ((line(1:1).eq.'B').or.(line(1:1).eq.'b')) return
      backspace(ifl)
      read(ifl,*)i,j,tmp
      if (icmpl.eq.0) then
         arr(i,j) = dcmplx(tmp,dimag(arr(i,j)))
      else
         arr(i,j) = dcmplx(dble(arr(i,j)),tmp)
      end if
      read_mat_el = .false.
 100  return
      end

      subroutine sflav_defaults
c     Initialize "default" set of parameters
      implicit double precision (a-h,o-z)
      double complex ls,ks,ds,es,us,ws
      double complex lms,rms,ums,dms,qms
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/msoft/lms(3,3),rms(3,3),ums(3,3),dms(3,3),qms(3,3)
      common/sflav_data/control_data(3),sm_data(31),ckm_data(4),
     $     extpar_data(49),extim_data(49),hadron_data(65),tb_min(6)
c     defaults for control variables
      control_data(1) = 1       ! sfermion input parameters in SLHA2 conventions
      control_data(2) = 2       ! fermion soft terms given as absolute values
      control_data(3) = 2       ! iterative numerical resummation of chiral corrections
c     SM default data
      sm_data(1) = 127.934d0    ! 1/alpha_em(M_Z)
      sm_data(3) = 0.1172d0     ! alpha_s(M_Z)
      sm_data(4) = 91.1876d0    ! M_Z (pole)
c     Fermion mass initialization, input: MSbar running quark masses
      sm_data(5) = 4.18d0       ! mb(mb) SM MSbar
      sm_data(6) = 173.5d0      ! mtop(pole)
      sm_data(7) = 1.77684d0    ! mtau(pole)
      sm_data(11) = 5.109989d-4 ! me(pole)
      sm_data(13) = 0.105658d0  ! mmu(pole)
      sm_data(21) = 4.7d-03     ! md(2 GeV) MSbar
      sm_data(22) = 2.1d-03     ! mu(2 GeV) MSbar
      sm_data(23) = 9.34d-02    ! ms(2 GeV) MSbar
      sm_data(24) = 1.279d0     ! mc(mc) MSbar
      sm_data(30) = 80.398d0    ! MW (pole), not a standard SLHA2 entry !!!
      sm_data(31) = 0.23116d0   ! s_W^2 (MSbar), not a standard SLHA2 entry !!!
c     CKM data
      ckm_data(1) = 0.2258d0    ! lambda
      ckm_data(2) = 0.808d0     ! A
      ckm_data(3) = 0.177d0     ! rho bar
      ckm_data(4) = 0.36d0      ! eta bar
c     clear MINPAR/EXTPAR and soft terms data
      do i=1,6
         tb_min(i) = 0.d0
      end do
      do i=1,49
         extpar_data(i) = 0.d0
         extim_data(i)  = 0.d0
      end do
      do i=1,3
         do j=1,3
            lms(i,j) = (0.d0,0.d0)
            rms(i,j) = (0.d0,0.d0)
            qms(i,j) = (0.d0,0.d0)
            dms(i,j) = (0.d0,0.d0)
            ums(i,j) = (0.d0,0.d0)
            ls(i,j) = (0.d0,0.d0)
            ds(i,j) = (0.d0,0.d0)
            us(i,j) = (0.d0,0.d0)
            ks(i,j) = (0.d0,0.d0) ! non-holomorphic A terms set to zero
            es(i,j) = (0.d0,0.d0) ! non-holomorphic A terms set to zero
            ws(i,j) = (0.d0,0.d0) ! non-holomorphic A terms set to zero
         end do
      end do
      hadron_data(1) = 0.1561d0 ! f_K
      hadron_data(2) = 0.2d0    ! f_D
      hadron_data(3) = 0.193d0  ! f_B_d
      hadron_data(4) = 0.232d0  ! f_B_s
      hadron_data(5) = 0.724d0  ! B_K for SM contribution to KKbar
      hadron_data(6) = 1.87d0   ! eta_cc in KK mixing (SM)
      hadron_data(7) = 0.496d0  ! eta_ct in KK mixing (SM)
      hadron_data(8) = 0.577d0  ! eta_tt in KK mixing (SM)
      hadron_data(9) = 2.d0     ! scale for B_K (non-SM)    
      hadron_data(10) = 0.61d0  ! B_K for VLL (non-SM)
      hadron_data(11) = 0.76d0  ! B_K for SLL1
      hadron_data(12) = 0.51d0  ! B_K for SLL2
      hadron_data(13) = 0.96d0  ! B_K for LR1 
      hadron_data(14) = 1.30d0  ! B_K for LR2 
      hadron_data(15) = 1.d0    ! B_D for SM contribution 
      hadron_data(16) = 2.d0    ! scale for B_D (non-SM)
      do i=1,5
         hadron_data(16+i) = 1.d0 ! B_D for VLL, SLL(2), LR(2) 
      end do 
      hadron_data(22) = 1.22d0 ! B_Bd for SM contribution 
      hadron_data(23) = 4.6d0   ! scale for B_B (non-SM, both Bd and Bs)
      hadron_data(24) = 0.87d0  ! B_Bd for VLL (non-SM)
      hadron_data(25) = 0.8d0   ! B_Bd for SLL1
      hadron_data(26) = 0.71d0  ! B_Bd for SLL2
      hadron_data(27) = 1.71d0  ! B_Bd for LR1 
      hadron_data(28) = 1.16d0  ! B_Bd for LR2 
      hadron_data(29) = 1.22d0  ! B_Bs for SM contribution 
      hadron_data(30) = 0.55d0  ! eta_b for BsBs (SM)
      hadron_data(31) = 0.87d0  ! B_Bs for VLL (non-SM)
      hadron_data(32) = 0.8d0   ! B_Bs for SLL1
      hadron_data(33) = 0.71d0  ! B_Bs for SLL2
      hadron_data(34) = 1.71d0  ! B_Bs for LR1 
      hadron_data(35) = 1.16d0  ! B_Bs for LR2 
      hadron_data(36) = 1.519d-12 ! Bd lifetime (experimental)
      hadron_data(37) = 1.512d-12 ! Bs lifetime (experimental)
      hadron_data(38) = 5.2792d0 ! Bd mass (experimental)
      hadron_data(39) = 5.3668d0 ! Bs mass (experimental)
      hadron_data(40) = 3.337d-13 ! Delta Bd (experimental)
      hadron_data(41) = 1.17d-11 ! Delta Bs (experimental)
      hadron_data(42) = 0.497614d0 ! K0 mass (experimental)
      hadron_data(43) = 3.483d-15 ! Delta mK (experimental)
      hadron_data(44) = 2.228d-3 ! eps_K (experimental)
      hadron_data(45) = 1.8645d0 ! D0 mass (experimental)
      hadron_data(46) = 1.56d-14 ! Delta mD (experimental)
      hadron_data(47) = 2.231d-10 ! parameter kappa0 in K^0->pi^0vv calculations
      hadron_data(48) = 5.173d-11 ! parameter kappa+ in K^+->pi^+vv calculations
      hadron_data(49) = 0.41d0  ! parameter P_c in K->pivv calculations
      hadron_data(50) = 0.013d-10 ! error of kappa0
      hadron_data(51) = 0.024d-11 ! error of kappa+
      hadron_data(52) = 0.03d0  ! error of P_c 
      hadron_data(53) = 0.79d0  ! neutron EDM_d QCD coefficient
      hadron_data(54) = -.2d0   ! neutron EDM_u QCD coefficient
      hadron_data(55) = 0.59d0  ! neutron CDM_d QCD coefficient
      hadron_data(56) = 0.3d0   ! neutron CDM_u QCD coefficient
      hadron_data(57) = 3.4d0   ! neutron CDM_g QCD coefficient
      hadron_data(58) = 1.18d0  ! neutron EDM chiral symmetry breaking scale
      hadron_data(59) = 1.5d0   ! pole c quark mass (in B-->X_s gamma and t->cH)
      hadron_data(60) = 0.1872d0 ! Br(tau->evv)
      hadron_data(61) = 5.27917d0 ! M_B_u
      hadron_data(62) = 0.297d0 ! Br(B->D tau nu)/Br(B->D l nu) in SM
      hadron_data(63) = 0.017d0 ! error of Br(B->D tau nu)/Br(B->D l nu) in SM
      hadron_data(64) = 0.252d0 ! Br(B->D* tau nu)/Br(B->D* l nu) in SM
      hadron_data(65) = 0.003d0 ! error of Br(B->D* tau nu)/Br(B->D* l nu) in SM
      return
      end

      subroutine sflav_sm
c     set default SM/hadronic parameters
      implicit double precision (a-h,o-z)
      common/edm_qcd/eta_ed,eta_eu,eta_cd,eta_cu,eta_g,alamx
      common/kpivv/ak0,del_ak0,akp,del_akp,pc,del_pc,alam_k
      common/sm_4q/eta_cc,eta_ct,eta_tt,eta_b,bk_sm,bd_sm,bb_sm(2)
      common/meson_data/dmk,amk,epsk,fk,dmd,amd,fd,
     $    amb(2),dmb(2),tau_b(2),fb(2)
      common/bx_4q/bk(5),bd(5),bb(2,5),amu_k,amu_d,amu_b
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
      common/fmass/em(3),um(3),dm(3)
      common/qmass_pole/ump(3),dmp(3)
      common/tau_gam/br_tau_evv
      common/dtau_data/dmbp,rd,del_rd,rds,del_rds
      common/sflav_data/control_data(3),sm_data(31),ckm_data(4),
     $     extpar_data(49),extim_data(49),hadron_data(65),tb_min(6)
      external init_phys,init_4q,init_const,init_control,init_units
c     set default SM variables, clear SUSY parameters
      call sflav_defaults
      alpha_z = 1/sm_data(1)    ! 1/alpha_em(M_Z)
      zm0 = sm_data(4)          ! M_Z
      wm0 = sm_data(30)         ! M_W
      st2_new = sm_data(31)     ! s_W^2 (here MSbar value)
      call vpar_update(zm0,wm0,alpha_z,st2_new)
c     QCD parameters and fermion mass initialization
c     input: MSbar running quark masses
      alpha_s = sm_data(3)      ! alpha_s(MZ)
      em(3) = sm_data(7)        ! m_tau (pole)
      em(1) = sm_data(11)       ! m_e (pole)
      em(2) = sm_data(13)       ! m_mu (pole)
      dml(1) = sm_data(21)      ! md(2 GeV)
      uml(1) = sm_data(22)      ! mu(2 GeV)
      dml(2) = sm_data(23)      ! ms(2 GeV)
      uml(2) = sm_data(24)      ! mc(mc)
      amuu(2) = uml(2)          ! charm mass assumed to be given as mc(mc) 
      bottom = sm_data(5)       ! mb(mb)
      bot_scale = bottom
      top_pole = sm_data(6)     ! mt(pole)
      top = 0.94d0*top_pole     ! mt(mt_pole) - running mt at mt(pole)
      top_scale = top
      call init_fermion_sector(alpha_s,top,top_scale,bottom,bot_scale)
c     read the CKM data
      alam = ckm_data(1)        ! lambda
      apar = ckm_data(2)        ! A
      rhobar = ckm_data(3)      ! rho bar
      etabar = ckm_data(4)      ! eta bar
      call ckm_wolf(alam,apar,rhobar,etabar)
c     read the hadronic data
      fk = hadron_data(1)       ! f_K
      fd = hadron_data(2)       ! f_D
      fb(1) = hadron_data(3)    ! f_B_d
      fb(2) = hadron_data(4)    ! f_B_s
      bk_sm = hadron_data(5)    ! B_K for SM contribution to KKbar
      eta_cc = hadron_data(6)   ! eta_cc in KK mixing (SM)
      eta_ct = hadron_data(7)   ! eta_cc in KK mixing (SM)
      eta_tt = hadron_data(8)   ! eta_cc in KK mixing (SM)
      amu_k = hadron_data(9)    ! scale for B_K (non-SM)    
      do i=1,5
         bk(i) = hadron_data(9+i) ! B_K for VLL, SLL(2), LR(2) 
      end do 
      bd_sm = hadron_data(15)   ! B_D for SM contribution 
      amu_d = hadron_data(16)   ! scale for B_D (non-SM)
      do i=1,5
         bd(i) = hadron_data(16+i) ! B_D for VLL, SLL(2), LR(2) 
      end do 
      bb_sm(1) = hadron_data(22) ! B_Bd for SM contribution 
      amu_b = hadron_data(23)   ! scale for B_B (non-SM, both Bd and Bs)
      do i=1,5
         bb(1,i) = hadron_data(23+i) ! B_Bd for VLL, SLL(2), LR(2) 
      end do 
      bb_sm(2) = hadron_data(29) ! B_Bs for SM contribution 
      eta_b = hadron_data(30)   ! eta_b for  BsBsbar (SM)
      do i=1,5
         bb(2,i) = hadron_data(30+i) ! B_Bs for VLL, SLL(2), LR(2) 
      end do 
      tau_b(1) = hadron_data(36) ! Bd lifetime (experimental)
      tau_b(2) = hadron_data(37) ! Bs lifetime (experimental)
      amb(1) = hadron_data(38)  ! Bd mass (experimental)
      amb(2) = hadron_data(39)  ! Bs mass (experimental)
      dmb(1) = hadron_data(40)  ! Delta Bd (experimental)
      dmb(2) = hadron_data(41)  ! Delta Bs (experimental)
      amk = hadron_data(42)     ! K0 mass (experimental)
      dmk = hadron_data(43)     ! Delta mK (experimental)
      epsk = hadron_data(44)    ! eps_K (experimental)
      amd = hadron_data(45)     ! D0 mass (experimental)
      dmd = hadron_data(46)     ! Delta mD (experimental)
      ak0 = hadron_data(47)     ! parameter kappa0 in K^0->pi^0vv calculations
      akp = hadron_data(48)     ! parameter kappa+ in K^+->pi^+vv calculations
      pc = hadron_data(49)      ! parameter P_c in K->pivv calculations
      del_ak0 = hadron_data(50) ! error of kappa0
      del_akp = hadron_data(51) ! error of kappa+
      del_pc = hadron_data(52)  ! error of P_c 
      eta_ed = hadron_data(53)  ! neutron EDM_d QCD coefficient
      eta_eu = hadron_data(54)  ! neutron EDM_u QCD coefficient
      eta_cd = hadron_data(55)  ! neutron CDM_d QCD coefficient
      eta_cu = hadron_data(56)  ! neutron CDM_u QCD coefficient
      eta_g = hadron_data(57)   ! neutron CDM_g QCD coefficient
      alamx = hadron_data(58)   ! neutron EDM chiral symmetry breaking scale
      ump(2) = hadron_data(59)  ! pole c quark mass (in B->X_s gamma and t->cH)
      br_tau_evv = hadron_data(60) ! Br(tau->evv)
      dmbp = hadron_data(61)    ! M_B_u
      rd = hadron_data(62)      ! Br(B->D tau nu)/Br(B->D l nu) in SM
      del_rd = hadron_data(63)  ! error of Br(B->D tau nu)/Br(B->D l nu) in SM
      rds = hadron_data(64)     ! Br(B->D* tau nu)/Br(B->D* l nu) in SM
      del_rds = hadron_data(65) ! error of Br(B->D* tau nu)/Br(B->D* l nu) in SM
      return
      end

      subroutine sflav_input(ilev,ierr)
c     Input read from file susy_flavor.in
      implicit double precision (a-h,o-z)
      character*18 block
      logical sldiag,sqdiag,find_block,read_var,read_mat_el
      double complex amg,amgg,amue
      double complex ls,ks,ds,es,us,ws
      double complex lms,rms,ums,dms,qms
      double complex cz,co,ci
      common/num/cz,co,ci,zero,one
      common/edm_qcd/eta_ed,eta_eu,eta_cd,eta_cu,eta_g,alamx
      common/kpivv/ak0,del_ak0,akp,del_akp,pc,del_pc,alam_k
      common/sm_4q/eta_cc,eta_ct,eta_tt,eta_b,bk_sm,bd_sm,bb_sm(2)
      common/meson_data/dmk,amk,epsk,fk,dmd,amd,fd,
     $    amb(2),dmb(2),tau_b(2),fb(2)
      common/dtau_data/dmbp,rd,del_rd,rds,del_rds
      common/bx_4q/bk(5),bd(5),bb(2,5),amu_k,amu_d,amu_b
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
      common/fmass/em(3),um(3),dm(3)
      common/qmass_pole/ump(3),dmp(3)
      common/tau_gam/br_tau_evv
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/msoft/lms(3,3),rms(3,3),ums(3,3),dms(3,3),qms(3,3)
      common/sf_cont/eps,indx(3,3),iconv
      common/sflav_data/control_data(3),sm_data(31),ckm_data(4),
     $     extpar_data(49),extim_data(49),hadron_data(65),tb_min(6)
c     reset error code
      ierr = 0
c     set default SM variables, clear SUSY parameters
      call sflav_defaults
c     open input file
      ifl = 1
      open(ifl,file='susy_flavor.in',status='old')
c     Find non-standard Block SOFTINP, check the status of soft parameters
      block = 'SOFTINP           '
      ibl = 3
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_var(ifl,control_data,ibl)) goto 11
         end do
      end if
c     iconv = 1: sfermion input parameters in SLHA2 conventions
c     iconv = 2: sfermion input parameters in conventions of hep-ph/9511250
 11   iconv = idnint(control_data(1)) ! soft parameters convention
      if ((iconv.ne.1).and.(iconv.ne.2)) stop
     $     'susy_flavor.in: incorrect iconv value'
c     input_type = 1: sfermion off-diagonal terms given as mass insertions
c                     LR diagonal terms given as dimensionless parameters
c     input_type = 2: fermion soft terms given as absolute values
c     for more information read comments in susy_flavor.in file
      input_type = idnint(control_data(2)) ! dimension of soft parameters
      if ((input_type.ne.1).and.(input_type.ne.2)) stop
     $     'susy_flavor.in: incorrect input_type value'
      ilev = idnint(control_data(3)) ! chiral resummation level
c     ilev = 0: no resummation
c     ilev = 1: analytical resummation in the decoupling limit
c     ilev = 2: iterative numerical resummation (default)
      if ((ilev.ne.0).and.(ilev.ne.1).and.(ilev.ne.2)) stop
     $     'susy_flavor.in: incorrect ilev value'
c     find SM input block, read the SM data
      block = 'SMINPUTS          '
      ibl = 31
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_var(ifl,sm_data,ibl)) goto 12
         end do
      end if
 12   alpha_z = 1/sm_data(1)    ! 1/alpha_em(M_Z)
      alpha_s = sm_data(3)      ! alpha_s(MZ)
      zm0 = sm_data(4)          ! M_Z
c     Fermion mass initialization, input: MSbar running quark masses
      bottom = sm_data(5)       ! mb(mb)
      bot_scale = bottom
      top_pole = sm_data(6)     ! mt(pole)
      top = 0.94d0*top_pole     ! mt(mt_pole) - running mt at mt(pole)
      top_scale = top
      em(3) = sm_data(7)        ! m_tau (pole)
      em(1) = sm_data(11)       ! m_e (pole)
      em(2) = sm_data(13)       ! m_mu (pole)
      dml(1) = sm_data(21)      ! md(2 GeV)
      uml(1) = sm_data(22)      ! mu(2 GeV)
      dml(2) = sm_data(23)      ! ms(2 GeV)
      uml(2) = sm_data(24)      ! mc(mc)
      amuu(2) = uml(2)          ! charm mass assumed to be given as mc(mc) 
      wm0 = sm_data(30)         ! M_W (not standard SLHA2!)
      st2 = sm_data(31)         ! sin^2(theta_W) in MSbar (not standard SLHA2!)
c     Electroweak and strong parameter initialization
      call vpar_update(zm0,wm0,alpha_z,st2) ! sets electroweak parameters
      call init_fermion_sector(alpha_s,top,top_scale,bottom,bot_scale)
c     find V_CKM input block, read the CKM data
      block = 'VCKMIN            '
      ibl = 4
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_var(ifl,ckm_data,ibl)) goto 13
         end do
      end if
 13   alam = ckm_data(1)        ! lambda
      apar = ckm_data(2)        ! A
      rhobar = ckm_data(3)      ! rho bar
      etabar = ckm_data(4)      ! eta bar
      call ckm_wolf(alam,apar,rhobar,etabar)
c     find MINPAR input block, read tan(beta) if defined (other entries
c     of MINPAR ignored)
      block = 'MINPAR            '
      ibl = 6
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_var(ifl,tb_min,ibl)) goto 14
         end do
      end if
 14   extpar_data(25) = tb_min(3) ! tan(beta) initialized in MINPAR
c     find EXTPAR input block, read the real Higgs and gaugino data
      block = 'EXTPAR            '
      ibl = 49
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_var(ifl,extpar_data,ibl)) goto 15
         end do
      end if
c 14   if (idnint(extpar_data(0)).ne.-1) stop
c     $     'susy_flavor.in: input has to be defined at EW scale!'
 15   amgg  = extpar_data(1)    ! M1 (bino mass, complex)
      amg   = extpar_data(2)    ! M2 (wino mass, complex)
      amglu = extpar_data(3)    ! M3 (gluino mass)
      amue  = extpar_data(23)   ! mu (complex)
      tanbe = extpar_data(25)   ! tan(beta)
      pm    = extpar_data(26)   ! M_A pole
      us(3,3)  = extpar_data(11) ! A_t
      ds(3,3)  = extpar_data(12) ! A_b
      ls(3,3)  = extpar_data(13) ! A_tau
      do i=1,3
         lms(i,i) = extpar_data(30+i)**2 ! M_eL, M_muL, M_tauL 
         rms(i,i) = extpar_data(33+i)**2 ! M_eR, M_muR, M_tauR 
         qms(i,i) = extpar_data(40+i)**2 ! M_q1L, M_q2L, M_q3L
         ums(i,i) = extpar_data(43+i)**2 ! M_uR, M_cR, M_tR
         dms(i,i) = extpar_data(46+i)**2 ! M_dR, M_sR, M_bR
      end do
c     SLHA2 EXTPAR extension (e.g. SuSpect): separate entries for 1st/2nd generation A-terms (optional)
      us(1,1) = extpar_data(14) ! A_u 
      ds(1,1) = extpar_data(15) ! A_d 
      ls(1,1) = extpar_data(16) ! A_e
      us(2,2) = extpar_data(17) ! A_c 
      ds(2,2) = extpar_data(18) ! A_s 
      ls(2,2) = extpar_data(19) ! A_mu
c     find optional IMEXTPAR input block, read the imaginary data
      block = 'IMEXTPAR          '
      ibl = 49
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_var(ifl,extim_data,ibl)) goto 16
         end do
 16      amgg = amgg + ci*extim_data(1) ! M1 (bino mass, complex)
         amg  = amg  + ci*extim_data(2) ! M2 (wino mass, complex)
         amue = amue + ci*extim_data(23) ! mu (complex)
         us(3,3) = us(3,3) + ci*extim_data(11) ! A_t
         ds(3,3) = ds(3,3) + ci*extim_data(12) ! A_b
         ls(3,3) = ls(3,3) + ci*extim_data(13) ! A_tau
c     SLHA2 IMEXTPAR extension: optional separate entries for 1st/2nd generation A-terms
         us(1,1) = us(1,1) + ci*extim_data(14) ! A_u
         ds(1,1) = ds(1,1) + ci*extim_data(15) ! A_d 
         ls(1,1) = ls(1,1) + ci*extim_data(16) ! A_e
         us(2,2) = us(2,2) + ci*extim_data(17) ! A_c
         ds(2,2) = ds(2,2) + ci*extim_data(18) ! A_s 
         ls(2,2) = ls(2,2) + ci*extim_data(19) ! A_mu
      end if
c     Higgs sector
      call init_higgs_sector(pm,tanbe,amue,ierr)
      if (ierr.ne.0) write(*,*) 'negative tree level Higgs mass^2?'
c     SUSY fermion sector
      call init_ino_sector(amgg,amg,amglu,amue,tanbe,ierr)
      if (ierr.ne.0) write(*,*) '-ino mass below M_Z/2?'
c     find sfermion soft input blocks (all are optional), read the sfermion data
      block = 'MSL2IN            '
      ibl = 6
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_mat_el(ifl,lms,0)) goto 17 ! left slepton mass M_SL^2, real part
         end do
      end if
 17   block = 'IMMSL2IN          '
      ibl = 3
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_mat_el(ifl,lms,1)) goto 18 ! left slepton mass M_SL^2, imaginary part
         end do
      end if
 18   block = 'MSE2IN            '
      ibl = 6
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_mat_el(ifl,rms,0)) goto 19 ! right slepton mass M_ER^2, real part
         end do
      end if
 19   block = 'IMMSE2IN          '
      ibl = 3
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_mat_el(ifl,rms,1)) goto 20 ! right slepton mass M_ER^2, imaginary part
         end do
      end if
 20   block = 'TEIN              '
      ibl = 9
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_mat_el(ifl,ls,0)) goto 21 ! slepton LR mixing A_L, real part
         end do
      end if
 21   block = 'IMTEIN            '
      ibl = 9
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_mat_el(ifl,ls,1)) goto 22 ! slepton LR mixing A_L, imaginary part
         end do
      end if
 22   block = 'MSQ2IN            '
      ibl = 6
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_mat_el(ifl,qms,0)) goto 23 ! left squark mass M_QL^2, real part
         end do
      end if
 23   block = 'IMMSQ2IN          '
      ibl = 3
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_mat_el(ifl,qms,1)) goto 24 ! left squark mass M_QL^2, imaginary part
         end do
      end if
 24   block = 'MSU2IN            '
      ibl = 6
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_mat_el(ifl,ums,0)) goto 25 ! right squark mass M_UR^2, real part
         end do
      end if
 25   block = 'IMMSU2IN          '
      ibl = 3
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_mat_el(ifl,ums,1)) goto 26 ! right squark mass M_UR^2, imaginary part
         end do
      end if
 26   block = 'MSD2IN            '
      ibl = 6
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_mat_el(ifl,dms,0)) goto 27 ! right squark mass M_DR^2, real part
         end do
      end if
 27   block = 'IMMSD2IN          '
      ibl = 3
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_mat_el(ifl,dms,1)) goto 28 ! right squark mass M_DR^2, imaginary part
         end do
      end if
 28   block = 'TUIN              '
      ibl = 9
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_mat_el(ifl,us,0)) goto 29 ! up-squark LR mixing A_U, real part
         end do
      end if
 29   block = 'IMTUIN            '
      ibl = 9
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_mat_el(ifl,us,1)) goto 30 ! up-squark LR mixing A_U, imaginary part
         end do
      end if
 30   block = 'TDIN              '
      ibl = 9
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_mat_el(ifl,ds,0)) goto 31 ! down-squark LR mixing A_D, real part
         end do
      end if
 31   block = 'IMTDIN            '
      ibl = 9
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_mat_el(ifl,ds,1)) goto 32 ! down-squark LR mixing A_D, imaginary part
         end do
      end if
 32   block = 'TEINH             '
      ibl = 9
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_mat_el(ifl,ks,0)) goto 33 ! slepton LR mixing A_L', real part
         end do
      end if
 33   block = 'IMTEINH           '
      ibl = 9
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_mat_el(ifl,ds,1)) goto 34 ! slepton LR mixing A_L', imaginary part
         end do
      end if
 34   block = 'TDINH             '
      ibl = 9
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_mat_el(ifl,es,0)) goto 35 ! down-squark LR mixing A_D', real part
         end do
      end if
 35   block = 'IMTDINH           '
      ibl = 9
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_mat_el(ifl,es,1)) goto 36 ! down-squark LR mixing A_D', imaginary part
         end do
      end if
 36   block = 'TUINH             '
      ibl = 9
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_mat_el(ifl,ws,0)) goto 37 ! up-squark LR mixing A_U', real part
         end do
      end if
 37   block = 'IMTUINH           '
      ibl = 9
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_mat_el(ifl,ws,1)) goto 38 ! up-squark LR mixing A_U', imaginary part
         end do
      end if
c     remove instabilities by adding tiny mass splitting
 38   do i=1,3
         lms(i,i) = (1 + eps*i)*lms(i,i)
         rms(i,i) = (1 - eps*i)*rms(i,i)
         qms(i,i) = (1 + eps*i)*dble(qms(i,i))
         ums(i,i) = (1 - eps*i)*dble(ums(i,i))
         dms(i,i) = (1 - eps*i)*dble(dms(i,i))
      end do
c     expand delta parameters to full mass entries
      if (input_type.eq.1) then
         do k=1,2
            do l=k+1,3
               lms(k,l) = lms(k,l)*sqrt(lms(k,k)*lms(l,l))
               rms(k,l) = rms(k,l)*sqrt(rms(k,k)*rms(l,l))
               qms(k,l) = qms(k,l)*sqrt(qms(k,k)*qms(l,l))
               ums(k,l) = ums(k,l)*sqrt(ums(k,k)*ums(l,l))
               dms(k,l) = dms(k,l)*sqrt(dms(k,k)*dms(l,l))
            end do
         end do
      end if
c     initialize h.c. LL and RR entries
      do k=1,2
         do l=k+1,3
            lms(l,k) = dconjg(lms(k,l))
            rms(l,k) = dconjg(rms(k,l))
            qms(l,k) = dconjg(qms(k,l))
            ums(l,k) = dconjg(ums(k,l))
            dms(l,k) = dconjg(dms(k,l))
         end do
      end do
c     expand LR mixing from dimensionless to dimensionful
      if (input_type.eq.1) then
         do k=1,3
            do l=1,3
               ls(k,l) = ls(k,l)*abs(lms(k,k)*rms(l,l))**0.25d0
               ds(k,l) = ds(k,l)*abs(qms(k,k)*dms(l,l))**0.25d0
               us(k,l) = us(k,l)*abs(qms(k,k)*ums(l,l))**0.25d0
               ks(k,l) = ks(k,l)*abs(lms(k,k)*rms(l,l))**0.25d0
               es(k,l) = es(k,l)*abs(qms(k,k)*dms(l,l))**0.25d0
               ws(k,l) = ws(k,l)*abs(qms(k,k)*ums(l,l))**0.25d0
            end do
         end do
      end if
c     if slepton input data in SLHA format, rewrite them to
c     hep-ph/9511250 convention
      if (iconv.eq.1) then
         call sl_slha_to_jr
         call sq_slha_to_jr
      end if
c     slepton diagonalization routine
      if (sldiag()) then
         write(*,*) 'negative tree level slepton mass^2?'
         ierr = 7
      end if
c     squark diagonalization routine
      if (sqdiag()) then
         write(*,*) 'negative tree level squark mass^2?'
         ierr = 14
      end if
c     reset status of physical Higgs mass after parameter changes
      call reset_phys_data
c     find optional hadronic data input block, read the data
      block = 'SFLAV_HADRON      '
      ibl = 65
      if (find_block(ifl,block)) then
         do i=1,ibl
            if (read_var(ifl,hadron_data,ibl)) goto 39
         end do
      end if
 39   fk = hadron_data(1)       ! f_K
      fd = hadron_data(2)       ! f_D
      fb(1) = hadron_data(3)    ! f_B_d
      fb(2) = hadron_data(4)    ! f_B_s
      bk_sm = hadron_data(5)    ! B_K for SM contribution to KKbar
      eta_cc = hadron_data(6)   ! eta_cc in KK mixing (SM)
      eta_ct = hadron_data(7)   ! eta_cc in KK mixing (SM)
      eta_tt = hadron_data(8)   ! eta_cc in KK mixing (SM)
      amu_k = hadron_data(9)    ! scale for B_K (non-SM)    
      do i=1,5
         bk(i) = hadron_data(9+i) ! B_K for VLL, SLL(2), LR(2) 
      end do 
      bd_sm = hadron_data(15)   ! B_D for SM contribution 
      amu_d = hadron_data(16)   ! scale for B_D (non-SM)
      do i=1,5
         bd(i) = hadron_data(16+i) ! B_D for VLL, SLL(2), LR(2) 
      end do 
      bb_sm(1) = hadron_data(22) ! B_Bd for SM contribution 
      amu_b = hadron_data(23)   ! scale for B_B (non-SM, both Bd and Bs)
      do i=1,5
         bb(1,i) = hadron_data(23+i) ! B_Bd for VLL, SLL(2), LR(2) 
      end do 
      bb_sm(2) = hadron_data(29) ! B_Bs for SM contribution 
      eta_b = hadron_data(30)   ! eta_b for  BsBsbar (SM)
      do i=1,5
         bb(2,i) = hadron_data(30+i) ! B_Bs for VLL, SLL(2), LR(2) 
      end do 
      tau_b(1) = hadron_data(36) ! Bd lifetime (experimental)
      tau_b(2) = hadron_data(37) ! Bs lifetime (experimental)
      amb(1) = hadron_data(38)  ! Bd mass (experimental)
      amb(2) = hadron_data(39)  ! Bs mass (experimental)
      dmb(1) = hadron_data(40)  ! Delta Bd (experimental)
      dmb(2) = hadron_data(41)  ! Delta Bs (experimental)
      amk = hadron_data(42)     ! K0 mass (experimental)
      dmk = hadron_data(43)     ! Delta mK (experimental)
      epsk = hadron_data(44)    ! eps_K (experimental)
      amd = hadron_data(45)     ! D0 mass (experimental)
      dmd = hadron_data(46)     ! Delta mD (experimental)
      ak0 = hadron_data(47)     ! parameter kappa0 in K^0->pi^0vv calculations
      akp = hadron_data(48)     ! parameter kappa+ in K^+->pi^+vv calculations
      pc = hadron_data(49)      ! parameter P_c in K->pivv calculations
      del_ak0 = hadron_data(50) ! error of kappa0
      del_akp = hadron_data(51) ! error of kappa+
      del_pc = hadron_data(52)  ! error of P_c 
      eta_ed = hadron_data(53)  ! neutron EDM_d QCD coefficient
      eta_eu = hadron_data(54)  ! neutron EDM_u QCD coefficient
      eta_cd = hadron_data(55)  ! neutron CDM_d QCD coefficient
      eta_cu = hadron_data(56)  ! neutron CDM_u QCD coefficient
      eta_g = hadron_data(57)   ! neutron CDM_g QCD coefficient
      alamx = hadron_data(58)   ! neutron EDM chiral symmetry breaking scale
      ump(2) = hadron_data(59)  ! pole c quark mass (in B->X_s gamma and t->cH)
      br_tau_evv = hadron_data(60) ! Br(tau->evv)
      dmbp = hadron_data(61)    ! M_B_u
      rd = hadron_data(62)      ! Br(B->D tau nu)/Br(B->D l nu) in SM
      del_rd = hadron_data(63)  ! error of Br(B->D tau nu)/Br(B->D l nu) in SM
      rds = hadron_data(64)     ! Br(B->D* tau nu)/Br(B->D* l nu) in SM
      del_rds = hadron_data(65) ! error of Br(B->D* tau nu)/Br(B->D* l nu) in SM

c     Neutral CP-even Higgs masses with the approximate 2-loop formulae 
c     (hep-ph/9903404, hep-ph/9307201)
      call mhcorr_app2(ierr)
      if (ierr.ne.0) write(*,*) 'error in 2-loop CP-even Higgs mass^2?'
      close(ifl)
      return
      end

      subroutine sflav_output(ilev,ierr)
c     writes SUSY_FLAVOR results to susy_flavor.out
      implicit double precision (a-h,o-z)
      double complex gm2,gm3
      double complex zu,zd,zn,zpos,zneg,zv,zl
      dimension corr_l(3),corr_d(3),corr_u(3),corr_ckm(3,3)
      common/hmass/cm(2),rm(2),ppm(2),zr(2,2),zh(2,2)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/neut/fnm(4),zn(4,4)
      common/fmass/em(3),um(3),dm(3)
      common/gmass/gm1,gm2,gm3
      common/hmass_EPA/pm_epa,hm1,hm2,sa,ca,sb,cb
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sflav_df0/edmn,edml(3),gminus2(3)
      common/sflav_df1/br_mu_egamma,br_tau_egamma,br_tau_mugamma,
     $     br_k0,br_kp,br_taunu,dtaunu_ratio,dstaunu_ratio,bxgamma,
     $     br_bdll(3,3),br_bsll(3,3),br_tuh,br_tch
      common/sflav_df2/eps_k,delta_mk,delta_md,
     $     delta_mbd,dmbd_re,dmbd_im,
     $     delta_mbs,dmbs_re,dmbs_im

      ifl = 1                   ! output file number
      open(ifl,file='susy_flavor.out',status='unknown')

      write(ifl,*)'#'
      write(ifl,*)'#         ***************************'
      write(ifl,*)'#         * SUSY_FLAVOR 2.50 output *'
      write(ifl,*)'#         ***************************'
      write(ifl,*)'#'


      write(ifl,*)'BLOCK SFLAV_CONTROL'
      write(ifl,101) 1, ilev, 'resummation level of chiral corrections'
      write(ifl,101) 2, ierr, 
     $     'error code (0 if all calculations were correct)'

      write(ifl,*)'BLOCK SFLAV_MASS        # Mass Spectrum'
      write(ifl,*)'#     code          mass          # particle'
      write(ifl,100) 24, wm, 'W+'
      write(ifl,100) 25, hm2, 'h (simple 2-loop approximation only)' 
      write(ifl,100) 35, hm1, 'H (simple 2-loop approximation only)'
      write(ifl,100) 36, ppm(1), 'A'
      write(ifl,100) 37, cm(1), 'H+'
      write(ifl,100) 41, em(1), 'e (pole)'
      write(ifl,100) 42, em(2), 'mu (pole)'
      write(ifl,100) 43, em(3), 'tau (pole)'
      write(ifl,100) 44, dm(1), 'md(mt) (running)'
      write(ifl,100) 45, dm(2), 'ms(mt) (running)'
      write(ifl,100) 46, dm(3), 'mb(mt) (running)'
      write(ifl,100) 47, um(1), 'mu(mt) (running)'
      write(ifl,100) 48, um(2), 'mc(mt) (running)'
      write(ifl,100) 49, um(3), 'mt(mt) (running)'
      write(ifl,100) 1000021, gm1, '~g'
      write(ifl,100) 1000022, fnm(1), '~chi_10'
      write(ifl,100) 1000023, fnm(2), '~chi_20'
      write(ifl,100) 1000025, fnm(3), '~chi_30'
      write(ifl,100) 1000035, fnm(4), '~chi_40'
      write(ifl,100) 1000024, fcm(1), '~chi_1+'
      write(ifl,100) 1000037, fcm(2), '~chi_2+'
      write(ifl,*)'# sfermion mass eigenstates'
      write(ifl,100) 101, sdm(1), '~d(1)'
      write(ifl,100) 102, sdm(2), '~d(2)'
      write(ifl,100) 103, sdm(3), '~d(3)'
      write(ifl,100) 104, sdm(4), '~d(4)'
      write(ifl,100) 105, sdm(5), '~d(5)'
      write(ifl,100) 106, sdm(6), '~d(6)'
      write(ifl,100) 111, sum(1), '~u(1)'
      write(ifl,100) 112, sum(2), '~u(2)'
      write(ifl,100) 113, sum(3), '~u(3)'
      write(ifl,100) 114, sum(4), '~u(4)'
      write(ifl,100) 115, sum(5), '~u(5)'
      write(ifl,100) 116, sum(6), '~u(6)'
      write(ifl,100) 121, slm(1), '~l(1)'
      write(ifl,100) 122, slm(2), '~l(2)'
      write(ifl,100) 123, slm(3), '~l(3)'
      write(ifl,100) 124, slm(4), '~l(4)'
      write(ifl,100) 125, slm(5), '~l(5)'
      write(ifl,100) 126, slm(6), '~l(6)'
      write(ifl,100) 131, vm(1), '~nu(1)'
      write(ifl,100) 132, vm(2), '~nu(2)'
      write(ifl,100) 133, vm(3), '~nu(3)'

c     estimate chiral corrections size
      call chiral_corr_size(corr_l,corr_d,corr_u,corr_ckm)

      write(ifl,*)'BLOCK SFLAV_CHIRAL_YUKAWA     '//
     $     '# Chiral corrections size to Yukawa couplings'
      write(ifl,100) 1, corr_l(1), 'correction to Y_e'
      write(ifl,100) 2, corr_l(2), 'correction to Y_mu'
      write(ifl,100) 3, corr_l(3), 'correction to Y_tau'
      write(ifl,100) 4, corr_d(1), 'correction to Y_d'
      write(ifl,100) 5, corr_d(2), 'correction to Y_s'
      write(ifl,100) 6, corr_d(3), 'correction to Y_b'
      write(ifl,100) 7, corr_u(1), 'correction to Y_u'
      write(ifl,100) 8, corr_u(2), 'correction to Y_c'
      write(ifl,100) 9, corr_u(3), 'correction to Y_t'

      write(ifl,*)'BLOCK SFLAV_CHIRAL_CKM     '//
     $     '# Chiral corrections size to CKM matrix'
      write(ifl,102) 1,1, corr_ckm(1,1), 'correction to V_11'
      write(ifl,102) 1,2, corr_ckm(1,2), 'correction to V_12'
      write(ifl,102) 1,3, corr_ckm(1,3), 'correction to V_13'
      write(ifl,102) 2,1, corr_ckm(2,1), 'correction to V_21'
      write(ifl,102) 2,2, corr_ckm(2,2), 'correction to V_22'
      write(ifl,102) 2,3, corr_ckm(2,3), 'correction to V_23'
      write(ifl,102) 3,1, corr_ckm(3,1), 'correction to V_31'
      write(ifl,102) 3,2, corr_ckm(3,2), 'correction to V_32'
      write(ifl,102) 3,3, corr_ckm(3,3), 'correction to V_33'

      write(ifl,*)'BLOCK SFLAV_DELTA_F0     # Delta F = 0 processes'
      write(ifl,100) 1, edml(1), 'EDM_e'
      write(ifl,100) 2, edml(2), 'EDM_mu'
      write(ifl,100) 3, edml(3), 'EDM_tau'
      write(ifl,100) 4, edmn, 'neutron EDM'
      write(ifl,100) 5, gminus2(1), '(g-2)_e/2, SUSY contribution'
      write(ifl,100) 6, gminus2(2), '(g-2)_mu/2, SUSY contribution'
      write(ifl,100) 7, gminus2(3), '(g-2)_tau/2, SUSY contribution'

      write(ifl,*)'BLOCK SFLAV_DELTA_F1     # Delta F = 1 processes'
      write(ifl,100) 1, br_mu_egamma, 'Br(mu-> e gamma)'
      write(ifl,100) 2, br_tau_egamma, 'Br(tau-> e gamma)'
      write(ifl,100) 3, br_tau_mugamma, 'Br(tau-> mu gamma)'
      write(ifl,100) 4, br_k0, 'Br(K0 -> pi0 nu nu)'
      write(ifl,100) 5, br_kp, 'Br(K+ -> pi+ nu nu)'
      write(ifl,100) 6, br_taunu, 'BR(B -> tau nu)'
      write(ifl,100) 7, dtaunu_ratio, 
     $     'BR(B -> D tau nu)/BR(B -> D l nu)'
      write(ifl,100) 8, dstaunu_ratio, 
     $     'BR(B -> D* tau nu)/BR(B -> D* l nu)'
      write(ifl,100) 9, bxgamma, 'BR(B -> X_s gamma)'
      write(ifl,100) 10, br_tuh, 'BR(t -> u h)'
      write(ifl,100) 11, br_tch, 'BR(t -> c h)'
      write(ifl,100) 12, br_bdll(1,1), 'BR(B_d -> e e)'
      write(ifl,100) 13, br_bdll(2,2), 'BR(B_d -> mu mu)'
      write(ifl,100) 14, br_bdll(3,3), 'BR(B_d -> tau tau)'
      write(ifl,100) 15, 2*br_bdll(2,1), 'BR(B_d -> mu e)'
      write(ifl,100) 16, 2*br_bdll(3,1), 'BR(B_d -> tau e)'
      write(ifl,100) 17, 2*br_bdll(3,2), 'BR(B_d -> tau mu)'
      write(ifl,100) 18, br_bsll(1,1), 'BR(B_s -> e e)'
      write(ifl,100) 19, br_bsll(2,2), 'BR(B_s -> mu mu)'
      write(ifl,100) 20, br_bsll(3,3), 'BR(B_s -> tau tau)'
      write(ifl,100) 21, 2*br_bsll(2,1), 'BR(B_s -> mu e)'
      write(ifl,100) 22, 2*br_bsll(3,1), 'BR(B_s -> tau e)'
      write(ifl,100) 23, 2*br_bsll(3,2), 'BR(B_s -> tau mu)'

      write(ifl,*)'BLOCK SFLAV_DELTA_F2     # Delta F = 2 processes'
      write(ifl,100) 1, eps_k, 'epsilon_K'
      write(ifl,100) 2, delta_mk, 'Delta m_K (GeV)'
      write(ifl,100) 3, delta_md, 'Delta m_D (GeV)'
      write(ifl,100) 4, delta_mbd, 'Delta m_Bd (GeV)'
      write(ifl,100) 5, dmbd_re, 'Re(H_eff_Bd)'
      write(ifl,100) 6, dmbd_im, 'Im(H_eff_Bd)'
      write(ifl,100) 7, delta_mbs, 'Delta m_Bs (GeV)'
      write(ifl,100) 8, dmbs_re, 'Re(H_eff_Bs)'
      write(ifl,100) 9, dmbs_im, 'Im(H_eff_Bs)'

      close(ifl)

 100  format(2x,i7,5x,1pe16.9,4x,'#',1x,a)
 101  format(3x,i2,3x,i2,22x,'#',1x,a)
 102  format(3x,i2,2x,i2,5x,1pe16.9,4x,'#',1x,a)
 
      return
      end

      subroutine cr_mat_print(a,n,ifl)
c     print real part of complex matrix
      double complex a(n,n)
      do i=1,n
         write(ifl,'(100(1pe10.3,1x))')(dble(a(i,j)),j=1,n)
      end do
      write(ifl,*)
      return
      end

      subroutine ci_mat_print(a,n,ifl)
c     print imaginary part of complex matrix
      double complex a(n,n)
      do i=1,n
         write(ifl,'(100(1pe10.3,1x))')(dimag(a(i,j)),j=1,n)
      end do
      write(ifl,*)
      return
      end
