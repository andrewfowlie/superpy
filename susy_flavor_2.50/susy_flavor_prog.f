c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: SUSY_FLAVOR_PROG.F
c     Released: 25:10:2013(J.R.)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Driver file for SUSY_FLAVOR library                             c
c     Example of MSSM parameter initialization inside the program     c
c     Test output for SUSY spectrum and implemented rare decays       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program susy_flavor_prog
      implicit double precision (a-h,o-z)
      dimension sll(3),slr(3),amsq(3),amsu(3),amsd(3)
      double complex slmi_l(3),slmi_r(3),slmi_lr(3,3),slmi_lrp(3,3)
      double complex sqmi_l(3),sdmi_r(3),sumi_r(3)
      double complex sdmi_lr(3,3),sumi_lr(3,3)
      double complex sdmi_lrp(3,3),sumi_lrp(3,3)
      double complex amg,amgg,amue
      common/sf_cont/eps,indx(3,3),iconv
c     common/ph_units/hbar,gev_cm,gev_s
      dimension corr_l(3),corr_d(3),corr_u(3),corr_ckm(3,3)
      
c     Parameters defined inside the code.  

c     Input parameters convention choice
      iconv = 1                 ! SLHA2 input conventions
c     iconv = 2                 ! hep-ph/9511250 input conventions
      
c     fixes the treatment of enhanced chiral correction resummation
c     ilev = 0                  ! no resummation, SUSY corrections strictly 1-loop
c     ilev = 1                  ! resummation using the decoupling limit
      ilev = 2                  ! exact iterative solution, may not always converge
      
c     choose MSSM sectors to include
      ih = 1                    ! Higgs + gauge diagrams included
      ic = 1                    ! chargino diagrams included
      in = 1                    ! neutralino diagrams included
      ig = 1                    ! gluino diagrams included

      call set_active_sector(ih,ic,in,ig) ! set control variables

      call sflav_sm             ! initialize auxiliary SM parameters

c     SM basic input initialization
      zm0 = 91.1876d0           ! M_Z
      wm0 = 80.398d0            ! M_W
      alpha_z = 1/127.934d0     ! alpha_em(M_Z)
      st2_new = 0.23116d0       ! s_W^2 (here MSbar value)
      call vpar_update(zm0,wm0,alpha_z,st2_new)

c     QCD parameters and fermion mass initialization
      alpha_s = 0.1172d0        ! alpha_s(MZ)
      top_scale = 163.09d0
      top = 163.09d0            ! MSbar running m_t(top_scale)
      bot_scale = 4.18d0
      bot = 4.18d0              ! MSbar running m_b(bot_scale)
      call init_fermion_sector(alpha_s,top,top_scale,bot,bot_scale)

c     CKM matrix initialization
      alam = 0.2258d0           ! lambda
      apar = 0.808d0            ! A
      rhobar = 0.177d0          ! rho bar
      etabar = 0.36d0           ! eta bar
      call ckm_wolf(alam,apar,rhobar,etabar)

c     Higgs sector parameters
      pm    = 200               ! M_A
      tanbe = 4                 ! tan(beta)
      amue  = (200.d0,100.d0)   ! mu
      call init_higgs_sector(pm,tanbe,amue,ierr)
      if (ierr.ne.0) stop 'negative tree level Higgs mass^2?'
      
c     Gaugino sector parameters. CAUTION: if M1 is set to 0 here then
c     program sets M1 and M2 GUT-related, i.e. M1 = 5/3 s_W^2/c_W^2*M2
      amgg  = (200.d0,0.d0)     ! M1 (bino mass)
      amg   = (300.d0,0.d0)     ! M2 (wino mass)
      amglu = 600               ! M3 (gluino mass)
      call init_ino_sector(amgg,amg,amglu,amue,tanbe,ierr)
      if (ierr.ne.0) write(*,*) '-ino mass below M_Z/2?'
      
c     Slepton diagonal soft breaking parameters
      sll(1) = 300.d0           ! left selectron mass scale
      sll(2) = 300.d0           ! left smuon mass scale
      sll(3) = 300.d0           ! left stau mass scale
      slr(1) = 300.d0           ! right selectron mass scale
      slr(2) = 300.d0           ! right smuon mass scale
      slr(3) = 300.d0           ! right stau mass scale
c     Slepton LL and RR mass insertions (hermitian matrices)
c     slmi_x(1),slmi_x(2), slmi_x(3) are 12,23,31 entry respectively
      do i=1,3
         slmi_l(i) = dcmplx(0.d0,0.d0) ! slepton LL mass insertion
         slmi_r(i) = dcmplx(0.d0,0.d0) ! slepton RR mass insertion
      end do 
      slmi_l(1) = (2.d-2,3.d-2) ! example, non-vanishing LL 12 entry
c     Slepton LR mass insertions, non-hermitian in general
c     All entries dimensionless (normalized to diagonal masses)
      do i=1,3
         do j=1,3
c     holomorphic LR mixing terms
            slmi_lr(i,j) = (0.d0,0.d0)
c     non-holomorphic LR mixing terms
            slmi_lrp(i,j) = (0.d0,0.d0)
         end do
      end do 
c     Example: diagonal entries normalized to Y_l as in SUGRA  
      slmi_lr(1,1) = (1.d-4,0.d0) ! A_e
      slmi_lr(2,2) = (1.d-2,0.d0) ! A_mu
      slmi_lr(3,3) = (1.d-1,0.d0) ! A_tau
      slmi_lr(2,3) = (2.d-2,1.d-2) ! example, non-vanishing LR 23 entry
c     Calculate physical masses and mixing angles
      call init_slepton_sector(sll,slr,slmi_l,slmi_r,slmi_lr,slmi_lrp
     $     ,ierr)
      if (ierr.ne.0) stop 'negative tree level slepton mass^2?'
      
c     Squark diagonal soft breaking parameters
      amsq(1) = 500.d0          ! left squark mass, 1st generation
      amsq(2) = 450.d0          ! left squark mass, 2nd generation
      amsq(3) = 400.d0          ! left squark mass, 3rd generation
      amsd(1) = 550.d0          ! right down squark mass
      amsd(2) = 550.d0          ! right strange squark mass
      amsd(3) = 300.d0          ! right sbottom mass
      amsu(1) = 450.d0          ! right up squark mass
      amsu(2) = 450.d0          ! right charm squark mass
      amsu(3) = 200.d0          ! right stop mass
c     Squark LL and RR mass insertions (hermitian matrices)
c     sqmi_l(1),sqmi_l(2), sqmi_l(3) are 12,23,31 entry respectively, etc.
      do i=1,3
         sqmi_l(i) = (0.d0,0.d0) ! squark LL mass insertion
         sumi_r(i) = (0.d0,0.d0) ! up-squark RR mass insertion
         sdmi_r(i) = (0.d0,0.d0) ! down-squark RR mass insertion
      end do 
      sqmi_l(2) = (-1.d-2,1.d-2) ! example, non-vanishing LL 23 entry
c     Squark holomorphic LR mass insertions, non-hermitian in general
c     All entries dimensionless (normalized to masses)
      do i=1,3
         do j=1,3
c     holomorphic LR mixing terms
            sumi_lr(i,j) = (0.d0,0.d0) ! up-squark 
            sdmi_lr(i,j) = (0.d0,0.d0) ! down-squark 
c     non-holomorphic LR mixing terms
            sumi_lrp(i,j) = (0.d0,0.d0) ! up-squark
            sdmi_lrp(i,j) = (0.d0,0.d0) ! down-squark
         end do
      end do 
c     Example: diagonal entries normalized to Y_d,Y_u as in SUGRA  
      sumi_lr(1,1) = dcmplx(1.d-5,0.d0)
      sumi_lr(2,2) = dcmplx(4.d-3,0.d0)
      sumi_lr(3,3) = dcmplx(1.d0,0.d0)
      sdmi_lr(1,1) = dcmplx(-1.d-4,0.d0)
      sdmi_lr(2,2) = dcmplx(-2.d-3,0.d0)
      sdmi_lr(3,3) = dcmplx(-8.d-2,0.d0)
      sdmi_lr(2,3) = (1.d-2,-1.d-2) ! example, non-vanishing down LR 23 entry
c     Calculate physical masses and mixing angles
      call init_squark_sector(amsq,amsu,amsd,sqmi_l,sumi_r,sdmi_r,
     $     sumi_lr,sdmi_lr,sumi_lrp,sdmi_lrp,ierr)
      if (ierr.ne.0) stop 'negative tree level squark mass^2?'

c     reset status of physical Higgs mass after parameter changes
      call reset_phys_data

c     Neutral CP-even Higgs masses with the 2-loop approximate formula
c     (see SUS_FlAVOR v2.5 manual Section 6)
      call mhcorr_app2(ierr)
      if (ierr.ne.0) stop 'negative CP-even Higgs mass^2?'

c     !!! End of input section !!!

c     perform resummation of chirally enhanced corrections
      call set_resummation_level(ilev,ierr)
      if (ierr.ne.0) write(*,*)ierr,
     $     'Error in chiral corrections resummation!'

c     print output to susy_flavor.out

      call susy_flavor          ! main routine calculating physical observables
      call sflav_output(ilev,ierr) ! output written to susy_flavor.out

c     repeat calculations and print them to screen

      call chiral_corr_size(corr_l,corr_d,corr_u,corr_ckm)
 
      write(*,96)'Corrections to lepton Yukawa couplings:     ',
     $     (100*corr_l(i),i=1,3)
      write(*,96)'Corrections to down-quark Yukawa couplings: ',
     $     (100*corr_d(i),i=1,3)
      write(*,96)'Corrections to up-quark Yukawa couplings:   ',
     $     (100*corr_u(i),i=1,3)
      write(*,96)'Corrections to CKM matrix elements: '
      do i=1,3
         write(*,96)' ',(100*corr_ckm(i,j),j=1,3)
      end do
      write(*,*)
 96   format(a,3(f8.1,"%",1x))

c     Results for implemented observables:
      write(*,99)'Electric dipole moments:'
      write(*,99)'Electron EDM = ',edm_l(1)
      write(*,99)'Muon EDM     = ',edm_l(2)
      write(*,99)'Tau EDM      = ',edm_l(3)
      write(*,99)'Neutron EDM  = ',edm_n()
      write(*,*)

      write(*,99)'g-2 anomaly, SUSY contribution:'
      write(*,99)'Electron (g-2)/2 = ',g_minus_2_anomaly(1)
      write(*,99)'Muon (g-2)/2     = ',g_minus_2_anomaly(2)
      write(*,99)'Tau (g-2)/2      = ',g_minus_2_anomaly(3)
      write(*,*)

      write(*,99)'l^J->l^I gamma decays:'
      write(*,99)'Br(mu-> e gamma)   = ',br_llg(2,1)
      write(*,99)'Br(tau-> e gamma)  = ',br_llg(3,1)
      write(*,99)'Br(tau-> mu gamma) = ',br_llg(3,2)
      write(*,*)

      write(*,99)'Neutrino K decays:'
      call k_pivv(br_k0,br_kp)
      write(*,99)'BR(K_L^0 -> pi^0 vv) = ',br_k0
      write(*,99)'BR(K^+   -> pi^+ vv) = ',br_kp
      write(*,*)

      write(*,99)'Leptonic B decays:'
      write(*,99)'BR(B_d -> mu^+ mu^-) = ',b_ll(3,1,2,2)
      write(*,99)'BR(B_s -> mu^+ mu^-) = ',b_ll(3,2,2,2)
      
      write(*,99)'BR(B_s -> mu^+ e^-)  = ',b_ll(3,2,2,1)
      call b_taunu(br_taunu,dtaunu_ratio,dstaunu_ratio)
      write(*,99)'BR(B -> tau nu)    = ',br_taunu
      write(*,99)'BR(B -> D tau nu)/BR(B -> D l nu) = ',dtaunu_ratio
      write(*,99)'BR(B -> D* tau nu)/BR(B -> D* l nu) = ',dstaunu_ratio
      write(*,*)

c     Physical quantities for BR(B->X_s g) calculation
      delb = 0.99d0             ! Photon energy infrared cutoff
      amiu_b= 4.8d0             ! Renormalization scale miu_b
      write(*,99)'BR(B -> X_s gamma) = ',bxg_nl(delb,amiu_b)
      write(*,*)

      write(*,99)'KK mixing:'
      call dd_kaon(eps_k,delta_mk)
      write(*,99)'eps_K     = ',eps_k
      write(*,99)'Delta m_K = ',delta_mk
c     write(*,99)'Delta m_K = ',delta_mk*gev_s/1.d12 ! in ps^-1
      write(*,*)

      write(*,99)'DD mixing:'
      call uu_dmeson(delta_md)
      write(*,99)'Delta m_D = ',delta_md
c     write(*,99)'Delta m_D = ',delta_md*gev_s/1.d12 ! in ps^-1
      write(*,*)

      write(*,99)'BB mixing:'
      call dd_bmeson(1,delta_mbd,dmb_re,dmb_im)
      write(*,99)'Re(H_eff_Bd) = ',dmb_re
      write(*,99)'Im(H_eff_Bd) = ',dmb_im
      write(*,99)'Delta m_B_d  = ',delta_mbd
c     write(*,99)'Delta m_B_d  = ',delta_mbd*gev_s/1.d12 ! in ps^-1
      call dd_bmeson(2,delta_mbs,dmb_re,dmb_im)
      write(*,99)'Delta m_B_s  = ',delta_mbs
c     write(*,99)'Delta m_B_s  = ',delta_mbs*gev_s/1.d12 ! in ps^-1
      write(*,*)

 99   format(a,6(1pe11.4,1x))
      end

