c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: MH_INIT.F
c     Revised: 15:05:1996(J.R.)
c     block data init_mh splitted into init_phys, init_const
c     and init_control
c     Revised: 20:02:2008
c     Initialization routines for SUSY sectors added
c     Revised: 20:02:2012
c     Normalization of LR soft terms changed

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     File contains initialization of standard masses and couplings,    c
c     control variables and some auxiliary initialization procedures    c
c     It also contains initialization routines for various SUSY sectors c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine init_run_qmass
c     running s,b masses at m_t scale
      implicit double precision (a-h,o-z)
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
      umu(3) = uml(3) + 1.d-5
      do i=1,3
         dmu(i) = qmass_nlo(dml(i),amud(i),umu(3))
      end do
      do i=1,2
         umu(i) = qmass_nlo(uml(i),amuu(i),umu(3))
      end do
      return
      end

      subroutine init_fermion_sector(alpha_s,tm,tscale,bm,bscale)
      implicit double precision (a-h,o-z)
      logical higgs, fermion
      common/fmass/em(3),um(3),dm(3)
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
      common/qmass_pole/ump(3),dmp(3)
      common/required_init/higgs,fermion
      external init_phys,init_control
c     Initialization of the running fermion masses. Masses of the light
c     quarks at scale mu=2 GeV are set in block data init_sm.  mc is
c     given at the scale mc. mt(mt) and mb(mb) are initialized there and
c     overwritten here by setting mt(tscale) and mb(bscale). Routine
c     calculates all quark masses at mu=tscale and rewrites into
c     common/fmass/
      amuu(3) = tscale
      uml(3) = tm               ! top mass at tscale
      amud(3) = bscale
      dml(3) = bm               ! bottom mass at bscale
      do iter=1,3
         call lam_fit(alpha_s)  ! fits Lambda_QCD at 3 loop level
         call lam_fit_nlo(alpha_s) ! fits Lambda_QCD at NLO level
         call init_run_qmass
         do i=1,3
            dm(i) = dmu(i)
            um(i) = umu(i)
         end do
         call qstep_update      ! update ordered msbar quark mass array
      end do
      ump(3) = qm_pole_nlo(uml(3)) ! NLO pole top mass
      dmp(3) = qm_pole_nlo(dml(3)) ! NLO pole bottom mass
c     if Higgs sector is initialized, calculate also Yukawa couplings
      if (higgs) call init_tree_yukawa
      fermion = .true.
      return
      end

      subroutine init_higgs_sector(pm,tanb,mu,ierr)
c     Tree-level initialization of SUSY Higgs sector
c     ierr=1     One or more scalar mass squares <= 0 (unstable Higgs potential)
      implicit double precision (a-h,o-z)
      double complex hmu,mu
      logical hdiag,higgs,fermion
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/hpar/hm1,hm2,hm12,hmu
      common/vev/v1,v2
      common/required_init/higgs,fermion
      external init_control
c     store mu in common
      hmu = mu
c     calculate beta angle
      beta = atan(tanb)
c     calculate VEV's
      v1    = 2*zm*sct*cos(beta)/e
      v2    = 2*zm*sct*sin(beta)/e
c     calculate Yukawa couplings if fermion masses initialized
      if (fermion) call init_tree_yukawa
c     calculate Higgs mass parameters in the Lagrangian
      sg2   = pm*pm
      sum   = sg2 - 2*abs(mu*mu)
      diff  = - (sg2 + zm2)*cos(2*beta)
      hm1   = (sum + diff)/2
      hm2   = (sum - diff)/2
      hm12  = - sg2/2*sin(2*beta)
c     Higgs mass digonalization
      ierr = intlog(hdiag())
      higgs=.true.
      return
      end

      subroutine init_ino_sector(m1,m2,m3,mu,tanb,ierr)
c     Tree-level initialization of SUSY fermion sector
c     ierr=1,2  Low -ino mass <= M_Z/2
      implicit double precision (a-h,o-z)
      double complex mu,m1,m2
      double precision m3
      double complex gm2,gm3,hmu
      logical cdiag,ndiag
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
c     gm1 = M3 (gluino mass), gm2 = M2, gm3 = M1
      common/gmass/gm1,gm2,gm3
      common/vev/v1,v2
      common/hpar/hm1,hm2,hm12,hmu
c     store mu in common
      hmu = mu
c     calculate beta angle
      beta = atan(tanb)
c     calculate VEV's
      v1    = 2*zm*sct*cos(beta)/e
      v2    = 2*zm*sct*sin(beta)/e
c     store gaugino masses in common
      gm1 = m3                  ! SU(3) mass
      gm2 = m2                  ! SU(2) mass
c     if M1=0, then M1 becomes GUT related to M2
      if (abs(m1).ne.0.d0) then
         gm3 = m1                ! U(1) mass
      else
         gm3 = 5/3.d0*st2/ct2*m2 ! U(1) mass GUT related to SU(2) mass
      end if
c     calculate -ino masses
      ierr =  intlog(cdiag()) + intlog(ndiag())
      return
      end

      subroutine init_slepton_sector(all,arr,slmi_l,slmi_r,slmi_lr,
     $     slmi_lrp,ierr)
c     Tree-level initialization of slepton sector
c     Routine assumes A_l'=0
c     slmi_l, slmi_r, slmi_lr contain slepton LL, RR and LR mass insertions
c     Exit code ierr=1: one or more slepton mass squares below 0
c     Higgs and fermion initialization has to be called prior to this routine! 
      implicit double precision (a-h,o-z)
      dimension all(3),arr(3)
      double complex slmi_l(3),slmi_r(3),slmi_lr(3,3),slmi_lrp(3,3)
      double complex ls,ks,ds,es,us,ws
      double complex lms,rms,ums,dms,qms
      logical sldiag,higgs,fermion
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/msoft/lms(3,3),rms(3,3),ums(3,3),dms(3,3),qms(3,3)
      common/required_init/higgs,fermion
      common/sf_cont/eps,indx(3,3),iconv
c     Higgs and fermion initialization has to be called prior to this one:
      if (.not.(higgs.and.fermion)) 
     $     stop 'Initialize Higgs/fermion data before slepton sector'
      do i=1,3
c     set diagonal soft slepton parameters
         lms(i,i) = (1 + eps*i)*all(i)*abs(all(i))
         rms(i,i) = (1 - eps*i)*arr(i)*abs(arr(i))
      end do
c     expand slepton LL and RR mass insertions
      do i=1,2
         do j=i+1,3
            lms(i,j) = slmi_l(indx(i,j))*sqrt(abs(lms(i,i)*lms(j,j)))
            lms(j,i) = dconjg(lms(i,j))
            rms(i,j) = slmi_r(indx(i,j))*sqrt(abs(rms(i,i)*rms(j,j)))
            rms(j,i) = dconjg(rms(i,j))
         end do
      end do
c     expand slepton LR mass insertions
      do i=1,3
         do j=1,3
            ls(i,j) = slmi_lr(i,j)*abs(lms(i,i)*rms(j,j))**0.25d0
            ks(i,j) = slmi_lrp(i,j)*abs(lms(i,i)*rms(j,j))**0.25d0
         end do
      end do
c     if slepton input data in SLHA format, rewrite them to
c     hep-ph/9511250 convention
      if (iconv.eq.1) call sl_slha_to_jr

c      call cr_mat_print(ls,3,6)
c      call ci_mat_print(ls,3,6)
c      stop

c     call diagonalization routine
      ierr = intlog(sldiag())
      return
      end

      subroutine init_squark_sector(asq,asu,asd,sqmi_l,sumi_r,sdmi_r,
     $     sumi_lr,sdmi_lr,sumi_lrp,sdmi_lrp,ierr)
c     Tree-level initialization of squark sector
c     Routine assumes A_d'=A_u'=0
c     sqmi_l, sdmi_r, sumi_r, sdmi_lr, sumi_lr contain squark LL, RR and
c     LR mass insertions
c     Exit code ierr=1: one or more squark mass squares below 0
c     Higgs and fermion initialization has to be called prior to this routine! 
      implicit double precision (a-h,o-z)
      dimension asq(3),asu(3),asd(3)
      double complex sqmi_l(3),sumi_r(3),sdmi_r(3)
      double complex sumi_lr(3,3),sdmi_lr(3,3)
      double complex sumi_lrp(3,3),sdmi_lrp(3,3)
      double complex ls,ks,ds,es,us,ws
      double complex lms,rms,ums,dms,qms
      logical sqdiag,higgs,fermion
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/msoft/lms(3,3),rms(3,3),ums(3,3),dms(3,3),qms(3,3)
      common/required_init/higgs,fermion
      common/sf_cont/eps,indx(3,3),iconv
c     Higgs and fermion initialization has to be called prior to this one:
      if (.not.(higgs.and.fermion)) 
     $     stop 'Initialize Higgs/fermion data before squark sector'
      do i=1,3
c     set diagonal soft squark parameters
         ums(i,i) = (1 - eps*i)*asu(i)*abs(asu(i))
         dms(i,i) = (1 - eps*i)*asd(i)*abs(asd(i))
         qms(i,i) = (1 + eps*i)*asq(i)*abs(asq(i))
      end do
c     expand LL and RR squark mass insertions
      do i=1,2
         do j=i+1,3
            qms(i,j) = sqmi_l(indx(i,j))*sqrt(qms(i,i)*qms(j,j))
            qms(j,i) = dconjg(qms(i,j))
            dms(i,j) = sdmi_r(indx(i,j))*sqrt(dms(i,i)*dms(j,j))
            dms(j,i) = dconjg(dms(i,j))
            ums(i,j) = sumi_r(indx(i,j))*sqrt(ums(i,i)*ums(j,j))
            ums(j,i) = dconjg(ums(i,j))
         end do
      end do
c     expand squark LR mass insertions
      do i=1,3
         do j=1,3
            ds(i,j) = sdmi_lr(i,j)*abs(qms(i,i)*dms(j,j))**0.25d0
            us(i,j) = sumi_lr(i,j)*abs(qms(i,i)*ums(j,j))**0.25d0
            es(i,j) = sdmi_lrp(i,j)*abs(qms(i,i)*dms(j,j))**0.25d0
            ws(i,j) = sumi_lrp(i,j)*abs(qms(i,i)*ums(j,j))**0.25d0
         end do
      end do

c     if squark input data in SLHA format, rewrite them to
c     hep-ph/9511250 convention
      if (iconv.eq.1) call sq_slha_to_jr

c     call diagonalization routine
      ierr = intlog(sqdiag())
      return
      end

      subroutine reset_phys_data()
c     After parameter changes some quantities has to be recalculated
      implicit double precision (a-h,o-z)
      logical zstat,hstat
      common/zwidth/z_gam,zfact,zstat
      common/hm_stat/hstat
c     Reset Z0 width calculations
      zstat = .false.
c     Physical Higgs masses have to be recalculated after parameter change
      hstat = .false.
      return
      end

      integer function intlog(x)
      logical x
      if (x) then 
         intlog = 1
      else
         intlog = 0 
      end if
      return
      end
 
      subroutine vpar_update(zm_new,wm_new,alpha_new,st2_new)
c     Consistent calculation of data related to M_Z, M_W
      implicit double precision (a-h,o-z)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/crdat/pbarn,ae,ve
      common/fvert/qf(4),vf(4),af(4),nc(4)
      common/stand/hm,vev
      external init_phys,init_const
      alpha = alpha_new
      e2 = 4*pi*alpha
      e = sqrt(e2)
      zm = zm_new
      wm = wm_new
c      ct = wm/zm                ! on shell value of c_W
c      ct2 = ct*ct
c      st2 = 1 - ct2
      st2 = st2_new             ! external s_W^2, usually MSbar
      ct2 = 1 - st2
      ct = sqrt(ct2)       
      st = sqrt(st2)
      sct = st*ct
      sct2 = sct*sct
      wm2 = wm*wm
      zm2 = zm*zm
c     calculate vev in SM
      vev = 2*wm*st/e
c     fermion coupling actualization for (v,l,u,d)
      iso = 1
      do 10 ig=1,4
        af(ig) = iso/4.d0/sct
        vf(ig) = af(ig)*(1 - 4*iso*qf(ig)*st2)
10      iso = - iso
      ae = af(2)
      ve = vf(2)
      return
      end

      subroutine ckm_update(v1,v2,v3,phi)
c     Kobayashi-Maskawa matrix initialization
      implicit double precision (a-h,o-z)
      double complex v_ckm,ckm_herm,cphi
      common/km_par/vkm1,vkm2,vkm3,phi0
      common/km_mat/ckm_herm(3,3)
      common/ckm/v_ckm(3,3)
      vkm1 = v1
      vkm2 = v2
      vkm3 = v3
      phi0 = phi
      sith = v1/sqrt(v1*v1 + v2*v2)
      coth = v2/sqrt(v1*v1 + v2*v2)
      cobe = sqrt(v1*v1 + v2*v2)
      sibe = sqrt(1 - v1*v1 - v2*v2)
      siga = v3/cobe
      coga = sqrt(1 - siga*siga)
      cphi = exp(dcmplx(0.d0,phi))
c     unitary ckm initialization
      v_ckm(1,1) = cobe*coth
      v_ckm(1,2) = cobe*sith
      v_ckm(1,3) = sibe/cphi
      v_ckm(2,1) = - siga*coth*sibe*cphi - sith*coga
      v_ckm(2,2) = coga*coth - siga*sibe*sith*cphi
      v_ckm(2,3) = siga*cobe
      v_ckm(3,1) = - sibe*coga*coth*cphi + siga*sith
      v_ckm(3,2) = - coga*sibe*sith*cphi - siga*coth
      v_ckm(3,3) = coga*cobe
      do i=1,3
         do j=1,3
            ckm_herm(i,j) = dconjg(v_ckm(j,i))
         end do
      end do
      call ckm_eff_init
      return
      end

      subroutine ckm_init(s12,s23,s13,phi)
      implicit double precision (a-h,o-z)
      double complex v_ckm,ckm_herm,cphi
      common/km_par/vkm1,vkm2,vkm3,phi0
      common/km_mat/ckm_herm(3,3)
      common/ckm/v_ckm(3,3)
      sith = s12
      coth = sqrt(1 - s12*s12)
      sibe = s13
      cobe = sqrt(1 - s13*s13)
      siga = s23      
      coga = sqrt(1 - s23*s23)
      vkm1 = sith*cobe
      vkm2 = coth*cobe
      vkm3 = siga*cobe
      phi0 = phi
      cphi = exp(dcmplx(0.d0,phi))
c     unitary ckm initialization
      v_ckm(1,1) = cobe*coth
      v_ckm(1,2) = cobe*sith
      v_ckm(1,3) = sibe/cphi
      v_ckm(2,1) = - siga*coth*sibe*cphi - sith*coga
      v_ckm(2,2) = coga*coth - siga*sibe*sith*cphi
      v_ckm(2,3) = siga*cobe
      v_ckm(3,1) = - sibe*coga*coth*cphi + siga*sith
      v_ckm(3,2) = - coga*sibe*sith*cphi - siga*coth
      v_ckm(3,3) = coga*cobe
      call correct_unitarity(v_ckm,3)
      do i=1,3
         do j=1,3
            ckm_herm(i,j) = dconjg(v_ckm(j,i))
         end do
      end do
      call ckm_eff_init
      return
      end

      subroutine ckm_wolf(al,a,rhobar,etabar)
c     KM matrix initialization in the (approximate) Wolfenstein
c     parameterization 
      implicit double precision (a-h,o-z)
      double complex tmp
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/kpivv/ak0,del_ak0,akp,del_akp,pc,del_pc,alam
      al2 = al*al
      al3 = al2*al
      al4 = al3*al
      tmp = dcmplx(rhobar,etabar)
      tmp = a*al3*sqrt((1 - a*a*al4)/(1 - al2))*tmp/(1 - a*a*al4*tmp)
c     CKM sin(theta_ij) and phase delta
      s12 = al
      s23 = a*al2
      s13 = abs(tmp)
      if (dble(tmp).eq.0) then
         del = pi/2
      else
         del = atan(dimag(tmp)/dble(tmp))
      end if
      if (dble(tmp).lt.0) del = pi + del
      call ckm_init(s12,s23,s13,del)
c     store Lambda for K->pivv calculations
      alam = al
      return
      end

      subroutine ckm_eff_init
      implicit double precision (a-h,o-z)
      double complex ckm,ckm_phys,ckm0,udl,udr,uul,uur
      common/ckm/ckm(3,3)
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      do i=1,3
         do j=1,3
            ckm_phys(i,j) = ckm(i,j)
            ckm0(i,j)     = ckm(i,j)
         end do
      end do
      return
      end

      subroutine vert_stat(vstatus,fstatus)
c     If vstatus=.false. vertex formfactors are not calculated and
c     equal to zero - significantly speeds up calculations
      logical vstat,vstatus,fstat,fstatus
      common/vswitch/vstat,fstat
      vstat = vstatus
      fstat = fstatus
      return
      end

      subroutine set_active_sector(jh,jc,jn,jg)
c     Integers ih,ic,in,ig switches on/off contributions of
c     (respectively) gauge+Higgs bosons, charginos, neutralinos and
c     gluinos to the Green's functions.  
      common/debug_4q/ih,ic,in,ig
      external init_control
      if ((jh.eq.0).or.(jh.eq.1)) then 
         ih = jh
      else
         stop 'incorrect IH in set_active_sector'
      end if
      if ((jc.eq.0).or.(jc.eq.1)) then 
         ic = jc
      else 
         stop 'incorrect IC in set_active_sector'
      end if
      if ((jn.eq.0).or.(jn.eq.1)) then 
         in = jn
      else 
         stop 'incorrect IN in set_active_sector'
      end if
      if ((jg.eq.0).or.(jg.eq.1)) then 
         ig = jg
      else 
         stop 'incorrect IG in set_active_sector'
      end if
      return
      end

      subroutine sl_slha_to_jr
c     translation of slepton input data from SLHA2 to hep-ph/9511250
c     format
      implicit double precision (a-h,o-z)
      double complex ls,ks,ds,es,us,ws
      double complex lms,rms,ums,dms,qms
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/msoft/lms(3,3),rms(3,3),ums(3,3),dms(3,3),qms(3,3)
      call transpose_cmatrix(ls,3)
      call transpose_cmatrix(rms,3)
      return
      end

      subroutine sq_slha_to_jr
c     translation of squark input data from SLHA2 to hep-ph/9511250
c     format
      implicit double precision (a-h,o-z)
      double complex ls,ks,ds,es,us,ws
      double complex lms,rms,ums,dms,qms
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/msoft/lms(3,3),rms(3,3),ums(3,3),dms(3,3),qms(3,3)
      call transpose_cmatrix(ds,3)
      call transpose_cmatrix(us,3)
      call transpose_cmatrix(dms,3)
      call transpose_cmatrix(ums,3)
      do i=1,3
         do j=1,3
            us(i,j) = - us(i,j)
      end do
      end do
      return
      end

      subroutine transpose_cmatrix(x,n)
c     transpose double complex matrix of size (n x n)
      implicit double precision (a-h,o-z)
      double complex x(n,n),tmp
      do i=1,n-1
         do j=i+1,n
            tmp = x(i,j)
            x(i,j) = x(j,i)
            x(j,i) = tmp
         end do
      end do      
      return
      end

      block data init_phys
c     Physical quantities initialization
      implicit double precision (a-h,o-z)
      double complex ckm
      double complex lms,rms,ums,dms,qms
      double complex ls,ks,ds,es,us,ws
      double complex gm2,gm3
      logical zstat
      common/zwidth/z_gam,zfact,zstat
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/msoft/lms(3,3),rms(3,3),ums(3,3),dms(3,3),qms(3,3)
      common/gmass/gm1,gm2,gm3
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/km_mat/ckm(3,3)
      common/km_par/vkm1,vkm2,vkm3,phi
      common/fmass/em(3),um(3),dm(3)
      common/qmass_pole/ump(3),dmp(3)
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
      common/tau_gam/br_tau_evv
      common/crdat/pbarn,ae,ve
      common/fvert/qf(4),vf(4),af(4),nc(4)
      common/nc_exp/z_inv,z_vis,z_width,cr_peak
      common/stand/hm,vev
      common/fermi/g_fermi
      data g_fermi/1.16639d-5/
      data nc/1,1,3,3/
      data qf/0.d0,-1.d0,0.66666666666666d0,-0.33333333333333d0/
      data af/0.5988406186d0,-0.5988406186d0,
     1        0.5988406186d0,-0.5988406186d0/
      data vf/0.5988406186d0,-0.0602802702d0,
     1        0.2398003864d0,-0.4193205025d0/
      data pbarn,ae,ve/3.8937966d8,-0.605059d0,-0.0764d0/
      data z_gam,zfact,zstat/2.496d0,1.d0,.false./
c     z_width, cr_peak are HADRONIC width and cross section
      data z_inv,z_vis,z_width,cr_peak/7.6d-3,23.d-3,1.746d0,41.51d3/
c     Input parameters for quark mass calculations
c     uml contains u-quark masses at scales amuu, similarly for d-quarks
c     umu, dmu are quark masses at mt=umu(3),
c     to be calculated in init_run_qmass
c     Light fermion masses at mu = 2GeV used by Ciuchini et al.
      data umu,uml,amuu/2.1d-3,1.279d0,163.1d0,
     $                  2.1d-3,1.279d0,163.1d0,
     $                  2.d0,   1.279d0,163.1d0/
      data dmu,dml,amud/4.7d-3,9.34d-2,4.19d0,
     $                  4.7d-3,9.34d-2,4.19d0,
     $                  2.d0,  2.d0,   4.19d0/
      data em/5.109989d-4,1.05658d-1,1.77684d0/
      data um/2.1d-3,1.279d0,163.1d0/
      data dm/4.73d-3,9.34d-2,4.18d0/
      data ump/4d-3,1.5d0,173.5d0/
      data dmp/7d-3,0.110d0,4.77d0/
      data vkm1,vkm2,vkm3,phi/0.222d0,0.975d0,0.044d0,0.d0/
      data ckm/(1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0),
     $         (0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0),
     $         (0.d0,0.d0),(0.d0,0.d0),(1.d0,0.d0)/
      data zm,zm2,wm,wm2/9.11876D+01,8.315251D+03,8.0398D+01,6.4963D+03/
      data st2,ct2/2.18433072D-01,7.81566928D-01/
      data st,ct/4.67368240D-01,8.84062740D-01/
      data sct,sct2/4.13182847D-01,1.70720065D-01/
      data e2,e/9.1701236276D-02,3.0282211986D-01/
      data alpha/7.2973525205D-03/
      data sq2,pi/1.414213562373095d0,3.1415926536d0/
      data hm,vev/100.d0,248.663d0/
      data br_tau_evv/0.1782d0/
      data lms/9*(0.d0,0.d0)/,rms/9*(0.d0,0.d0)/
      data dms/9*(0.d0,0.d0)/,ums/9*(0.d0,0.d0)/,qms/9*(0.d0,0.d0)/
      data ls/9*(0.d0,0.d0)/,ks/9*(0.d0,0.d0)/
      data ds/9*(0.d0,0.d0)/,es/9*(0.d0,0.d0)/
      data us/9*(0.d0,0.d0)/,ws/9*(0.d0,0.d0)/
      data gm1,gm2,gm3/600.d0,(200.d0,0.d0),(100.d0,0.d0)/
      end

      block data init_4q
      implicit double precision (a-h,o-z)
      logical init_eta,init_alpha_susy,eff_yuk,init_eff_yuk
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sm_4q/eta_cc,eta_ct,eta_tt,eta_b,bk_sm,bd_sm,bb_sm(2)
      common/ev_mat_4q/vxx(4),sxx(4,2,2),sxy(4,2,2),init_eta
      common/meson_data/dmk,amk,epsk,fk,dmd,amd,fd,
     $    amb(2),dmb(2),tau_b(2),fb(2)
      common/dtau_data/dmbp,rd,del_rd,rds,del_rds
      common/bx_4q/bk(5),bd(5),bb(2,5),amu_k,amu_d,amu_b
      common/eff_yuk_4q/eps_d(3),eff_yuk,init_eff_yuk
      common/diag_type/ihpeng,izpeng,ifpeng,ibox
c     include Higgs penguins
      data ihpeng,izpeng,ifpeng,ibox/4*1/
c     inclusion of eps,eps' terms in 4q box calculations
      data eff_yuk,init_eff_yuk/.false.,.true./
      data eps_d/3*0.d0/
c     K mesons
      data dmk,amk,fk,epsk/3.483d-15,0.497614d0,0.156d0,2.228d-3/
c     B_K ordering: VLL, SLL(2), LR(2)
      data bk/0.61d0,0.76d0,0.51d0,0.96d0,1.30d0/
      data bk_sm/0.724d0/
c     D mesons
      data dmd,amd,fd/1.56d-14,1.8645d0,0.2d0/
      data bd/5*1.d0/
      data bd_sm/1.d0/
c     B_d,B_s mesons
      data amb/5.2792d0,5.3668d0/
      data dmb/3.337d-13,1.17d-11/
      data tau_b/1.519d-12,1.516d-12/
      data fb/0.192d0,0.228d0/
      data bb/2*0.87d0,2*0.8d0,2*0.71d0,2*1.71d0,2*1.16d0/
      data bb_sm/1.305d0,1.345d0/
c     Scales for B_X calculations
      data amu_k,amu_d,amu_b/2.d0,2.d0,4.6d0/
c     Debug/status variablesq
      data init_eta/.true./
      data init_alpha_susy/.true./
c     For SM-like formfactor (vector-LL) in kaon and (partially) BB
c     mixing much more refined calculation of QCD factors exist.
c     Lets replace our simplified formulae by the exact numerical values,
c     taken from arXiv: 0909.1333 [hep-ph]
c     KK mixing:
      data eta_cc,eta_ct,eta_tt/1.87d0,0.496d0,0.5765d0/
c     BB mixing:
      data eta_b/0.55d0/
c     B_u meson and  Br(B->D(D^*) tau nu)/Br(B->D(D^*) l nu) in SM
      data dmbp/5.27917d0/      ! M_B_u
      data rd,del_rd/0.297d0,0.017d0/ 
      data rds,del_rds/0.252d0,0.003d0/
      end


      block data init_const
c     common numerical constants initialization
      implicit double precision (a-h,o-z)
      double complex cz,co,ci
      common/delta/delta(6,6)
      common/eps/eps(2,2)
      common/num/cz,co,ci,zero,one
      data cz,co,ci/(0.d0,0.d0),(1.d0,0.d0),(0.d0,1.d0)/
      data zero,one/0.d0,1.d0/
      data delta/1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
     $           0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,
     $           0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,
     $           0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,
     $           0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,
     $           0.d0,0.d0,0.d0,0.d0,0.d0,1.d0/
      data eps/0.d0,-1.d0,
     $         1.d0,0.d0/
      end

      block data init_control
c     control variables initialization
      implicit double precision (a-h,o-z)
      logical vstat,fstat,hstat,bstat
      logical higgs,fermion
      common/required_init/higgs,fermion
      common/vswitch/vstat,fstat
      common/bswitch/bstat
      common/hm_stat/hstat
      common/dimreg/idflag
      common/debug_4q/ih,ic,in,ig
      common/sf_cont/eps,indx(3,3),iconv
      common/mssm_charged_higgs_min_index/mhmin
c     status of Higgs and fermion sector initialization
      data higgs,fermion/2*.false./
      data idflag/1/
      data ih,ig,ic,in/4*1/
      data vstat,fstat,hstat,bstat/2*.true.,2*.false./
      data eps/1.d-5/
      data indx/0,1,3,
     $          1,0,2,
     $          3,2,0/
      data iconv/2/
      data mhmin/1/
      end

      block data init_units
      implicit double precision (a-h,o-z)
      common/ph_units/hbar,gev_cm,gev_s
c     1 GeV in cm^(-1)
      data gev_cm/5.067689d+13/
c     1 GeV in s^(-1)
      data gev_s/1.519255d+24/
c     Planck constant in GeV sec 
      data hbar/6.58211915d-25/
      end

