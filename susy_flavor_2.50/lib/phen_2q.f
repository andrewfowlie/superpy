c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: PHEN_2Q.F
c     Released: 29:08:2007 (J.R.)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expressions for QCD evolution of          c
c     Wilson coefficients of the operators present in 2-d quark    c
c     mixing and for the phenomenological quantities like B->ll,   c
c     B->Xll, K/B-> pivv (to be added gradually).                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine eta_2q_evol
c     QCD evolution of Wilson coefficients of the effective
c     2-quark operators
      implicit double precision (a-h,o-z)
      logical init_eta_2q,init_alpha_susy
      dimension qm(3)
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
c     put into the common below evolution coefficients...
      common/ev_mat_2q/vx(2),sx(2),tx(2),init_eta_2q
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
      common/bx_4q/bk(5),bd(5),bb(2,5),amu_k,amu_d,amu_b
c     running s,b masses at m_t scale
      call init_run_qmass
c     NLO alpha_s at M_SUSY scale(s)
c     M_S_{U,D} = (M_gluino + average M_{U,D} mass)/2
      if (init_alpha_susy) call init_alpha_s_susy
c     running quark masses (c,b,t)
      qm(1) = uml(2)
      qm(2) = dml(3)
      qm(3) = umu(3)
c     evolution from M_SUSY to mu_t
c     start from M_SUSY = (M_gluino + averaged M_D)/2
      al_s = alfas_nlo(qm(3))/4/pi
      eta = g3d*g3d/al_s/16/pi/pi
c     define evolution below for V, S and T sectors:
      vx(2) = 1
      sx(2) = 1
      tx(2) = 1
c     evolution from mu_t to mu_B
      al_s = alfas_nlo(amu_b)/4/pi
      eta = alfas_nlo(qm(3))/alfas_nlo(amu_b)
c     define evolution below for V, S and T sectors:
      vx(1) = 1
      sx(1) = 1
      tx(1) = 1
c     initialization of QCD factors finished:
      init_eta_2q = .false.
      return
      end

      subroutine dl_wil_run(i,j,k,l)
      implicit double precision (a-h,o-z)
      logical init_eta_2q
      double complex dl_vll,dl_vrr,dl_vlr,dl_vrl,dl_sll,dl_srr,
     $     dl_slr,dl_srl,dl_tl,dl_tr
      double complex dls_vll,dls_vrr,dls_vlr,dls_vrl,dls_sll,dls_srr,
     $     dls_slr,dls_srl,dls_tl,dls_tr
      common/dl_wil_coeff/dls_vll,dls_vrr,dls_vlr,dls_vrl,
     $     dls_sll,dls_srr,dls_slr,dls_srl,dls_tl,dls_tr
      common/ev_mat_2q/vx(2),sx(2),tx(2),init_eta_2q
      if (init_eta_2q) call eta_2q_evol()
c     NLO evolution of the MSSM part
c     First step: from mu = M_SUSY = (M_gluino + aver M_D)/2 to mu = m_t
      dls_vll  = vx(2)*dl_vll(i,j,k,l)
      dls_vrr  = vx(2)*dl_vrr(i,j,k,l)
      dls_vlr  = vx(2)*dl_vlr(i,j,k,l)
      dls_vrl  = vx(2)*dl_vrl(i,j,k,l)
      dls_sll  = sx(2)*dl_sll(i,j,k,l)
      dls_srr  = sx(2)*dl_srr(i,j,k,l)
      dls_slr  = sx(2)*dl_slr(i,j,k,l)
      dls_srl  = sx(2)*dl_srl(i,j,k,l)
      dls_tl   = tx(2)*dl_tl(i,j,k,l)
      dls_tr   = tx(2)*dl_tr(i,j,k,l)
c     Second step: from mu = m_t to mu = m_b
      dls_vll  = vx(1)*dls_vll
      dls_vrr  = vx(1)*dls_vrr
      dls_vlr  = vx(1)*dls_vlr
      dls_vrl  = vx(1)*dls_vrl
      dls_sll  = sx(1)*dls_sll
      dls_srr  = sx(1)*dls_srr
      dls_slr  = sx(1)*dls_slr
      dls_srl  = sx(1)*dls_srl
      dls_tl   = tx(1)*dls_tl
      dls_tr   = tx(1)*dls_tr
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Br(B->l^+l^-) calculation                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function b_ll(i,j,k,l)
c     Decay of B_d(s)-> l_K l_L pair
c     i=3 denotes B decay, j=3 is \bar B decay
c     second index j,i=1 or 2 defines B meson, B_d or B_s respectively
      implicit double precision (a-h,o-z)
      double complex fs,fp,fv,fa
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/meson_data/dmk,amk,epsk,fk,dmd,amd,fd,
     $    amb(2),dmb(2),tau_b(2),fb(2)
      common/fmass/em(3),um(3),dm(3)
      double complex dls_vll,dls_vrr,dls_vlr,dls_vrl,dls_sll,dls_srr,
     $     dls_slr,dls_srl,dls_tl,dls_tr
      common/dl_wil_coeff/dls_vll,dls_vrr,dls_vlr,dls_vrl,
     $     dls_sll,dls_srr,dls_slr,dls_srl,dls_tl,dls_tr
      common/ph_units/hbar,gev_cm,gev_s
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
      common/bx_4q/bk(5),bd(5),bb(2,5),amu_k,amu_d,amu_b
      external init_4q,init_2q,init_units
      if (max(i,j).ne.3) stop 'No b quark index in b_ll?'
      if (min(i,j).eq.3) stop 'B_b decay not implemented!'
      call dl_wil_run(i,j,k,l)
c     form factors
      ii = min(i,j)      
      amb2 = amb(ii)*amb(ii)
      rm = amb2/(qmass_nlo(dml(ii),amud(ii),amu_b)
     $     + qmass_nlo(dml(3),amud(3),amu_b))
      del = em(k) - em(l)
      sum = em(k) + em(l)
      fv = dls_vrr - dls_vll + dls_vrl - dls_vlr 
      fa = dls_vll + dls_vrr - dls_vlr - dls_vrl
      fs = rm*(dls_sll - dls_srr + dls_slr - dls_srl)
      fp = rm*(dls_slr + dls_srl - dls_sll - dls_srr)
c     matrix element
      b_ll = (amb2 - sum*sum)*abs(fs)**2 + (amb2 - del*del)*abs(fp)**2
     $     + sum*sum*(amb2 - del*del)*abs(fa)**2
     $     + 2*sum*(amb2 - del*del)*dble(fp*dconjg(fa))
      if (k.ne.l) b_ll = b_ll + del*del*(amb2 - sum*sum)*abs(fv)**2
     $     - 2*del*(amb2 - sum*sum)*dble(fs*dconjg(fv))
c     branching ratio itself
      b_ll = b_ll*fb(ii)*fb(ii)*tau_b(ii)/128/pi/amb(ii)/hbar
     $     * sqrt(1 - sum*sum/amb2)*sqrt(1 - del*del/amb2)
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Br(K_L -> pi^0 \bar v v) and  Br(K^+ -> pi^+ \bar v v) calculation  c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine k_pivv(br_k0,br_kp)
c     Decays K^0_L -> pi^0 \bar v v and K^+ -> pi^+ \bar v v
c     compare hep-ph/0408142
      implicit double precision (a-h,o-z)
      double complex xx
      double complex dd_vv_l,dd_vv_r
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/kpivv/ak0,del_ak0,akp,del_akp,pc,del_pc,alam
      external init_2q
      xx = (0.d0,0.d0)
      do k=1,3
         do l=1,3
            xx = xx + dd_vv_l(2,1,k,l) + dd_vv_r(2,1,k,l)
         end do 
      end do
      xx = - (4*pi*st2*wm/e2)**2*xx/3/alam**5
c     Br(K^0_L -> pi^0 \bar v v)
      br_k0 = ak0*dimag(xx)**2
c     Br(K^+ -> pi^+ \bar v v)
      br_kp = akp*(dimag(xx)**2 
     $     + dble(ckm_phys(2,2)*dconjg(ckm_phys(2,1))*pc/alam + xx)**2)
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Neutron EDM calculation                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function edm_n()
c     complete neutron EDM. Included:
c       - electric dipole moment of quarks
c       - chromoelectric dipole moment of quarks
c       - chromoelectric dipole moment of gluon
c       - approximate QCD coefficients
      implicit double precision (a-h,o-z)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/edm_qcd/eta_ed,eta_eu,eta_cd,eta_cu,eta_g,alamx
      common/ph_units/hbar,gev_cm,gev_s
      external init_units
      external init_2q
      edm_n = eta_ed*edm_d(1) + eta_eu*edm_u(1) 
     $     + e*(eta_cd*cdm_d(1) + eta_cu*cdm_u(1))
     $     + e/4.d0/pi*eta_g*alamx*cdm_g()
      edm_n = edm_n/e/gev_cm
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Br(B^+ -> tau^+ v) and  Br(B^0 -> D tau v)/Br(B^0 -> Dlv)           c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine b_taunu(br_taunu,dtaunu_rat,dstaunu_rat)
c     Decays B_u -> tau v and B_u -> D (D^*) tau v (charged B_u meson)
      implicit double precision (a-h,o-z)
      double complex c_l,c_r
      double complex yh_eff_r,yhl_eff_r
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/meson_data/dmk,amk,epsk,fk,dmd,amd,fd,
     $    amb(2),dmb(2),tau_b(2),fb(2)
      common/dtau_data/bm,rd,del_rd,rds,del_rds
      common/fmass/em(3),um(3),dm(3)
      common/fermi/g_fermi
      common/ph_units/hbar,gev_cm,gev_s
      external init_4q
c     New Physics contributions to D->tau nu
      c_l = - yh_eff_l(3,1,1)*dconjg(yhl_eff_r(3,3,1))/sq2/2/cm(1)
     $     /cm(1)/g_fermi/ckm_phys(1,3)
      c_r = - yh_eff_r(3,1,1)*dconjg(yhl_eff_r(3,3,1))/sq2/2/cm(1)
     $     /cm(1)/g_fermi/ckm_phys(1,3)
      br_taunu = tau_b(1)*gev_s*bm/8/pi*abs(g_fermi*em(3)*fb(1)
     $     * ckm_phys(1,3)*(1 - em(3)*em(3)/bm/bm)
     $     * (1 + bm*bm/dm(3)/em(3)*(c_r - c_l)))**2
c     New Physics contributions to D->tau nu
      c_l = - yh_eff_l(3,2,1)*dconjg(yhl_eff_r(3,3,1))/sq2/2/cm(1)
     $     /cm(1)/g_fermi/ckm_phys(2,3)
      c_r = - yh_eff_r(3,2,1)*dconjg(yhl_eff_r(3,3,1))/sq2/2/cm(1)
     $     /cm(1)/g_fermi/ckm_phys(2,3)
      dtaunu_rat = rd*(1 + 1.5d0*dble(c_r + c_l) + abs(c_r + c_l)**2)
      dstaunu_rat = rds*(1 + 0.12d0*dble(c_r - c_l) 
     $     + 0.05*abs(c_r - c_l)**2)
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Br(t->ch(dh))                                                        c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function gam_suu(i,k)
c     Decays width u^J -> u^K H^0_i
      implicit double precision (a-h,o-z)
      double complex uus(2),uug(2),form(2)
      double complex zu,zd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/hmass_EPA/pm,hm1,hm2,sa,ca,sb,cb
      common/qmass_pole/um(3),dm(3)
      common/fmass/emr(3),umr(3),dmr(3)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/debug_4q/ih,ic,in,ig 
      tm = um(3)                ! top mass pole
      qm = um(k)                ! lighter quark mass
      if (i.eq.1) then          ! Higgs mass 
         hm = hm1
      else
         hm = hm2
      end if
      if (tm.le.(qm + hm)) then
         gam_suu = 0.d0
         return
      end if
c     store control variables
      ih0 = ih
      ic0 = ic
      in0 = in
      ig0 = ig
c     W+Higgs contribution to h/H-uu vertex (at MA)
      ic = 0
      in = 0
      ig = 0
      call suu_vert(i,3,k,form)
c     SUSY contribution to h/H-uu vertex (at MSUSY)
      ih = 0
      ic = ic0                  ! ic0...ig0 at initial values
      in = in0
      ig = ig0
      call suu_vert(i,3,k,uus)
      ih = ih0                  ! restore also ih0
c     gluon-uu vertex (at MSUSY, only for i=2 - lighter CP-even Higgs)
      do ii =1,2
         uug(ii) = (0.d0,0.d0)
      end do
      if (i.eq.2) call gluu_vert(3,k,uug)
c     evaluate SUSY scale as average squark mass scale
      susy = 0.d0
      do ii=1,6
         susy = susy + sum(ii) + sdm(ii)
      end do
      susy = susy/12            ! average SUSY scale
      tmr = umr(3)              ! running mt(mt)
      r = alfas(susy)/alfas(tmr) 
      b3 = 7.d0                 ! 11 - 2 N_f/3 with N_f=6
c     add all formfactors including RGE running (uug = - 2 C_g!)
      do ii =1,2
c     full C_h(mt)
         form(ii) = form(ii) + (uus(ii) + 6/7.d0*e/st/wm*tmr**2*(1 -
     $        r**(14.d0/3/b3))*uug(ii))*r**(-4/b3) 
c     full C_g(mt) (-1/2 included here)
         uug(ii) = - uug(ii)/2*r**(2.d0/3/b3) 
      end do
      f1 = (((tm - qm)**2 - hm**2)*((tm + qm)**2 - hm**2))**0.5d0
      f2 = (tm**2 + qm**2 - hm**2)/2/tm
      gam_suu_h = f1*(f2*(abs(form(1))**2 + abs(form(2))**2) 
     $     + 2*qm*dble(form(1)*dconjg(form(2))))/16/pi/tm**2
      gam_suu_g = e*tm**3/32/pi/st/wm*(1 - (hm/tm)**2)**2
     $     *dble(dconjg(form(1))*uug(2) + dconjg(form(2))*uug(1))
      gam_suu = 1.018d0*gam_suu_h + 0.049d0*gam_suu_g
      return
      end

      double precision function br_suu(i,k)
c     Br(u^J -> u^K H^0_i) normalized to t->bW
      implicit double precision (a-h,o-z)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/qmass_pole/um(3),dm(3)
      xt = um(3)*um(3)/wm2
      gam_tbw = e2/st2/64/pi*um(3)*(2 + xt)*(1 - 1/xt)**2
      gam_s   = gam_suu(i,k)
      br_suu = gam_s/(gam_tbw + gam_s)
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Data initialization block                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      block data init_2q
      implicit double precision (a-h,o-z)
      logical init_eta_2q
c     Neutron EDM QCD correction factors and chiral symmetry breaking scale
      common/edm_qcd/eta_ed,eta_eu,eta_cd,eta_cu,eta_g,alamx
      common/ev_mat_2q/vx(2),sx(2),tx(2),init_eta_2q
      common/kpivv/ak0,del_ak0,akp,del_akp,pc,del_pc,alam
      data init_eta_2q/.true./
      data ak0,akp/2.231d-10,5.173d-11/
      data del_ak0,del_akp/0.013d-10,0.024d-11/
      data pc,del_pc/0.41d0,0.03d0/
      data alam/0.225d0/
      data eta_ed,eta_eu,eta_cd,eta_cu/0.79d0,-0.2d0,0.59d0,0.3d0/
      data eta_g,alamx/3.4d0,1.18d0/
      end 
