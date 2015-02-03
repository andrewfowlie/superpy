c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: PHEN_4Q.F
c     Released: 21:12:1999 (P.Ch.)
c     Revised: 26:02:2001 (J.R.)
c     Improved calculation of Wilson coefficient evolution
c     added a la Buras et al, hep-ph/0102316

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expressions for QCD evolution of          c
c     Wilson coefficients of the operators present in 4-d quark    c
c     mixing and for the phenomenological quantities like eps_K,   c
c     Delta_m_Bs and Delta_m_Bd etc.                               c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine eta_4q_evol
c      QCD evolution of Wilson coefficients of the effective
c      4-quark operators
      implicit double precision (a-h,o-z)
      logical init_eta,init_alpha_susy
      dimension qm(3)
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/ev_mat_4q/vxx(4),sxx(4,2,2),sxy(4,2,2),init_eta
      common/bx_4q/bk(5),bd(5),bb(2,5),amu_k,amu_d,amu_b
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
c     running s,b masses at m_t scale
      call init_run_qmass
c      NLO alpha_s at M_SUSY scale(s)
c      M_S_{U,D} = (M_gluino + average M_{U,D} mass)/2
      if (init_alpha_susy) call init_alpha_s_susy
c      running quark masses (c,b,t)
      qm(1) = uml(2)
      qm(2) = dml(3)
      qm(3) = umu(3)
c      evolution from M_SUSY to mu_t
c      start from M_SUSY = (M_gluino + averaged M_D)/2
      al_s = alfas_nlo(qm(3))/4/pi
      eta = g3d*g3d/al_s/16/pi/pi
      eta1 = eta**(6/21.d0)
c      "VLL" sector
      vxx(3) = eta1*(1 + al_s*1.3707d0*(1 - eta))
c      "LR" sector
      eta1 = eta**(3/21.d0)
      eta2 = eta**(24/21.d0)
      sxy(3,1,1) = eta1 + al_s*(0.9219d0/eta2 + ( - 2.2194d0
     $    + 1.2975d0*eta)*eta1)
      sxy(3,1,2) = al_s*1.3828d0*(eta2 - 1/eta2)
      sxy(3,2,1) = 2/3.d0*(eta1 - 1/eta2)
     $    + al_s*(( - 10.1463d0 + 0.8650d0*eta)*eta1
     $    + ( - 6.4603d0 + 15.7415d0*eta)/eta2)
      sxy(3,2,2) = 1/eta2 + al_s*(0.9219d0*eta2
     $    + (9.6904d0 - 10.6122d0*eta)/eta2)
c      "SLL" sector
      eta1 = 1/eta**0.6916d0
      eta2 = eta**0.7869d0
      sxx(3,1,1) = 1.0153d0*eta1 - 0.0153d0*eta2
     $    + al_s*(eta1*(5.6478d0 - 6.0350d0*eta)
     $    + eta2*(0.3272d0 + 0.0600d0*eta))
      sxx(3,1,2) = 1.9325d0*(eta1 - eta2)
     $    + al_s*(eta1*(10.7494d0 - 37.9209d0*eta)
     $    + eta2*(41.2256d0 - 14.0841d0*eta))
      sxx(3,2,1) = - 0.0081d0*(eta1 - eta2)
     $    + al_s*(eta1*(0.0454d0 + 0.0479d0*eta)
     $    + eta2*(-0.0618d0 - 0.0315d0*eta))
      sxx(3,2,2) = - 0.0153d0*eta1 + 1.0153d0*eta2
     $    + al_s*(eta1*(0.0865d0 + 0.3007d0*eta)
     $    + eta2*(-7.7870d0 + 7.3999d0*eta))
c      repeat for M_SUSY = (M_gluino + averaged M_U)/2
      al_s = alfas_nlo(qm(3))/4/pi
      eta = g3u*g3u/al_s/16/pi/pi
      eta1 = eta**(6/21.d0)
c      "VLL" sector
      vxx(4) = eta1*(1 + al_s*1.3707d0*(1 - eta))
c      "LR" sector
      eta1 = eta**(3/21d0)
      eta2 = eta**(24/21d0)
      sxy(4,1,1) = eta1 + al_s*(0.9219d0/eta2 + ( - 2.2194d0
     $    + 1.2975d0*eta)*eta1)
      sxy(4,1,2) = al_s*1.3828d0*(eta2 - 1/eta2)
      sxy(4,2,1) = 2/3.d0*(eta1 - 1/eta2)
     $    + al_s*(( - 10.1463d0 + 0.8650d0*eta)*eta1
     $    + ( - 6.4603d0 + 15.7415d0*eta)/eta2)
      sxy(4,2,2) = 1/eta2 + al_s*(0.9219d0*eta2
     $    + (9.6904d0 - 10.6122d0*eta)/eta2)
c      "SLL" sector
      eta1 = 1/eta**0.6916d0
      eta2 = eta**0.7869d0
      sxx(4,1,1) = 1.0153d0*eta1 - 0.0153d0*eta2
     $    + al_s*(eta1*(5.6478d0 - 6.0350d0*eta)
     $    + eta2*(0.3272d0 + 0.0600d0*eta))
      sxx(4,1,2) = 1.9325d0*(eta1 - eta2)
     $    + al_s*(eta1*(10.7494d0 - 37.9209d0*eta)
     $    + eta2*(41.2256d0 - 14.0841d0*eta))
      sxx(4,2,1) = - 0.0081d0*(eta1 - eta2)
     $    + al_s*(eta1*(0.0454d0 + 0.0479d0*eta)
     $    + eta2*(-0.0618d0 - 0.0315d0*eta))
      sxx(4,2,2) = - 0.0153d0*eta1 + 1.0153d0*eta2
     $    + al_s*(eta1*(0.0865d0 + 0.3007d0*eta)
     $    + eta2*(-7.7870d0 + 7.3999d0*eta))
c      evolution from mu_t to mu_B
      al_s = alfas_nlo(amu_b)/4/pi
      eta = alfas_nlo(qm(3))/alfas_nlo(amu_b)
      eta1 = eta**(6/23.d0)
c      "VLL" sector
      vxx(1) = eta1*(1 + al_s*1.6273d0*(1 - eta))
c      "LR" sector
      eta1 = eta**(3/23d0)
      eta2 = eta**(24/23.d0)
      eta3 = eta**(26/23.d0)
      sxy(1,1,1) = eta1 + al_s*(0.9250d0/eta2 + ( - 2.0994d0
     $    + 1.1744d0*eta)*eta1)
      sxy(1,1,2) = al_s*1.3875d0*(eta3 - 1/eta2)
      sxy(1,2,1) = 2/3.d0*(eta1 - 1/eta2)
     $    + al_s*(( - 11.7329d0 + 0.7829d0*eta)*eta1
     $    + ( - 5.3048d0 + 16.2548d0*eta)/eta2)
      sxy(1,2,2) = 1/eta2 + al_s*(0.9250d0*eta3
     $    + (7.9572d0 - 8.8822d0*eta)/eta2)
c      "SLL" sector
      eta1 = 1/eta**0.6315d0
      eta2 = eta**0.7184d0
      sxx(1,1,1) = 1.0153d0*eta1 - 0.0153d0*eta2
     $    + al_s*(eta1*(4.8177d0 - 5.2272d0*eta)
     $    + eta2*(0.3371d0 + 0.0724d0*eta))
      sxx(1,1,2) = 1.9325d0*(eta1 - eta2)
     $    + al_s*(eta1*(9.1696d0 - 38.8778d0*eta)
     $    + eta2*(42.5021d0 - 12.7939d0*eta))
      sxx(1,2,1) = - 0.0081d0*(eta1 - eta2)
     $    + al_s*(eta1*(0.0531d0 + 0.0415d0*eta)
     $    + eta2*( - 0.0566d0 - 0.0380d0*eta))
      sxx(1,2,2) = - 0.0153d0*eta1 + 1.0153d0*eta2
     $    + al_s*(eta1*(0.1011d0 + 0.3083d0*eta)
     $    + eta2*( - 7.1314d0 + 6.7219d0*eta))
c      evolution from mu_t to mu_K
      al_s = alfas_nlo(amu_k)/4/pi
      eta = alfas_nlo(qm(3))/alfas_nlo(qm(2))
      etap = alfas_nlo(qm(2))/alfas_nlo(amu_k)
      eta1 = eta**(6/23.d0)
      eta1p = etap**(6/25.d0)
c      "VLL" sector
      vxx(2) = eta1*eta1p*(1 + al_s*(1.7917d0 - 0.1644d0*etap
     $    - 1.6273d0*eta*etap))
c      "LR" sector
      eta1 = eta**(3/23d0)
      eta2 = eta**(24/23.d0)
      eta3 = eta**(26/23.d0)
      eta1p = etap**(3/25d0)
      eta2p = etap**(24/25.d0)
      eta3p = etap**(28/25.d0)
      sxy(2,1,1) = eta1*eta1p + al_s*(0.9279d0/eta2/eta2p
     $    - 0.0029d0*eta3p/eta2 + ( - 2.0241d0 - 0.0753d0*etap
     $    + 1.1744d0*eta*etap)*eta1*eta1p)
      sxy(2,1,2) = al_s*( -1.3918d0/eta2/eta2p + 0.0043d0*eta3p/eta2
     $    + 1.3875d0*eta3*eta3p)
      sxy(2,2,1) = 2/3.d0*(eta1*eta1p - 1/eta2/eta2p)
     $    + al_s*( -0.0019d0*eta3p/eta2 + 5.0000d0*etap**(1.d0/25)*eta1
     $    + (- 16.6828d0 - 0.0502d0*etap + 0.7829d0*eta*etap)*eta1*eta1p
     $    + ( - 4.4701d0 - 0.8327d0*etap
     $    + 16.2548d0*eta*etap)/eta2/eta2p)
      sxy(2,2,2) = 1/eta2/eta2p + al_s*(0.0029d0*eta3p/eta2
     $    + 0.9250d0*eta3*eta3p  + (6.7052d0 + 1.2491d0*etap
     $    - 8.8822d0*eta*etap)/eta2/eta2p)
c      "SLL" sector
      eta1 = 1/eta**0.6315d0
      eta2 = eta**0.7184d0
      eta1p = 1/etap**0.5810d0
      eta2p = etap**0.6610d0
      sxx(2,1,1) = 1.0153d0*eta1*eta1p - 0.0153d0*eta2*eta2p
     $    + al_s*(etap*(0.0020d0*eta1*eta2p - 0.0334d0*eta1p*eta2)
     $    + eta1*eta1p*(4.2458d0 + 0.5700d0*etap - 5.2272d0*eta*etap)
     $    + eta2*eta2p*(0.3640d0 + 0.0064d0*etap + 0.0724d0*eta*etap))
      sxx(2,1,2) = 1.9325d0*(eta1*eta1p - eta2*eta2p)
     $    + al_s*(etap*(0.0038d0*eta2p*eta1 - 4.2075d0*eta1p*eta2)
     $    + eta1*eta1p*(8.0810d0 + 1.0848d0*etap - 38.8778d0*eta*etap)
     $    + eta2*eta2p*(45.9008d0 + 0.8087d0*etap - 12.7939d0*eta*etap))
      sxx(2,2,1) = - 0.0081d0*(eta1*eta1p - eta2*eta2p)
     $    + al_s*(etap*( - 0.0011d0*eta2p*eta1 + 0.0003d0*eta1p*eta2)
     $    + eta2*eta2p*( - 0.0534d0 - 0.0034d0*etap - 0.0380d0*eta*etap)
     $    + eta1*eta1p*(0.0587d0 - 0.0045d0*etap + 0.0415d0*eta*etap))
      sxx(2,2,2) = - 0.0153d0*eta1*eta1p + 1.0153d0*eta2*eta2p
     $    + al_s*(etap*(- 0.0020d0*eta2p*eta1 + 0.0334d0*eta1p*eta2)
     $    + eta1*eta1p*(0.1117d0 - 0.0086d0*etap + 0.3083d0*eta*etap)
     $    + eta2*eta2p*(- 6.7398d0 - 0.4249d0*etap + 6.7219d0*eta*etap))
c     QCD formfactors for SM-like formfactor vector-LL initialized
c     separately - see routine BLOCK DATA init_4q
c     initialization of QCD factors finished:
      init_eta = .false.
      return
      end

      double complex function dd_vll_sm_wil(i,j)
c     Full A^V_LL SM formfactor at low energy
      implicit double precision (a-h,o-z)
      double complex vc,vt
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
      common/sm_4q/eta_cc,eta_ct,eta_tt,eta_b,bk_sm,bd_sm,bb_sm(2)
      common/fermi/g_fermi
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      xt = umu(3)*umu(3)/wm2
      xc = uml(2)*uml(2)/wm2
      vc = dconjg(ckm_phys(2,j))*ckm_phys(2,i)
      vt = dconjg(ckm_phys(3,j))*ckm_phys(3,i)
c      check if we calculate K or B meson mixing, i.e. if we are
c      at m_s or at m_b mass scale
      if (max(i,j).lt.3) then
         dd_vll_sm_wil = eta_cc*vc*vc*seq_dd(xc)
     $        + 2*eta_ct*vc*vt*sneq_dd(xc,xt) + eta_tt*vt*vt*seq_dd(xt)
      else
        dd_vll_sm_wil = eta_b*(vc*vc*seq_dd(xc) + 2*vc*vt*sneq_dd(xc,xt)
     $      + vt*vt*seq_dd(xt))
      end if
      dd_vll_sm_wil = (g_fermi*wm/2/pi)**2*dd_vll_sm_wil
      return
      end

      double precision function seq_dd(x)
      implicit double precision (a-h,o-z)
      seq_dd = x/4*(1 + 9/(1 - x) - 6/(1 - x)/(1 - x))
     $    - 3.d0/2*x*x*x/(1 - x)/(1 - x)/(1 - x)*log(x)
      return
      end

      double precision function sneq_dd(x,y)
      implicit double precision (a-h,o-z)
      sneq_dd = x/4*((y*y - 8*y + 4)/(1 - y)/(1 - y)*log(y)
     $    - 3*y/(1 - y) - 4*log(x))
      return
      end

      subroutine dd_wil_run(i,j)
c     CAUTION: SM contribution calculated SEPARATELY in dd_vll_sm_wil
      implicit double precision (a-h,o-z)
      logical init_eta
      double complex tmp
      double complex dd_vll,dd_vrr,dd_sll,dd_srr,dd_vlr,dd_slr,dd_tl,
     $     dd_tr
c      double complex dd_vll_sm_wil
      double complex dd_sll_yuk,dd_srr_yuk,dd_slr_yuk
      double complex dds_vll,dds_vrr,dds_lr,dds_sll,dds_srr
      common/dd_wil_coeff/dds_vll,dds_vrr,dds_lr(2),
     $    dds_sll(2),dds_srr(2)
      common/ev_mat_4q/vxx(4),sxx(4,2,2),sxy(4,2,2),init_eta
      common/diag_type/ihpeng,izpeng,ifpeng,ibox
      if (init_eta) call eta_4q_evol()
c      check if we calculate K or B meson mixing, i.e. if we are
c      at mu=2 GeV (k=1) or at mu=m_b (k=3) mass scale
      if (max(i,j).lt.3) then
        k = 2
      else
        k = 1
      end if
c      transition to the Buras et al. operator basis:
c      almost identical to ours, apart from the sign of the
c      tensor operator: ours O_{L,R}^T = (-1) x their O_2^{SLL,SRR}
      dds_vll    = dd_vll(i,j)
      dds_vrr    = dd_vrr(i,j)
      dds_sll(1) = dd_sll(i,j)
      dds_sll(2) = - dd_tl(i,j)
      dds_srr(1) = dd_srr(i,j)
      dds_srr(2) = - dd_tr(i,j)
      dds_lr(1)  = dd_vlr(i,j)
      dds_lr(2)  = dd_slr(i,j)
      if (ihpeng.eq.1) then
         dds_sll(1) = dds_sll(1) + ihpeng*dd_sll_yuk(i,j)
         dds_srr(1) = dds_srr(1) + ihpeng*dd_srr_yuk(i,j)
         dds_lr(2)  = dds_lr(2) + ihpeng*dd_slr_yuk(i,j)
      end if
c      NLO evolution of the MSSM part
c      First step: from mu = M_SUSY = (M_gluino + aver M_D)/2 to mu = m_t
      dds_vll  = vxx(3)*dds_vll
      dds_vrr  = vxx(3)*dds_vrr
      tmp = dds_sll(1)
      dds_sll(1) = sxx(3,1,1)*dds_sll(1) + sxx(3,1,2)*dds_sll(2)
      dds_sll(2) = sxx(3,2,1)*tmp        + sxx(3,2,2)*dds_sll(2)
      tmp = dds_srr(1)
      dds_srr(1) = sxx(3,1,1)*dds_srr(1) + sxx(3,1,2)*dds_srr(2)
      dds_srr(2) = sxx(3,2,1)*tmp        + sxx(3,2,2)*dds_srr(2)
      tmp = dds_lr(1)
      dds_lr(1)  = sxy(3,1,1)*dds_lr(1)  + sxy(3,1,2)*dds_lr(2)
      dds_lr(2)  = sxy(3,2,1)*tmp        + sxy(3,2,2)*dds_lr(2)
c      Second step: from mu = m_t to mu = m_b or mu = 2GeV
      dds_vll  = vxx(k)*dds_vll
      dds_vrr  = vxx(k)*dds_vrr
      tmp = dds_sll(1)
      dds_sll(1) = sxx(k,1,1)*dds_sll(1) + sxx(k,1,2)*dds_sll(2)
      dds_sll(2) = sxx(k,2,1)*tmp        + sxx(k,2,2)*dds_sll(2)
      tmp = dds_srr(1)
      dds_srr(1) = sxx(k,1,1)*dds_srr(1) + sxx(k,1,2)*dds_srr(2)
      dds_srr(2) = sxx(k,2,1)*tmp        + sxx(k,2,2)*dds_srr(2)
      tmp = dds_lr(1)
      dds_lr(1)  = sxy(k,1,1)*dds_lr(1)  + sxy(k,1,2)*dds_lr(2)
      dds_lr(2)  = sxy(k,2,1)*tmp        + sxy(k,2,2)*dds_lr(2)
c      NLO evolution of the MSSM part
c      full VLL coefficient: SM + MSSM. One can add it here, but
c      SM part goes with different B_K(B_B) factor, so we do it
c      later, in dd_kaon and dd_bmeson routines.
c      If for some reason you wish to add it here, uncomment the
c      line below:
c      dds_vll = dds_vll + dd_vll_sm_wil(i,j)
      return
      end

      subroutine uu_wil_run(i,j)
      implicit double precision (a-h,o-z)
      logical init_eta
      double complex uu_vll,uu_vrr,uu_sll,uu_srr,uu_vlr,uu_slr,uu_tl,
     $     uu_tr
      double complex tmp
      double complex uus_vll,uus_vrr,uus_lr,uus_sll,uus_srr
      common/uu_wil_coeff/uus_vll,uus_vrr,uus_lr(2),
     $    uus_sll(2),uus_srr(2)
      common/ev_mat_4q/vxx(4),sxx(4,2,2),sxy(4,2,2),init_eta
      if (init_eta) call eta_4q_evol()
c      We can consider only D=(c \bar u) mesons, so that we are
c      at mu = 2 GeV mass scale:
      k=2
c      Just in case, check arguments
      if (i*j.ne.2) stop 'Wrong indices in uu_wil_run'
c      transition to the Buras et al. operator basis:
c      almost identical to ours, apart from the sign of the
c      tensor operator: ours O_{L,R}^T = (-1) x their O_2^{SLL,SRR}
      uus_vll    = uu_vll(i,j)
      uus_vrr    = uu_vrr(i,j)
      uus_sll(1) = uu_sll(i,j)
      uus_sll(2) = - uu_tl(i,j)
      uus_srr(1) = uu_srr(i,j)
      uus_srr(2) = - uu_tr(i,j)
      uus_lr(1)  = uu_vlr(i,j)
      uus_lr(2)  = uu_slr(i,j)
c      NLO evolution of the MSSM part
c      First step: from mu = M_SUSY = (M_gluino + aver M_U)/2 to mu = m_t
      uus_vll  = vxx(4)*uus_vll
      uus_vrr  = vxx(4)*uus_vrr
      tmp = uus_sll(1)
      uus_sll(1) = sxx(4,1,1)*uus_sll(1) + sxx(4,1,2)*uus_sll(2)
      uus_sll(2) = sxx(4,2,1)*tmp        + sxx(4,2,2)*uus_sll(2)
      tmp = uus_srr(1)
      uus_srr(1) = sxx(4,1,1)*uus_srr(1) + sxx(4,1,2)*uus_srr(2)
      uus_srr(2) = sxx(4,2,1)*tmp        + sxx(4,2,2)*uus_srr(2)
      tmp = uus_lr(1)
      uus_lr(1)  = sxy(4,1,1)*uus_lr(1)  + sxy(4,1,2)*uus_lr(2)
      uus_lr(2)  = sxy(4,2,1)*tmp        + sxy(4,2,2)*uus_lr(2)
c      Second step: from mu = m_t to mu = 2 GeV
      uus_vll  = vxx(k)*uus_vll
      uus_vrr  = vxx(k)*uus_vrr
      tmp = uus_sll(1)
      uus_sll(1) = sxx(k,1,1)*uus_sll(1) + sxx(k,1,2)*uus_sll(2)
      uus_sll(2) = sxx(k,2,1)*tmp        + sxx(k,2,2)*uus_sll(2)
      tmp = uus_srr(1)
      uus_srr(1) = sxx(k,1,1)*uus_srr(1) + sxx(k,1,2)*uus_srr(2)
      uus_srr(2) = sxx(k,2,1)*tmp        + sxx(k,2,2)*uus_srr(2)
      tmp = uus_lr(1)
      uus_lr(1)  = sxy(k,1,1)*uus_lr(1)  + sxy(k,1,2)*uus_lr(2)
      uus_lr(2)  = sxy(k,2,1)*tmp        + sxy(k,2,2)*uus_lr(2)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Here are the expressions for the phenomenological quantities c
c     eps_K, Delta_M(K), Delta_M(D), Delta_M(Bs), Delta_M(Bd)      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dd_kaon(eps_k,delta_mk)
c      value of eps_K (with exp(i pi/4) removed) and of delta_mk
      implicit double precision (a-h,o-z)
      double complex heff_mat
      double complex dd_vll_sm_wil
      double complex dds_vll,dds_vrr,dds_lr,dds_sll,dds_srr
      common/dd_wil_coeff/dds_vll,dds_vrr,dds_lr(2),
     $    dds_sll(2),dds_srr(2)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/meson_data/dmk,amk,epsk,fk,dmd,amd,fd,
     $    amb(2),dmb(2),tau_b(2),fb(2)
      common/bx_4q/bk(5),bd(5),bb(2,5),amu_k,amu_d,amu_b
      common/sm_4q/eta_cc,eta_ct,eta_tt,eta_b,bk_sm,bd_sm,bb_sm(2)
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
      call dd_wil_run(1,2)
      rm = (amk/(qmass_nlo(dml(1),amud(1),amu_k)
     $     + qmass_nlo(dml(2),amud(2),amu_k)))**2
      heff_mat = amk*fk*fk/24*(8*bk_sm*dd_vll_sm_wil(1,2)
     $    + 8*bk(1)*(dds_vll + dds_vrr)
     $    - 5*bk(2)*rm*(dds_sll(1) + dds_srr(1))
     $    - 12*bk(3)*rm*(dds_sll(2) + dds_srr(2))
     $    - 4*bk(4)*rm*dds_lr(1) + 6*bk(5)*rm*dds_lr(2))
      eps_k = - dimag(heff_mat)/dmk/sq2
      delta_mk = 2*dble(heff_mat)
      return
      end

      subroutine uu_dmeson(delta_md)
      implicit double precision (a-h,o-z)
      double complex heff_mat
      double complex uus_vll,uus_vrr,uus_lr,uus_sll,uus_srr
      common/uu_wil_coeff/uus_vll,uus_vrr,uus_lr(2),
     $    uus_sll(2),uus_srr(2)
      common/meson_data/dmk,amk,epsk,fk,dmd,amd,fd,
     $    amb(2),dmb(2),tau_b(2),fb(2)
      common/bx_4q/bk(5),bd(5),bb(2,5),amu_k,amu_d,amu_b
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
      call uu_wil_run(1,2)
      rm = (amd/(qmass_nlo(uml(1),amuu(1),amu_d)
     $     + qmass_nlo(uml(2),amuu(2),amu_d)))**2
      heff_mat = amd*fd*fd/24*(8*bd(1)*(uus_vll + uus_vrr)
     $    - 5*bd(2)*rm*(uus_sll(1) + uus_srr(1))
     $    - 12*bd(3)*rm*(uus_sll(2) + uus_srr(2))
     $    - 4*bd(4)*rm*uus_lr(1) + 6*bd(5)*rm*uus_lr(2))
      delta_md = 2*dble(heff_mat)
      return
      end

      subroutine dd_bmeson(i,delta_mb,dmb_re,dmb_im)
c      i = 1: Bd mesons; i=2: Bs mesons
      implicit double precision (a-h,o-z)
      double complex heff_mat,heff_mat_sm,heff_mat_np
      double complex dd_vll_sm_wil
      double complex dds_vll,dds_vrr,dds_lr,dds_sll,dds_srr
      common/dd_wil_coeff/dds_vll,dds_vrr,dds_lr(2),
     $    dds_sll(2),dds_srr(2)
      common/meson_data/dmk,amk,epsk,fk,dmd,amd,fd,
     $    amb(2),dmb(2),tau_b(2),fb(2)
      common/bx_4q/bk(5),bd(5),bb(2,5),amu_k,amu_d,amu_b
      common/sm_4q/eta_cc,eta_ct,eta_tt,eta_b,bk_sm,bd_sm,bb_sm(2)
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
c     in the common /heff_bmeson/ below:
c       heff_mat is full matrix element, \Delta M_12 
c       heff_mat_sm is SM matrix element, \Delta M_12^SM 
c       heff_mat_np is "new physics" matrix element, \Delta M_12^NP
c       iq = 1: Bd mesons; iq=2: Bs mesons
      common/heff_bmeson/heff_mat,heff_mat_sm,heff_mat_np,iq
      iq = i
      call dd_wil_run(i,3)
      rm = (amb(i)/(qmass_nlo(dml(i),amud(i),amu_b)
     $     + qmass_nlo(dml(3),amud(3),amu_b)))**2
      heff_mat_sm = amb(i)*fb(i)*fb(i)/3*bb_sm(i)*dd_vll_sm_wil(i,3)
      heff_mat_np = amb(i)*fb(i)*fb(i)/24*(8*bb(i,1)*(dds_vll + dds_vrr)
     $     - 5*bb(i,2)*rm*(dds_sll(1) + dds_srr(1)) 
     $     - 12*bb(i,3)*rm*(dds_sll(2) + dds_srr(2)) 
     $     - 4*bb(i,4)*rm*dds_lr(1) + 6*bb(i,5)*rm*dds_lr(2))
      heff_mat = heff_mat_sm + heff_mat_np
      delta_mb = 2*abs(heff_mat)
      dmb_re = dble(heff_mat)
      dmb_im = dimag(heff_mat)
      return
      end

ccccccccccccccccccccccccccc
c      SM expressions:    c
ccccccccccccccccccccccccccc

      subroutine dd_kaon_sm(eps_k,delta_mk)
c      value of eps_K (with exp(i pi/4) removed) and of delta_mk
c      SM part only!
      implicit double precision (a-h,o-z)
      double complex heff_mat,dd_vll_sm_wil
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/meson_data/dmk,amk,epsk,fk,dmd,amd,fd,
     $    amb(2),dmb(2),tau_b(2),fb(2)
      common/sm_4q/eta_cc,eta_ct,eta_tt,eta_b,bk_sm,bd_sm,bb_sm(2)
      heff_mat = amk*fk*fk*bk_sm/3*dd_vll_sm_wil(1,2)
      eps_k = - dimag(heff_mat)/dmk/sq2
      delta_mk = 2*dble(heff_mat)
      return
      end

      subroutine uu_dmeson_sm(delta_md)
      implicit double precision (a-h,o-z)
      logical init_eta
      double complex uu_vll
      common/meson_data/dmk,amk,epsk,fk,dmd,amd,fd,
     $    amb(2),dmb(2),tau_b(2),fb(2)
      common/ev_mat_4q/vxx(4),sxx(4,2,2),sxy(4,2,2),init_eta
      common/sm_4q/eta_cc,eta_ct,eta_tt,eta_b,bk_sm,bd_sm,bb_sm(2)
      if (init_eta) call eta_4q_evol()
      delta_md = 2/3.d0*amd*fd*fd*bd_sm*vxx(2)*dble(uu_vll(1,2))
      return
      end

      subroutine dd_bmeson_sm(i,delta_mb,dmb_re,dmb_im)
c      i = 1: Bd mesons; i=2: Bs mesons
c      SM part only!
      implicit double precision (a-h,o-z)
      double complex dd_vll_sm_wil,heff_mat
      common/meson_data/dmk,amk,epsk,fk,dmd,amd,fd,
     $    amb(2),dmb(2),tau_b(2),fb(2)
      common/sm_4q/eta_cc,eta_ct,eta_tt,eta_b,bk_sm,bd_sm,bb_sm(2)
      heff_mat = amb(i)*fb(i)*fb(i)*bb_sm(i)/3*dd_vll_sm_wil(i,3)
      delta_mb = 2*abs(heff_mat)
      dmb_re = dble(heff_mat)
      dmb_im = dimag(heff_mat)
      return
      end

      subroutine init_alpha_s_susy
      implicit double precision (a-h,o-z)
      double complex zu,zd,gm2,gm3
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      external init_4q
      asusy = 0
      do i=1,6
         asusy = asusy + sdm(i)
      end do
      asusy = gm1/2 + asusy/12
      g3d = 2*sqrt(pi*alfas_nlo(asusy))
      asusy = 0
      do i=1,6
         asusy = asusy + sum(i)
      end do
      asusy = gm1/2 + asusy/12
      g3u = 2*sqrt(pi*alfas_nlo(asusy))
      init_alpha_susy = .false.
      return
      end

