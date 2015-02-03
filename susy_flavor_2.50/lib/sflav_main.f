c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: SFLAV_MAIN.F
c     Released: 26:10:2013(J.R.)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Main routine of SUSY_FLAVOR library calculating available       c
c     physical observables                                            c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine susy_flavor
c     calculate available SUSY_FLAVOR observables
c     store in common blocks for output routine
      implicit double precision (a-h,o-z)
      common/sflav_df0/edmn,edml(3),gminus2(3)
      common/sflav_df1/br_mu_egamma,br_tau_egamma,br_tau_mugamma,
     $     br_k0,br_kp,br_taunu,dtaunu_ratio,dstaunu_ratio,bxgamma,
     $     br_bdll(3,3),br_bsll(3,3),br_tuh,br_tch
      common/sflav_df2/eps_k,delta_mk,delta_md,
     $     delta_mbd,dmbd_re,dmbd_im,
     $     delta_mbs,dmbs_re,dmbs_im

c     Delta F = 0 processes
      do i=1,3
         edml(i) = edm_l(i)     ! lepton EDMs
         gminus2(i) = g_minus_2_anomaly(i) ! lepton (g-2)/2, SUSY contribution
      end do
      edmn = edm_n()            ! neutron EDM

c     Delta F = 1 processes
c     l^J->l^I gamma decays:
      br_mu_egamma   = br_llg(2,1) ! Br(mu-> e gamma)
      br_tau_egamma  = br_llg(3,1) ! Br(tau-> e gamma)
      br_tau_mugamma = br_llg(3,2) ! Br(tau-> mu gamma)
c     K->neutrino decays:
      call k_pivv(br_k0,br_kp)  ! BR(K_L^0 -> pi^0 vv), BR(K^+ -> pi^+ vv)
c     Leptonic B decays:
      do i=1,3
         do j=1,3
            br_bdll(i,j) = b_ll(3,1,i,j) ! BR(B_d -> l^I+ l^J-)
            br_bsll(i,j) = b_ll(3,2,i,j) ! BR(B_s -> l^I+ l^J-)
         end do 
      end do
c     B -> tau decays: BR(B -> tau nu), 
c     BR(B -> D(D*) tau nu)/BR(B -> D(D*) l nu)
      call b_taunu(br_taunu,dtaunu_ratio,dstaunu_ratio) 
c     B->X_s gamma decay:
      delb = 0.99d0             ! Photon energy infrared cutoff
      amiu_b= 4.8d0             ! Renormalization scale miu_b
      bxgamma = bxg_nl(delb,amiu_b) ! BR(B -> X_S gamma)

c     top -> Higgs+quark decays
      br_tuh = br_suu(2,1)      ! BR(t -> u h)
      br_tch = br_suu(2,2)      ! BR(t -> c h)
c      gam_tuh = gam_suu(2,1)    ! Gamma(t -> u h)
c      gam_tch = gam_suu(2,2)    ! Gamma(t -> c h)

c     Delta F = 2 processes
c     KK mixing:
      call dd_kaon(eps_k,delta_mk) ! epsilon_K, Delta m_K
c     Result in ps^-1: delta_mk*gev_s/1.d12
c     DD mixing:
      call uu_dmeson(delta_md)  ! Delta m_D
c     BB mixing:
      call dd_bmeson(1,delta_mbd,dmbd_re,dmbd_im) ! ,Delta m_Bd, Re(H_eff_Bd), Im(H_eff_Bd)
      call dd_bmeson(2,delta_mbs,dmbs_re,dmbs_im) ! ,Delta m_Bs, Re(H_eff_Bs), Im(H_eff_Bs)
      return
      end

