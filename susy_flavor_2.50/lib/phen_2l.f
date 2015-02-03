c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: PHEN_2L.F
c     Released: 11:02:2012 (J.R.)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expressions for the phenomenological      c
c     quantities involving two external leptons like mu->e gamma,  c
c     lepton EDM, g-2                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function gam_llg(j,i)
c     l^J->l^I decay width
      implicit double precision (a-h,o-z)
      double complex cfl(5),cfr(5)
      double complex cl,cr
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
c     fermion masses
      common/fmass/em(3),um(3),dm(3)
      if (em(j).le.em(i)) then
         gam_llg = 0
         return
      end if
      call ll_gam(i,j,cfl,cfr)
      cl = cfl(1) + em(i)*cfl(4) - em(j)*cfr(4)
      cr = cfr(1) - em(j)*cfl(4) + em(i)*cfr(4)
      gam_llg = e2*em(j)**3/4/pi*(abs(cl)**2 + abs(cr)**2)
      return
      end

      double precision function br_llg(j,i)
c     BR(l^J->l^I) decay
      implicit double precision (a-h,o-z)
      double complex cfl(5),cfr(5)
      double complex cl,cr
      double precision em,um,dm,g_fermi
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/fmass/em(3),um(3),dm(3)
      common/fermi/g_fermi
      common/tau_gam/br_tau_evv
      if (em(j).le.em(i)) then
         br_llg = 0
         return
      end if
      call ll_gam(i,j,cfl,cfr)
      cl = cfl(1) + em(i)*cfl(4) - em(j)*cfr(4)
      cr = cfr(1) - em(j)*cfl(4) + em(i)*cfr(4)
c
c      write(*,*)i,j
c      write(*,*)32*pi*pi*st2/em(j)/e2*cl ! c7l
c      write(*,*)32*pi*pi*st2/em(j)/e2*cr ! c7r
c      write(*,*)
c
      br_llg = 3*(4*pi*e/em(j)/g_fermi)**2*(abs(cl)**2 + abs(cr)**2)
      if (j.eq.3) br_llg = br_llg*br_tau_evv
      return
      end

      double precision function g_minus_2_anomaly(i)
c     magnetic moment anomaly of lepton, a=(g-2)/2
      implicit double precision (a-h,o-z)
      double complex cfl(5),cfr(5)
      common/fmass/em(3),um(3),dm(3)
      call ll_gam(i,i,cfl,cfr)
      g_minus_2_anomaly = - 4*em(i)*dble(cfl(1) 
     $     + em(i)*(cfl(4) - cfr(4)))
      return
      end

      double precision function edm_l(i)
c     electric dipole moment of lepton
      implicit double precision (a-h,o-z)
      double complex cfl(5),cfr(5)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/ph_units/hbar,gev_cm,gev_s
      external init_units
      call ll_gam(i,i,cfl,cfr)
      edm_l = - 2*e*dimag(cfl(1))/e/gev_cm
      return
      end

