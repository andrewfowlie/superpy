c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: QCD_FUN.F
c     Released: 28:10:1992(J.R.)
c     Revised:  21:02:2001(J.R.)
c     Error in calculations of Lambda QCD for 5 flavours corrected 
c     Revised:  02:07:2007(J.R.)
c     Fixed compatibility issues between FCNC and Higgs libraries

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Set of auxiliary functions for the QCD corrections to the widths    c
c     and branching ratios of Higgs boson decays.                         c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer function nflav(q)
c     Running number of active flavors
      implicit double precision (a-h,o-z)
      common/lamqcd/qcd4,qcd5,qcd6,qstep(6),qerr
      nflav = 0
      do i=1,6
         if (q.ge.qstep(i)) nflav = nflav + 1
      end do
      return
      end

      subroutine qstep_update
c     update squark masses in qstep array(s)
      implicit double precision (a-h,o-z)
      common/lamqcd/qcd4,qcd5,qcd6,qstep(6),qerr
      common/lamqcd_nlo/qcd4_nlo,qcd5_nlo,qcd6_nlo,qstep_nlo(6),qerr_nlo
      common/fmass/em(3),um(3),dm(3)
      qstep(1) = min(um(1),dm(1))
      qstep(2) = max(um(1),dm(1))
      qstep(3) = dm(2)
      qstep(4) = um(2)
      qstep(5) = dm(3)
      qstep(6) = um(3)
      do i=1,6
         qstep_nlo(i) = qstep(i)
      end do 
      return
      end

      double precision function alfas(q)
c     Alpha QCD calculation for nf flavors
c     Third order of perturbation expansion
      implicit double precision (a-h,o-z)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      nf = nflav(q)
      b0 = 11 - 2*nf/3.d0
      b1 = 102 - 38*nf/3.d0
      b2 = 1428.5d0 - 5033*nf/18.d0 + 325*nf*nf/54.d0
      al = 2*log(q/qcd_lam(nf))
      alfas = 4*pi/b0/al*(1 - b1/b0/b0*log(al)/al
     1      + (b1/b0/b0*log(al)/al)**2*(1 - 1/log(al))
     2      + (b0*b2 - b1*b1)/b0**4/al/al)
      return
      end

      double precision function qcd_lam(nf)
c     Lambda QCD for nf active flavors
      implicit double precision (a-h,o-z)
      common/lamqcd/qcd4,qcd5,qcd6,qstep(6),qerr
      common/fmass/em(3),um(3),dm(3)
      external init_qcd
      goto (1,1,1,1,2,3),nf
1     if (dm(3).ne.qstep(5)) call lambda4(dm(3))
      qcd_lam = qcd4
      return
2     qcd_lam = qcd5
      return
3     if (um(3).ne.qstep(6)) call lambda6(um(3))
      qcd_lam = qcd6
      return
      end

      subroutine lam_fit(alfas)
c     Lambda QCD calculation for 5 flavors and given alpha_s(Mz)
      implicit double precision (a-h,o-z)
      common/lamqcd/qcd4,qcd5,qcd6,qstep(6),qerr
      common/fmass/em(3),um(3),dm(3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      call qstep_update         !     update array qstep
c     calculate lambda for 4,5,6 flavours
      nf = 5
      b0 = 11 - 2*nf/3.d0
      b1 = 102 - 38*nf/3.d0
      b2 = 1428.5d0 - 5033*nf/18.d0 + 325*nf*nf/54.d0

      x_old = 2*log(zm/qcd5)
      do 10 i=1,30
        x_new = 4*pi/b0/alfas*(1 - b1/b0/b0*log(x_old)/x_old
     1      + (b1/b0/b0*log(x_old)/x_old)**2*(1 - 1/log(x_old))
     2      + (b0*b2 - b1*b1)/b0**4/x_old/x_old)
        if (abs(x_new - x_old).le.qerr) then
c     Update new qcd4,5,6 values
          qcd5 = zm*exp( - x_new/2)
          call lambda4(dm(3))
          call lambda6(um(3))
          return
        end if
10      x_old = x_new
      stop 'Error in lambda QCD calculations!'
      end

      subroutine lambda6(tm)
c     Lambda QCD calculation for 6 flavors and given top mass
      implicit double precision (a-h,o-z)
      common/lamqcd/qcd4,qcd5,qcd6,qstep(6),qerr
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      if (tm.le.zm) stop 'alpha_s is not calculated properly for mt<MZ!'
      nf = 5
      b05 = 11 - 2*nf/3.d0
      b15 = 102 - 38*nf/3.d0
      b25 = 1428.5d0 - 5033*nf/18.d0 + 325*nf*nf/54.d0
      al = 2*log(tm/qcd5)
      a5  = 1.d0/b05/al*(1 - b15/b05/b05*log(al)/al
     1      + (b15/b05/b05*log(al)/al)**2*(1 - 1/log(al))
     2      + (b05*b25 - b15*b15)/b05**4/al/al)
      nf = 6
      b06 = 11 - 2*nf/3.d0
      b16 = 102 - 38*nf/3.d0
      b26 = 1428.5d0 - 5033*nf/18.d0 + 325*nf*nf/54.d0

      x_old = al
      do 10 i=1,30
        x_new = 1/b06/a5*(1 - b16/b06/b06*log(x_old)/x_old
     1      + (b16/b06/b06*log(x_old)/x_old)**2*(1 - 1/log(x_old))
     2      + (b06*b26 - b16*b16)/b06**4/x_old/x_old)
        if (abs(x_new - x_old).le.qerr) then
c     Update new top mass and qcd6 values
          qstep(6) = tm
          qcd6 = tm*exp( - x_new/2)
          return
        end if
10      x_old = x_new
      stop 'Error in lambda6 QCD calculations!'
      end

      subroutine lambda4(bm)
c     Lambda QCD calculation for 4 flavors and given lambda5 flavors
      implicit double precision (a-h,o-z)
      common/lamqcd/qcd4,qcd5,qcd6,qstep(6),qerr
      nf = 5
      b05 = 11 - 2*nf/3.d0
      b15 = 102 - 38*nf/3.d0
      b25 = 1428.5d0 - 5033*nf/18.d0 + 325*nf*nf/54.d0
      al = 2*log(bm/qcd5)
      a5  = 1.d0/b05/al*(1 - b15/b05/b05*log(al)/al
     1      + (b15/b05/b05*log(al)/al)**2*(1 - 1/log(al))
     2      + (b05*b25 - b15*b15)/b05**4/al/al)
      nf = 4
      b04 = 11 - 2*nf/3.d0
      b14 = 102 - 38*nf/3.d0
      b24 = 1428.5d0 - 5033*nf/18.d0 + 325*nf*nf/54.d0

      x_old = al
      do 10 i=1,30
        x_new = 1/b04/a5*(1 - b14/b04/b04*log(x_old)/x_old
     1      + (b14/b04/b04*log(x_old)/x_old)**2*(1 - 1/log(x_old))
     2      + (b04*b24 - b14*b14)/b04**4/x_old/x_old)
        if (abs(x_new - x_old).le.qerr) then
c     Update new qcd4 values
          qstep(5) = bm
          qcd4 = bm*exp( - x_new/2)
          return
        end if
10      x_old = x_new
      stop 'Error in lambda4 QCD calculations!'
      end

      block data init_qcd
      implicit double precision (a-h,o-z)
      common/dzeta3/g0,dz3
      common/lamqcd/qcd4,qcd5,qcd6,qstep(6),qerr
      common/lamqcd_nlo/qcd4_nlo,qcd5_nlo,qcd6_nlo,qstep_nlo(6),qerr_nlo
      data qcd4,qcd5,qcd6/288.20d-3,208.36d-3,91.838d-3/
      data qstep/2.15d-3,4.7d-3,9.35d-2,1.275d0,4.18d0,163.2d0/
      data qerr/1.d-3/
      data qcd4_nlo,qcd5_nlo,qcd6_nlo/324.81d-3,226.23d-3,95.222d-3/
      data qstep_nlo/2.15d-3,4.7d-3,9.35d-2,1.275d0,4.18d0,163.2d0/
      data qerr_nlo/1.d-3/
      data g0,dz3/-8.d0,438.5513262d0/
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     routines for running quark mass calculations              c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function run_mult(q)
c     Coefficient multiplying m_hat and running quark mass
      implicit double precision (a-h,o-z)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/dzeta3/g0,dz3
      al = alfas(q)
      nf = nflav(q)
      b0 = 11 - 2*nf/3.d0
      b1 = 102 - 38*nf/3.d0
      b2 = 1428.5d0 - 5033*nf/18.d0 + 325*nf*nf/54.d0
      g1 = - 404/3.d0 + 40*nf/9.d0
      g2 = 2/3.d0*(140*nf*nf/27.d0 + dz3*nf - 3747)
      run_mult = (b0*al/2/pi)**(-g0/2/b0)*(1
     1         + (b1*g0 - b0*g1)/b0/b0*al/8/pi
     2         + (((b1*g0 - b0*g1)/b0/b0)**2/2 +g1*b1/b0/b0
     3         + g0*(b0*b2 - b1*b1)/b0/b0/b0 - g2/b0)*(al/8/pi)**2)
      return
      end

      double precision function runmass(q,qm)
c     running quark mass
      implicit double precision (a-h,o-z)
      common/lamqcd/qcd4,qcd5,qcd6,qstep(6),qerr
      x0 = qm/run_mult(qm)
      n0 = nflav(qm)
      nf = max(4,nflav(q))
      if (n0.ne.nf) then
        istep = sign(1,nf - n0)
        eps = istep/1.d6
        do 10 i=n0 + istep,nf,istep
          j = i - (istep - 1)/2
10        x0 = x0*run_mult(qstep(j) - eps)/run_mult(qstep(j) + eps)
      end if
      runmass = x0*run_mult(q)
      return
      end

      double precision function sp_del(iscal,beta)
c     QCD corrections to scalar/pseudoscalar -> qq pair decay
c     iscal = 1 for scalar decays
c     iscal = 2 for pseudoscalar decays
      implicit double precision (a-h,o-z)
      rbt   = (1 - beta)/(1 + beta)
      bsq   = beta*beta
      if (iscal.eq.1) then
        sp_del = aaa(beta)/beta + 3.d0/8/bsq*(7*bsq - 1)
     1        + (3 + 34*bsq - 13*bsq*bsq)/16/beta/bsq*log(1/rbt)
      else if (iscal.eq.2) then
        sp_del = aaa(beta)/beta + 3.d0/8*(7 - bsq)
     1       + (19 + 2*bsq + 3*bsq*bsq)/16/beta*log(1/rbt)
      else
        stop 'Argument in SP_DEL outside range!'
      end if
      return
      end

      double precision function aaa(beta)
c     Common part of QCD corrections to neutral Higgs -> qq pair decay
      implicit double precision (a-h,o-z)
      double complex li2,rbt
      rbt = (1 - beta)/(1 + beta)
      aaa = dble((4*li2(rbt) + 2*li2( - rbt) - 2*log(beta)*log(1/rbt)
     1    - 3*log(2/(1 + beta))*log(1/rbt))*(1 + beta*beta)
     2    - 3*beta*log(4/(1 - beta*beta)) - 4*beta*log(beta))
      return
      end

      double precision function run_gam(qm,hm,iscal)
c     Coefficient multiplying tree H->qq decay with due to the QCD corrections
c     iscal = 1 for scalar decays
c     iscal = 2 for pseudoscalar decays
      implicit double precision (a-h,o-z)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      beta = sqrt(1 - 4*qm*qm/hm/hm)
      run_gam = 1 + 4*alfas(hm)/pi*(sp_del(iscal,beta)/3 - log(qm/hm))
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     another set of procedures calculating alpha_s strictly at NLO   c
c     (to avoid some problems with mixing of expansion orders)        c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function alfas_nlo(q)
c     Alpha QCD calculation for nf flavors at NLO
c     Third order of perturbation expansion
      implicit double precision (a-h,o-z)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      nf = nflav(q)
      b0 = 11 - 2*nf/3.d0
      b1 = 102 - 38*nf/3.d0
      al = 2*log(q/qcd_lam_nlo(nf))
      alfas_nlo = 4*pi/b0/al*(1 - b1/b0/b0*log(al)/al)
      return
      end

      double precision function qcd_lam_nlo(nf)
c     Lambda QCD for nf active flavors
      implicit double precision (a-h,o-z)
      common/lamqcd_nlo/qcd4,qcd5,qcd6,qstep(6),qerr
      common/fmass/em(3),um(3),dm(3)
      external init_qcd
      goto (1,1,1,1,2,3),nf
1     if (dm(3).ne.qstep(5)) call lambda4_nlo(dm(3))
      qcd_lam_nlo = qcd4
      return
2     qcd_lam_nlo = qcd5
      return
3     if (um(3).ne.qstep(6)) call lambda6_nlo(um(3))
      qcd_lam_nlo = qcd6
      return
      end

      subroutine lam_fit_nlo(alfas)
c     NLO Lambda QCD calculation for 5 flavors and given alpha_s(Mz)
      implicit double precision (a-h,o-z)
      common/lamqcd_nlo/qcd4,qcd5,qcd6,qstep(6),qerr
      common/fmass/em(3),um(3),dm(3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      call qstep_update         ! update array qstep
c     calculate lambda for 4,5,6 flavours
      nf = 5
      b0 = 11 - 2*nf/3.d0
      b1 = 102 - 38*nf/3.d0

      x_old = 2*log(zm/qcd5)
      do 10 i=1,30
        x_new = 4*pi/b0/alfas*(1 - b1/b0/b0*log(x_old)/x_old)
        if (abs(x_new - x_old).le.qerr) then
c      Update new qcd4,5,6 values
          qcd5 = zm*exp( - x_new/2)
          call lambda4_nlo(dm(3))
          call lambda6_nlo(um(3))
c     Buras hep-ph/9806471: qcd4/5/6 = 325/226/92 MeV, mu_t = 166 GeV
          return
        end if
10      x_old = x_new
      stop 'Error in NLO lambda QCD calculations!'
      end

      subroutine lambda6_nlo(tm)
c     Lambda QCD calculation for 6 flavors and given top mass
      implicit double precision (a-h,o-z)
      common/lamqcd_nlo/qcd4,qcd5,qcd6,qstep(6),qerr
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      if (tm.le.zm) stop 'alpha_s is not calculated properly for mt<MZ!'
      nf = 5
      b05 = 11 - 2*nf/3.d0
      b15 = 102 - 38*nf/3.d0
      al = 2*log(tm/qcd5)
      a5  = 1.d0/b05/al*(1 - b15/b05/b05*log(al)/al)
      nf = 6
      b06 = 11 - 2*nf/3.d0
      b16 = 102 - 38*nf/3.d0

      x_old = al
      do 10 i=1,30
        x_new = 1/b06/a5*(1 - b16/b06/b06*log(x_old)/x_old)
        if (abs(x_new - x_old).le.qerr) then
c     Update new top mass and qcd6 values
          qstep(6) = tm
          qcd6 = tm*exp( - x_new/2)
          return
        end if
10      x_old = x_new
      stop 'Error in NLO lambda6 QCD calculations!'
      end

      subroutine lambda4_nlo(bm)
c     Lambda QCD calculation for 4 flavors and given lambda5 flavors
      implicit double precision (a-h,o-z)
      common/lamqcd_nlo/qcd4,qcd5,qcd6,qstep(6),qerr
      nf = 5
      b05 = 11 - 2*nf/3.d0
      b15 = 102 - 38*nf/3.d0
      al = 2*log(bm/qcd5)
      a5  = 1.d0/b05/al*(1 - b15/b05/b05*log(al)/al)

      nf = 4
      b04 = 11 - 2*nf/3.d0
      b14 = 102 - 38*nf/3.d0

      x_old = al
      do 10 i=1,30
        x_new = 1/b04/a5*(1 - b14/b04/b04*log(x_old)/x_old)
        if (abs(x_new - x_old).le.qerr) then
c     Update new qcd4 values
          qstep(5) = bm
          qcd4 = bm*exp( - x_new/2)
          return
        end if
10      x_old = x_new
      stop 'Error in NLO lambda4 QCD calculations!'
      end

      double precision function qmass_nlo(qm,amu0,amu)
c     NLO running quark mass
c     qmass_nlo returns running quark mass at scale amu for given qm(amu0)
      implicit double precision (a-h,o-z)
      common/lamqcd_nlo/qcd4,qcd5,qcd6,qstep(6),qerr
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      if (amu0.gt.amu) stop 'qmass_nlo works for only amu>amu0!'
      qmass_nlo = qm
      if (amu0.eq.amu) return
      nf0 = max(1,nflav(amu0))
      nf  = nflav(amu)
      g0 = 8
      do i=nf0,nf
         b0 = 11 - 2*i/3.d0
         b1 = 102 - 38*i/3.d0
         g1 = 4/3.d0*(101 - 10*i/3.d0)
         al0 = alfas_nlo(max(amu0,qstep(i)))
         if (i.lt.6) then
            al1 = alfas_nlo(min(amu,qstep(i+1)))
         else
            al1 = alfas(amu)
         end if
         qmass_nlo = qmass_nlo*(al1/al0)**(g0/2/b0)
     $        *(1 + (g1 - b1*g0/b0)*(al1 - al0)/8/pi/b0)
      end do
      return
      end

      double precision function qm_pole_nlo(qm)
c     calculates pole quark mass for given running MSbar mass qm(qm)
      implicit double precision (a-h,o-z)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/lamqcd/qcd4,qcd5,qcd6,qstep(6),qerr
      al = alfas(qm)/pi
      nl = nflav(qm) - 1
      delm = 0.d0
      do i=1,nl
         delm = delm + 1 - 4/3.d0*qstep(i)/qm
      end do
      qm_pole_nlo = qm*(1 + 4/3.d0*al 
     $     + (13.4434d0 - 1.0414d0*delm)*al*al)
c     $     + (0.6527*nl*nl - 26.655*nl + 190.595)*al*al*al) ! 3rd loop
      return
      end
