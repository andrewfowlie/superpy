c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: BSG_NL.F
c     Released: 12:02:1997 (J.R.)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expressions for the Br(B->X_s gamma)      c
c     in the non-leading approximation.                            c
c     See Chetyrkin,Misiak,Munz hep-ph/9612313                     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function bxg_nl(delta,bm)
c     MSSM result
c     (1 - delta) - relative minimal photon energy
c     bm          - renormalization scale mu_b 
      implicit double precision (a-h,o-z)
      double complex r
      logical bxg_init
      common/bsg_dat/al1,al2,br_cev,al_em,xi3,r(8),gam(8),bxg_init
      common/qmass_pole/um(3),dm(3)
      external init_bsg
      if (.not.bxg_init) call bsg_init
      z = um(2)*um(2)/dm(3)/dm(3)
      d_np = 6*al2*((1 - z)**3/g_bsg(z) - 1)
      bxg_nl = max(br_cev*(1 + d_np/dm(3)/dm(3))*rq_bsg(delta,bm),0.d0)
      return
      end

      double precision function rq_bsg(delta,bm)
c     Br(B->X_s gamma) in the non-leading order, r_quark function
      implicit double precision (a-h,o-z)
      double complex r
      logical bxg_init
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      common/vpar/st,ct,st2,ct2,sct,sct2,q,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/bsg_dat/al1,al2,br_cev,al_em,xi3,r(8),gam(8),bxg_init
      common/qmass_pole/um(3),dm(3)
      z = um(2)*um(2)/dm(3)/dm(3)
      rq_bsg = 6*al_em/pi/g_bsg(z)*f_bsg(z,dm(3))
     $     * abs(ckm_phys(3,2)*ckm_phys(3,3)/ckm_phys(2,3))**2
     $     * (d2_bsg(bm) + a_bsg(delta,bm) + del_bsg(bm))
      return
      end

      double precision function d2_bsg(bm)
c     Br(B->X_s gamma) in the non-leading order, D function squared
      implicit double precision (a-h,o-z)
      double complex r,tmp,c7
      double complex c0_eff_l,c1_eff_l,c0_eff_r,c1_eff_r
      logical bxg_init
      common/qmass_pole/um(3),dm(3)
      common/bsg_dat/al1,al2,br_cev,al_em,xi3,r(8),gam(8),bxg_init
      common/vpar/st,ct,st2,ct2,sct,sct2,q,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      zl = log(dm(3)/bm)
      c7 = c0_eff_l(7,bm) 
      tmp = c1_eff_l(7,bm)
      do i=1,8
         if (i.ne.7) then
            tmp = tmp + (r(i) + gam(i)*zl)*c0_eff_l(i,bm)
         else
            tmp = tmp + (r(7) + gam(7)*zl)*c7
         end if
      end do
      d2_bsg = abs(c7)**2 + alfas(bm)/2/pi*dble(tmp*dconjg(c7))
c     Additional contribution from the opposite chirality operators, 
c     negligible in SM
      c7 = c0_eff_r(bm)  
      d2_bsg = d2_bsg + abs(c7)**2 
     $     + alfas(bm)/2/pi*dble(c1_eff_r(bm)*dconjg(c7))
      return
      end

      double precision function a_bsg(delta,bm)
c     Br(B->X_s gamma) in the non-leading order, A function
      implicit double precision (a-h,o-z)
      double complex c0_eff_l
      logical fm_init
      common/bsg_fm/fmc(8,8),d_old,fm_init
      common/vpar/st,ct,st2,ct2,sct,sct2,q,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      al = alfas(bm)/pi
      dl = log(delta)
      a_bsg = 0
      if (fm_init.or.(delta.ne.d_old)) call fm_bsg_init(delta) 
      do j=1,8
         do i=1,j
            a_bsg = a_bsg 
     $           + dble(c0_eff_l(i,bm)*dconjg(c0_eff_l(j,bm)))*fmc(i,j)
         end do
      end do
      a_bsg = (exp( -al*dl*(7 + 2*dl)/3.d0) - 1)*abs(c0_eff_l(7,bm))**2 
     $     + al*a_bsg 
      return
      end

      double precision function del_bsg(bm)
c     Br(B->X_s gamma) in the non-leading order, Delta function
      implicit double precision (a-h,o-z)
      double complex r,c0_eff_l,c7
      logical bxg_init
      common/qmass_pole/um(3),dm(3)
      common/bsg_dat/al1,al2,br_cev,al_em,xi3,r(8),gam(8),bxg_init
      dg_np = (al1 - 9*al2)/2.d0
      dc_np =  - al2/9        
      c7 = c0_eff_l(7,bm) 
      del_bsg = dg_np/dm(3)/dm(3)*abs(c7)**2
     $     + dc_np/um(2)/um(2)*dble(dconjg(c7)*(c0_eff_l(2,bm) 
     $     - c0_eff_l(1,bm)/6))
      return
      end

      double precision function f_bsg(z,dm)
c     Br(B->X_s gamma) in the non-leading order, F function
      implicit double precision (a-h,o-z)
      common/vpar/st,ct,st2,ct2,sct,sct2,q,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      al = alfas(dm)/pi/3
      f_bsg = (1 - 8*al)/(1 - 2*al*h_bsg(z)/g_bsg(z))
      return
      end

      double precision function g_bsg(z)
c     Phase space factor
      implicit double precision (a-h,o-z)
      g_bsg = 1 + z*(- 8 - 12*z*log(z) + z*z*(8 - z))
      return
      end       

      double precision function h_bsg(z)
c     Phase space factor
      implicit double precision (a-h,o-z)
      double complex li2,x1,x2
      common/vpar/st,ct,st2,ct2,sct,sct2,q,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      z2 = z*z
      x1 = dcmplx(z,0.d0)
      x2 = dcmplx(sqrt(z),0.d0)
      pi2 = pi*pi
      h_bsg = dble(- (1 - z2)*(75 - 956*z + 75*z2)/12
     $     + z*log(z)*(60 + 270*z - 4*z2 + 17*z*z2)/3
     $     + z2*(log(z))**2*(36 + z2)
     $     + (1 - z2)*log(1 - z)*(17 - 64*z + 17*z2)/3
     $     - 4*log(z)*log(1 - z)*(1 + 30*z2 + z2*z2)
     $     - (1 + 16*z2 + z2*z2)*(6*li2(x1) - pi2)
     $     - 32*z**1.5d0*(1 + z)*(pi2 - 4*li2(x2) + 4*li2(-x2)
     $     - 2*log(z)*log((1 - x2)/(1 + x2))))
      return
      end       

      double complex function c0_eff_l(i,bm)
c     leading order c_i effective cofficients
      implicit double precision (a-h,o-z)
      double complex c7_bsg_l,c8_bsg_l
      common/vpar/st,ct,st2,ct2,sct,sct2,q,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/ceff_pow/p(6),cc(6,8),h(8),a(8),e(8),g(8),f(8)
      eta = alfas(wm)/alfas(bm)
      c0_eff_l = (0.d0,0.d0)
      if (i.ne.7) then
         do j=1,6
             c0_eff_l =  c0_eff_l + cc(j,i)*eta**p(j) 
         end do
         if (i.eq.8) c0_eff_l = 
     $        c0_eff_l + (c8_bsg_l() + 313063.d0/363036)*eta**(14.d0/23)
      else
         do j=1,8
            c0_eff_l = c0_eff_l + h(j)*eta**a(j)
         end do
         c0_eff_l = c0_eff_l + eta**(16.d0/23)*c7_bsg_l() 
     $        + 8.d0/3*(eta**(14.d0/23) - eta**(16.d0/23))*c8_bsg_l()
      end if
      return
      end

      double complex function c1_eff_l(i,bm)
c     non-leading order c_7 effective cofficients at bm scale
      implicit double precision (a-h,o-z)
      double complex c7_bsg_l,c8_bsg_l,z,li2
      common/vpar/st,ct,st2,ct2,sct,sct2,q,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/ceff_pow/p(6),cc(6,8),h(8),a(8),e(8),g(8),f(8)
      common/fmass/em(3),um(3),dm(3)
      if (i.ne.7) stop 'c1_eff_l(i,bm) defined for i=7 only'
      eta = alfas(wm)/alfas(bm)
c     MS bar top quark mass at M_W
      aw = alfas(wm)
      at = alfas(um(3))
      tm_w = um(3)*(aw/at)**(12.d0/23)
      x = tm_w*tm_w/wm2
      z = dcmplx(1 - 1/x,0.d0)
      ex = - x*(18 - 11*x - x*x)/12/(x - 1)**3 - 2.d0*log(x)/3
     $     + x*x*(15 - 16*x + 4*x*x)/6/(x - 1)**4*log(x)
      c7 = - x*(16*x*x*x + 122*x*x - 80*x + 8)/9/(x - 1)**4*dble(li2(z))
     $     + x*x*(6*x*x + 46*x - 28)/3/(x - 1)**5*(log(x))**2
     $     + ( - 102*x**5 - 588*x**4 - 2262*x**3 + 3244*x*x - 1364*x 
     $     + 208)/81/(x - 1)**5*log(x) 
     $     + (1646*x**4 + 12205*x**3 - 10740*x*x + 2509*x 
     $     - 436)/486/(x - 1)**4
      c8 =  x*(- 4*x*x*x + 40*x*x + 41*x + 1)/6/(x - 1)**4*dble(li2(z))
     $     - x*x*(17*x + 31)/2/(x - 1)**5*(log(x))**2
     $     + ( - 210*x**5 + 1086*x**4 + 4893*x**3 + 2857*x*x - 1994*x 
     $     + 280)/216/(x - 1)**5*log(x) 
     $     + (737*x**4 - 14102*x**3 - 28209*x*x + 610*x 
     $     - 508)/1296/(x - 1)**4
      c1_eff_l = eta**(39/23.d0)*c7 
     $     + 8.d0/3*(eta**(37.d0/23) - eta**(39.d0/23))*c8
     $     + (297664.d0/14283*eta**(16.d0/23)
     $     - 7164416.d0/357075*eta**(14.d0/23)
     $     + 256868.d0/14283*eta**(37.d0/23)
     $     - 6698884.d0/357075*eta**(39.d0/23))*c8_bsg_l()
     $     + 37208.d0/4761*(eta**(39.d0/23) 
     $     - eta**(16.d0/23))*c7_bsg_l()
      do j=1,8
         c1_eff_l = c1_eff_l  + (eta*(e(j)*ex + g(j)) + f(j))*eta**a(j)
      end do
      return
      end

      double complex function c0_eff_r(bm)
c     leading order c_7 effective cofficients
c     non SM contribution from opposite chirality coupling
      implicit double precision (a-h,o-z)
      double complex c7_bsg_r,c8_bsg_r
      common/vpar/st,ct,st2,ct2,sct,sct2,q,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      eta = alfas(wm)/alfas(bm)
      c0_eff_r = eta**(16.d0/23)*c7_bsg_r() 
     $     + 8.d0/3*(eta**(14.d0/23) - eta**(16.d0/23))*c8_bsg_r()
      return
      end

      double complex function c1_eff_r(bm)
c     non-leading order c_7 effective cofficients at bm scale
c     non SM contribution from opposite chirality coupling
      implicit double precision (a-h,o-z)
      double complex c7_bsg_r,c8_bsg_r
      common/vpar/st,ct,st2,ct2,sct,sct2,q,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      eta = alfas(wm)/alfas(bm)
      c1_eff_r = (297664.d0/14283*eta**(16.d0/23)
     $     - 7164416.d0/357075*eta**(14.d0/23)
     $     + 256868.d0/14283*eta**(37.d0/23)
     $     - 6698884.d0/357075*eta**(39.d0/23))*c8_bsg_r()
     $     + 37208.d0/4761*(eta**(39.d0/23) 
     $     - eta**(16.d0/23))*c7_bsg_r()
      return
      end

      double complex function c7_bsg_l()
c     Coefficient c7^(0) at MW scale
      implicit double precision (a-h,o-z)
      double complex cfl(5),cfr(5)
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/fmass/em(3),um(3),dm(3)
      common/fermi/g_ferm
      do ijk=1,5
         cfl(ijk) = (0.d0,0.d0)
         cfr(ijk) = (0.d0,0.d0)
      end do
      call dd_gam(3,2,cfl,cfr)
      c7_bsg_l = - 8*pi*pi/sq2/g_ferm/dm(3)/ckm_phys(3,2)
     $     /dconjg(ckm_phys(3,3))*(cfl(1) - cfl(4)*dm(3) + cfr(4)*dm(2))
      return
      end

      double complex function c8_bsg_l()
c     Coefficient c8^(0) at MW scale
      implicit double precision (a-h,o-z)
      double complex cfl(5),cfr(5)
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/fmass/em(3),um(3),dm(3)
      common/fermi/g_ferm
      do ijk=1,5
         cfl(ijk) = (0.d0,0.d0)
         cfr(ijk) = (0.d0,0.d0)
      end do
      call dd_gluon(3,2,cfl,cfr)
      c8_bsg_l = - 8*pi*pi/sq2/g_ferm/dm(3)/ckm_phys(3,2)
     $     /dconjg(ckm_phys(3,3))*(cfl(1) - cfl(4)*dm(3) + cfr(4)*dm(2))
      return
      end

      double complex function c7_bsg_r()
c     Coefficient c7^(0) at MW scale
      implicit double precision (a-h,o-z)
      double complex cfl(5),cfr(5)
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/fmass/em(3),um(3),dm(3)
      common/fermi/g_ferm
      do ijk=1,5
         cfl(ijk) = (0.d0,0.d0)
         cfr(ijk) = (0.d0,0.d0)
      end do
      call dd_gam(3,2,cfl,cfr)
      c7_bsg_r = - 8*pi*pi/sq2/g_ferm/dm(3)/ckm_phys(3,2)
     $     /dconjg(ckm_phys(3,3))*(cfr(1) - cfr(4)*dm(3) + cfl(4)*dm(2))
      return
      end

      double complex function c8_bsg_r()
c     Coefficient c8^(0) at MW scale
      implicit double precision (a-h,o-z)
      double complex cfl(5),cfr(5)
      double complex ckm_phys,ckm0,udl,udr,uul,uur
      common/ckm_switch/ckm_phys(3,3),ckm0(3,3),udl(3,3),udr(3,3),
     $     uul(3,3),uur(3,3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/fmass/em(3),um(3),dm(3)
      common/fermi/g_ferm
      do ijk=1,5
         cfl(ijk) = (0.d0,0.d0)
         cfr(ijk) = (0.d0,0.d0)
      end do
      call dd_gluon(3,2,cfl,cfr)
      c8_bsg_r = - 8*pi*pi/sq2/g_ferm/dm(3)/ckm_phys(3,2)
     $     /dconjg(ckm_phys(3,3))*(cfr(1) - cfr(4)*dm(3) + cfl(4)*dm(2))
      return
      end

      subroutine fm_bsg_init(d)
c     Photon cutoff coefficients - initialization
      implicit double precision (a-h,o-z)
      double precision gia_bsg,gir_bsg
      double complex li2,li2arg     
      logical fm_init
      common/vpar/st,ct,st2,ct2,sct,sct2,q,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/qmass_pole/um(3),dm(3)
      common/bsg_auxpar/z,ip
      common/bsg_fm/fmc(8,8),d_old,fm_init
      external gia_bsg,gir_bsg
c     Basic phase space integrals for three-loop bsgamma
      z = um(2)*um(2)/dm(3)/dm(3)
      a = (1 - d)/z
      b = 1/z
      errin = 1.d-3
      niter = 20
      ip = 0
      r0 = romberg(gir_bsg,0.d0,a,errin,err1,niter)
      ip = 1
      a1 = romberg(gia_bsg,0.d0,a,errin,err2,niter)
      r1 = romberg(gir_bsg,a,b,errin,err3,niter)
      ip = 2
      a2 = romberg(gia_bsg,a,b,errin,err4,niter)
      err = max(err1,err2,err3,err4)
      if (err/errin.gt.10) 
     $     stop 'Unexpectedly large error of integration in fm_bsg_init'
c     Common dilogarithm and logarithms
      li2arg = 1 - d
      dil = dble(li2(li2arg))
      dl1 = log(d)
      dl2 = log(1 - d)
c     f_ij functions initialization
c     Results stored in common/bsg_fm/fmc(8,8)
      do i=1,8
         do j=1,8
            fmc(i,j) = 0
         end do
      end do
      fmc(8,8) = ( - 2*log(dm(3)/dm(2))*(d*d + 2*d + 4*dl2) + 4*dil
     $   - 2.d0*pi*pi/3 - d*(2 + d)*dl1 + 8*dl2 - 2.d0*d*d*d/3 + 3*d*d
     $   + 7*d)/27
      fmc(7,7) = d*(10 + d - 2.d0*d*d/3 + (d - 4)*dl1)/3
      fmc(7,8) = 2.d0/9*(4*dil - 2.d0*pi*pi/3 - 4*d*dl1 + 9*d - d*d
     $     + d*d*d/3)
      fmc(2,2) = 16.d0*z/27*(d*a1 + a2)
      fmc(2,7) = - 8.d0*z*z/9*(d*r0 + r1)
      fmc(2,8) = - fmc(2,7)/3
      fmc(1,1) =   fmc(2,2)/36
      fmc(1,2) = - fmc(2,2)/3
      fmc(1,7) = - fmc(2,7)/6
      fmc(1,8) = - fmc(2,8)/6
      d_old = d
      fm_init = .false.
      return
      end

      double complex function gi_bsg(t)
      implicit double precision (a-h,o-z)
      common/vpar/st,ct,st2,ct2,sct,sct2,q,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      if (t.lt.0) stop 't<0 in gi_bsg'
      if (t.lt.4) then
         cl = atan(sqrt(t/(4 - t)))
         gi_bsg = - dcmplx(2*cl*cl,0.d0)
      else
         cl = log((sqrt(t) + sqrt(t - 4))/2)
         gi_bsg = 2*dcmplx(cl*cl - pi*pi/4, - pi*cl) 
      end if
      return
      end
      
      double precision function gia_bsg(t)
      implicit double precision (a-h,o-z)
      double complex gi_bsg
      common/bsg_auxpar/z,ip
      if (t.lt.0) stop 't<=0 in gia_bsg'
      if (t.eq.0) then
         gia_bsg = 0
      else
         gia_bsg = abs(1 - z*t)**ip*abs(gi_bsg(t)/t + 0.5d0)**2
      end if
      return
      end

      double precision function gir_bsg(t)
      implicit double precision (a-h,o-z)
      double complex gi_bsg
      common/bsg_auxpar/z,ip
      if (t.lt.0) stop 't<0 in gir_bsg'
      gir_bsg = abs(1 - z*t)**ip*dble(gi_bsg(t) + t/2)
      return
      end

      subroutine bsg_init
      implicit double precision (a-h,o-z)
      double complex r
      logical bxg_init
      dimension p0(2),a0(4)
      common/vpar/st,ct,st2,ct2,sct,sct2,q,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/bsg_dat/al1,al2,br_cev,al_em,xi3,r(8),gam(8),bxg_init
      common/ceff_pow/p(6),cc(6,8),h(8),a(8),e(8),g(8),f(8)
      common/qmass_pole/um(3),dm(3)
c     ceff_pow initialization
      data p0/-2.d0,1.d0/
      data a0/7.d0,8.d0,3.d0,-6.d0/
c     c_eff coefficients
      do i=1,2
         p(i) = 6.d0*p0(i)/23
      end do
      cc(1,2) = 1.d0/3
      cc(2,2) = 2.d0/3
      cc(1,3) = - 1.d0/27
      cc(2,3) = 2.d0/63
      cc(1,4) = 1.d0/9
      cc(2,4) = 1.d0/21
      cc(1,5) = 1.d0/108
      cc(2,5) = - 1.d0/126
      cc(1,6) = - 1.d0/36
      cc(2,6) = - 1.d0/84
      do i=1,4
        a(i) = 2.d0*a0(i)/23
      end do
      h(1) = 626126.d0/272277
      h(2) = - 56281.d0/51730
      h(3) = - 3.d0/7
      h(4) = - 1.d0/14
      e(1) = 4661194.d0/816831
      e(2) =  - 8516.d0/2217
c     r_i vector
      z = um(2)*um(2)/dm(3)/dm(3)
      zh = sqrt(z)
      z2 = z*z
      zl = log(z)
      zl2 = zl*zl
      pi2 = pi*pi
      r(2) = 2.d0/243*dcmplx(- 833 + 144*pi2*z*zh 
     $     + z*(1728 - 180*pi2 - 1296*xi3 + (1296 - 324*pi2)*zl
     $     + 108*zl2 + 36*zl*zl2)
     $     + 36*z2*(18 + 2*pi2 + (12 - 6*pi2)*zl + zl*zl2)
     $     + 6*z*z2*(- 9 - 14*pi2 + 182*zl - 126*zl2),
     $     24*pi*(- 5 + 3*z*(15 - pi2 + 3*zl + 3*zl2)
     $     + 3*z2*(- pi2 + 3*zl2) + 4*z*z2*(7 - 3*zl)))
      r(1) = - r(2)/6
      r(3) = (0.d0,0.d0)
      r(4) = (0.d0,0.d0)
      r(5) = (0.d0,0.d0)
      r(6) = (0.d0,0.d0)
      r(7) = - 2.d0/9*(15 + 4*pi2)
      r(8) = - 4.d0/27*dcmplx(- 33 + 2*pi2, - 6*pi)
c     gamma^(0)eff_i7 vector
      gam(1) = - 208/243.d0
      gam(2) =   416/81.d0
      gam(3) = - 176/81.d0
      gam(4) = - 152/243.d0
      gam(5) = - 6272/81.d0
      gam(6) =   4624/243.d0
      gam(7) =   32/3.d0
      gam(8) = - 32/9.d0
      bxg_init = .true.
      return
      end

      block data init_bsg
      implicit double precision (a-h,o-z)
      double complex r
      logical bxg_init,fm_init
      common/bsg_fm/fmc(8,8),d_old,fm_init
      common/bsg_dat/al1,al2,br_cev,al_em,xi3,r(8),gam(8),bxg_init
      common/ceff_pow/p(6),cc(6,8),h(8),a(8),e(8),g(8),f(8)
      data al1,al2,br_cev,al_em,xi3/0.5d0,0.122d0,0.108d0,0.00729735d0,
     $     1.20206d0/
      data bxg_init,fm_init/.false.,.true./
      data d_old/-1.d0/
c     ceff_pow initialization
      data p/0.d0,0.d0,0.4086d0,-0.4230d0,-0.8994d0,0.1456d0/
      data cc/-1.d0,1.d0,0.d0,0.d0,0.d0,0.d0,
     $        6*0.d0,
     $        2*0.d0,-0.0659d0, 0.0595d0,-0.0218d0, 0.0335d0,
     $        2*0.d0, 0.0237d0,-0.0173d0,-0.1336d0,-0.0316d0,
     $        2*0.d0, 0.0094d0,-0.0100d0, 0.0010d0,-0.0017d0,
     $        2*0.d0, 0.0108d0, 0.0163d0, 0.0103d0, 0.0023d0,
     $        6*0.d0,
     $        2*0.d0,-0.9135d0, 0.0873d0,-0.0571d0, 0.0209d0/
      data a/4*0.d0, 0.4086d0,-0.4230d0,-0.8994d0, 0.1456d0/
      data h/4*0.d0,-0.6494d0,-0.0380d0,-0.0186d0,-0.0057d0/
      data e/4*0.d0,-1.9043d0,-0.1008d0,0.1216d0,0.0183d0/
      data f/-17.3023d0,8.5027d0,4.5508d0,0.7519d0,2.0040d0,
     $     0.7476d0,-0.5385d0,0.0914d0/
      data g/14.8088d0,-10.8090d0,-0.8740d0,0.4218d0,-2.9347d0,
     $     0.3971d0,0.1600d0,0.0225d0/
      end
      


