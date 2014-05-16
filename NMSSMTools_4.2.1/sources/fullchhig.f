      subroutine getPiHpHm(g,gp,ll,kk,ht,hb,htau,v1,v2,xx,
     $     Al,At,Ab,Atau,p,Q,piHpHm)
* by courtesy of P. Slavich and K.H. Phan

      implicit none

      double precision g,gp,ll,kk,ht,hb,htau,v1,v2,xx,
     $     Al,At,Ab,Atau,Q,p,piHpHm
      
      double precision mstop(2),msbot(2),mstau(2),msup(2),msdown(2),
     $     msel(2),msnutau,msnue,Rt(2,2),Rb(2,2),Rtau(2,2),
     $     NN(5,5),UU(2,2),VV(2,2),mch(2),mne(5),
     $     mhh(3),maa(3),mhc(2),RS(3,3),RP(3,3),RC(2,2),
     $     myB0,myF,myA0,myG,pi,cb,sb,mw2,mz2,mt2,mb2,mtau2,
     $     gb2,sw2,cw2,aHpnech(5,2),bHpnech(5,2),
     $     lHpHmtt(2,2),lHpHmbb(2,2),lHpHmtata(2,2),lHpHmntnt,
     $     lHpHmuu(2,2),lHpHmdd(2,2),lHpHmee(2,2),lHpHmnn,
     $     lHptb(2,2),lHpntta(2),lHpud(2,2),lHpne(2),  
     $     lHpHmhh(3,3),lHpHmaa(3,3),lHphc(3,2),lHpac(3,2),lHpHmcc(2),
     $     higgs,weak,fermions,sfermions,sfermions3g,inos

      integer i,j
      
      common/TREE_MASSES/mhh,maa,mhc,RS,RP,RC,mch,UU,VV,mne,NN,
     $     mstop,msbot,mstau,Rt,Rb,Rtau,msup,msdown,msel,msnutau,msnue
      
      pi = 4d0*atan(1.d0)

      PiHpHm = 0d0

      cb = v1/sqrt(v1**2+v2**2)
      sb = v2/sqrt(v1**2+v2**2)
      gb2 = (g**2+gp**2)/2d0
      sw2 = gp**2/2d0/gb2
      cw2 = 1-sw2
      mw2 = g**2/2d0*(v1**2+v2**2)
      mz2 = gb2*(v1**2+v2**2)

*     the higgs+gauge and pure gauge contributions

      weak = 0d0

      do i = 1,3                

         weak = weak            ! scalar-W
     $        + g**2/4d0*(sb*RS(i,1)-cb*RS(i,2))**2
     $        *myF(p**2,mhh(i)**2,mw2,Q**2)
         
         weak = weak            ! pseudoscalar-W
     $        + g**2/4d0*(sb*RP(i,1)+cb*RP(i,2))**2
     $        *myF(p**2,maa(i)**2,mw2,Q**2)
         
      enddo

      weak = weak               ! pure gauge
     $     +g**2/4d0*(cw2-sw2)**2/cw2*myF(p**2,mhc(2)**2,mz2,Q**2)
     $     +sw2*g**2*myF(p**2,mhc(2)**2,0d0,Q**2)
     $     +2*g**2*myA0(mw2,Q**2)
     $     +g**2*(cw2-sw2)**2/cw2*myA0(mz2,Q**2)
      
*     the pure higgs contributions

      call coupl_Hp_hh(g,gp,ll,kk,v1,v2,xx,Al,RS,RP,
     $     lHpHmhh,lHpHmaa,lHphc,lHpac,lHpHmcc)

      higgs = 0d0

      do i = 1,3
         do j=1,2
            
            higgs = higgs       ! charged+scalar and charged+pseudo
     $           +lHphc(i,j)**2*myB0(p**2,mhh(i)**2,mhc(j)**2,Q**2)
     $           +lHpac(i,j)**2*myB0(p**2,maa(i)**2,mhc(j)**2,Q**2)
         enddo

         higgs = higgs          ! scalar and pseudoscalar bubbles
     $        +lHpHmhh(i,i)*myA0(mhh(i)**2,Q**2)
     $        +lHpHmaa(i,i)*myA0(maa(i)**2,Q**2)
      enddo

      higgs = higgs             ! charged bubbles (note the symmetry factor)
     $     +lHpHmcc(1)*myA0(mhc(1)**2,Q**2)
     $     +4*lHpHmcc(2)*myA0(mhc(2)**2,Q**2)
               
*     the fermion contributions

      mt2 = ht**2*v2**2
      mb2 = hb**2*v1**2
      mtau2 = htau**2*v1**2
      
      fermions = 
     $     +3*((ht**2*cb**2+hb**2*sb**2)*myG(p**2,mt2,mb2,Q**2)
     $     -2*ht*hb*Sqrt(mt2*mb2)*2*sb*cb*myB0(p**2,mt2,mb2,Q**2))
     $     +htau**2*sb**2*myG(p**2,mtau2,0d0,Q**2)

*     the sfermion contributions

      call coupl_Hp_sf(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
     $     Rt,Rb,Rtau,lHpHmtt,lHpHmbb,lHpHmtata,lHpHmntnt,lHpHmuu,
     $     lHpHmdd,lHpHmee,lHpHmnn,lHptb,lHpntta,lHpud,lHpne)

      sfermions = 0d0           ! first two generations
      
      do i=1,2
         do j = 1,2

            sfermions = sfermions ! up-down squarks
     $          +6*lHpud(i,j)**2*myB0(p**2,msup(i)**2,msdown(j)**2,Q**2)
         enddo
         
         sfermions = sfermions  ! up-down sleptons
     $        +2*lHpne(i)**2*myB0(p**2,msel(i)**2,msnue**2,Q**2)
         
         sfermions = sfermions  ! charged-sfermion bubble
     $        +6*lHpHmuu(i,i)*myA0(msup(i)**2,Q**2)
     $        +6*lHpHmdd(i,i)*myA0(msdown(i)**2,Q**2)
     $        +2*lHpHmee(i,i)*myA0(msel(i)**2,Q**2)

      enddo
      
      sfermions = sfermions     ! sneutrino bubble
     $     + 2*lHpHmnn*myA0(msnue**2,Q**2)
      
      sfermions3g = 0d0         ! third generation

      do i=1,2
         do j = 1,2

            sfermions3g = sfermions3g ! top-bottom squarks
     $          +3*lHptb(i,j)**2*myB0(p**2,mstop(i)**2,msbot(j)**2,Q**2)
         enddo
         
         sfermions3g = sfermions3g ! stau-sneutrino
     $        +lHpntta(i)**2*myB0(p**2,mstau(i)**2,msnutau**2,Q**2)
         
         sfermions3g = sfermions3g ! charged-sfermion bubble
     $        +3*lHpHmtt(i,i)*myA0(mstop(i)**2,Q**2)
     $        +3*lHpHmbb(i,i)*myA0(msbot(i)**2,Q**2)
     $        +lHpHmtata(i,i)*myA0(mstau(i)**2,Q**2)
         
      enddo

      sfermions3g = sfermions3g           ! sneutrino bubble       
     $     +lHpHmntnt*myA0(msnutau**2,Q**2)
      
*     the chargino/neutralino contribution

      inos = 0d0

      call coupl_Hp_ino(g,gp,ll,v1,v2,NN,UU,VV,aHpnech,bHpnech)

      do i = 1,5
         do j = 1,2
            
            inos = inos
     $           +(aHpnech(i,j)**2+bHpnech(i,j)**2)
     $           *myG(p**2,mne(i)**2,mch(j)**2,Q**2)
     $           -4*aHpnech(i,j)*bHpnech(i,j)*
     $           mne(i)*mch(j)*myB0(p**2,mne(i)**2,mch(j)**2,Q**2)
         enddo
      enddo

      PiHpHm = weak + higgs + fermions + sfermions + sfermions3g + inos

c      write(*,*) 'gauge =',weak
c      write(*,*) 'higgs =',higgs
c      write(*,*) 'fermions =',fermions
c      write(*,*) 'sfermions1+2g =',sfermions
c      write(*,*) 'sfermions3g =',sfermions3g
c      write(*,*) 'inos =',inos
c      write(*,*) 'total =',PiHpHm

      PiHpHm = PiHpHm/16d0/pi**2

      end
      
*      
***********************************************************************
*     

      subroutine coupl_Hp_ino(g,gp,ll,v1,v2,NN,UU,VV,aHpnech,bHpnech)
      
      implicit none

      double precision g,gp,ll,v1,v2,NN(5,5),UU(2,2),VV(2,2),
     $     aHpnech(5,2),bHpnech(5,2)

      double precision cb,sb,sq2

      integer i,j
      
      cb = v1/sqrt(v1**2+v2**2)
      sb = v2/v1*cb

      sq2 = sqrt(2d0)

      do i=1,5
         do j=1,2
            aHpnech(i,j) = 0d0
            bHpnech(i,j) = 0d0
         enddo
      enddo

      do i=1,5
         do j=1,2
            aHpnech(i,j) = aHpnech(i,j)
     $           +sb*(-g*NN(i,3)*UU(j,1)
     $           +(g*NN(i,2)+gp*NN(i,1))*UU(j,2)/sq2)
     $           -cb*ll*NN(i,5)*UU(j,2)
            
            bHpnech(i,j) = bHpnech(i,j)
     $           +cb*(-g*NN(i,4)*VV(j,1)
     $           -(g*NN(i,2)+gp*NN(i,1))*VV(j,2)/sq2)
     $           -sb*ll*NN(i,5)*VV(j,2)
         enddo
      enddo
         
      end

*      
***********************************************************************
*     

      subroutine coupl_Hp_sf(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
     $     Rt,Rb,Rtau,lHpHmtt,lHpHmbb,lHpHmtata,lHpHmntnt,lHpHmuu,
     $     lHpHmdd,lHpHmee,lHpHmnn,lHptb,lHpntta,lHpud,lHpne)

      implicit none

      double precision g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
     $     Rt(2,2),Rb(2,2),Rtau(2,2),lHpHmtt(2,2),lHpHmbb(2,2),
     $     lHpHmtata(2,2),lHpHmntnt,lHpHmuu(2,2),lHpHmdd(2,2),
     $     lHpHmee(2,2),lHpHmnn,lHptb(2,2),lHpntta(2),lHpud(2,2),
     $     lHpne(2),lHpHmtt_int(2,2),lHpHmbb_int(2,2),
     $     lHpHmtata_int(2,2),lHptb_int(2,2),lHpntta_int(2)

      double precision YuL,YuR,YdL,YdR,YeL,YeR,Ynu,cb,sb,c2b

      integer i,j,k,l

      cb = v1/sqrt(v1**2+v2**2)
      sb = v2/v1*cb
      c2b = cb**2-sb**2

      YuL = 1/3d0 
      YdL = 1/3d0
      Ynu = -1d0
      YeL = -1d0
      YuR = -4/3d0
      YdR = 2/3d0
      YeR = 2d0

*     FIRST TWO GENERATIONS

*     quartic

      do i=1,2
         do j=1,2
            lHpHmuu(i,j) = 0d0
            lHpHmdd(i,j) = 0d0
            lHpHmee(i,j) = 0d0
         enddo
      enddo

      lHpHmuu(1,1) = c2b/4d0*(g**2+YuL*gp**2)
      lHpHmuu(2,2) = c2b/4d0*YuR*gp**2
      
      lHpHmdd(1,1) = -c2b/4d0*(g**2-YdL*gp**2)
      lHpHmdd(2,2) = c2b/4d0*YdR*gp**2
      
      lHpHmee(1,1) = -c2b/4d0*(g**2-YeL*gp**2)
      lHpHmee(2,2) = c2b/4d0*YeR*gp**2

      lHpHmnn = c2b/4d0*(g**2+Ynu*gp**2)
      
*     trilinear

      do i=1,2
         do j=1,2
            lHpud(i,j) = 0d0
         enddo
            lHpne(i) = 0d0            
      enddo

      lHpud(1,1) = g**2/2d0*(v2*cb+v1*sb)
      
      lHpne(1) = g**2/2d0*(v2*cb+v1*sb)

*     THIRD GENERATION (add the Yukawa interactions)

*     quartic

      do i=1,2
         do j=1,2
            lHpHmtt_int(i,j) = lHpHmuu(i,j)
            lHpHmbb_int(i,j) = lHpHmdd(i,j)
            lHpHmtata_int(i,j) = lHpHmee(i,j)
         enddo
      enddo

      lHpHmtt_int(1,1) = lHpHmtt_int(1,1) + hb**2*sb**2
      lHpHmtt_int(2,2) = lHpHmtt_int(2,2) + ht**2*cb**2

      lHpHmbb_int(1,1) = lHpHmbb_int(1,1) + ht**2*cb**2
      lHpHmbb_int(2,2) = lHpHmbb_int(2,2) + hb**2*sb**2

      lHpHmtata_int(2,2) = lHpHmtata_int(2,2) + htau**2*sb**2

      lHpHmntnt = lHpHmnn + htau**2*sb**2

*     trilinear

      lHptb_int(1,1) = lHpud(1,1) - ht**2*v2*cb - hb**2*v1*sb
      lHptb_int(2,2) = lHpud(2,2) - ht*hb*(v1*cb + v2*sb)
      lHptb_int(1,2) = lHpud(1,2) - hb*(ll*xx*cb + Ab*sb)
      lHptb_int(2,1) = lHpud(2,1) - ht*(ll*xx*sb + At*cb)

      lHpntta_int(1) = lHpne(1) - htau**2*v1*sb
      lHpntta_int(2) = lHpne(2) - htau*(ll*xx*cb + Atau*sb)
      
*     now rotate the third-generation couplings

      do i = 1,2
         do j = 1,2
            lHpHmtt(i,j) = 0d0
            lHpHmbb(i,j) = 0d0
            lHpHmtata(i,j) = 0d0
            lHptb(i,j) = 0d0
            do k=1,2
               do l=1,2
                  lHpHmtt(i,j) = lHpHmtt(i,j)
     $                 +Rt(i,k)*Rt(j,l)*lHpHmtt_int(k,l)
                  lHpHmbb(i,j) = lHpHmbb(i,j)
     $                 +Rb(i,k)*Rb(j,l)*lHpHmbb_int(k,l)
                  lHpHmtata(i,j) = lHpHmtata(i,j)
     $                 +Rtau(i,k)*Rtau(j,l)*lHpHmtata_int(k,l)
                  lHptb(i,j) = lHptb(i,j)
     $                 +Rt(i,k)*Rb(j,l)*lHptb_int(k,l)
               enddo
            enddo
         enddo
      enddo
      
*     stau-sneutrino trilinear

      do i = 1,2
         lHpntta(i) = 0d0
         do j=1,2
            lHpntta(i) = lHpntta(i) + Rtau(i,j)*lHpntta_int(j)
         enddo
      enddo

      end

*      
***********************************************************************
*     

      subroutine coupl_Hp_hh(g,gp,ll,kk,v1,v2,xx,Al,RS,RP,
     $     lHpHmhh,lHpHmaa,lHphc,lHpac,lHpHmcc)

      implicit none

      double precision g,gp,ll,kk,v1,v2,xx,Al,RS(3,3),RP(3,3),
     $     lHpHmhh(3,3),lHpHmaa(3,3),lHphc(3,2),lHpac(3,2),lHpHmcc(2)

      double precision sq2,c2b,s2b,
     $     lHpHmss(3,3),lHpHmpp(3,3),lHpsc(3,2),lHppc(3,2)

      integer i,j,k,l
      
      c2b = (v1**2-v2**2)/(v1**2+v2**2)
      s2b = 2*v1*v2/(v1**2+v2**2)
      sq2 = sqrt(2d0)

      do i=1,3                  ! initialize
         do j=1,3
            lHpHmhh(i,j) = 0d0
            lHpHmaa(i,j) = 0d0
            lHpHmss(i,j) = 0d0
            lHpHmpp(i,j) = 0d0
            if(j.lt.3) then
               lHphc(i,j) = 0d0
               lHpac(i,j) = 0d0
               lHpsc(i,j) = 0d0
               lHppc(i,j) = 0d0   
            endif
         enddo
      enddo

*     quartic charged-neutral couplings
 
      lHpHmss(1,1) = (g**2-gp**2*c2b)/8d0
      lHpHmss(1,2) = -(2*ll**2-g**2)*s2b/8d0
      lHpHmss(2,1) = lHpHmss(1,2)
      lHpHmss(2,2) = (g**2+gp**2*c2b)/8d0 
      lHpHmss(3,3) = ll*(ll+kk*s2b)/2d0

      lHpHmpp(1,1) = (g**2-gp**2*c2b)/8d0
      lHpHmpp(1,2) = (2*ll**2-g**2)*s2b/8d0
      lHpHmpp(2,1) = lHpHmpp(1,2)
      lHpHmpp(2,2) = (g**2+gp**2*c2b)/8d0 
      lHpHmpp(3,3) = ll*(ll-kk*s2b)/2d0
      
*     rotate the quartics

      do i=1,3              
         do j=1,3
            do k=1,3
               do l=1,3
                  lHpHmhh(i,j) = lHpHmhh(i,j)
     $                 +RS(i,k)*RS(j,l)*lHpHmss(k,l)
                  lHpHmaa(i,j) = lHpHmaa(i,j)
     $                 +RP(i,k)*RP(j,l)*lHpHmpp(k,l)
               enddo
            enddo
         enddo
      enddo
      
*     trilinear charged-neutral couplings

      lHpsc(1,1) = (-v1*gp**2*s2b+v2*(2*ll**2-g**2)*c2b)/2d0/sq2
      lHpsc(1,2) = (v1*(g**2-gp**2*c2b)-v2*(2*ll**2-g**2)*s2b)/2d0/sq2
      lHpsc(2,1) = (v2*gp**2*s2b+v1*(2*ll**2-g**2)*c2b)/2d0/sq2
      lHpsc(2,2) = (v2*(g**2+gp**2*c2b)-v1*(2*ll**2-g**2)*s2b)/2d0/sq2
      lHpsc(3,1) = -ll/sq2*(Al+2*kk*xx)*c2b
      lHpsc(3,2) = ll/sq2*(2*ll*xx+(Al+2*kk*xx)*s2b)

      lHppc(1,1) = v2*(2*ll**2-g**2)/2d0/sq2               
      lHppc(2,1) = v1*(2*ll**2-g**2)/2d0/sq2               
      lHppc(3,1) = ll/sq2*(Al-2*kk*xx)

*     rotate the trilinears

      do i=1,3              
         do j=1,2
            do k=1,3
               lHphc(i,j) = lHphc(i,j) + RS(i,k)*lHpsc(k,j)
               lHpac(i,j) = lHpac(i,j) + RP(i,k)*lHppc(k,j)
            enddo
         enddo
      enddo

*     quartic charged-charged couplings

c      lHpHmcc(1) = -(g**2+gp**2)*(c2b**2-s2b**2)/4d0
c      lHpHmcc(2) = (g**2+gp**2)*c2b**2/8d0   
      lHpHmcc(1) = -(g**2+gp**2)*(c2b**2-s2b**2)/4d0 + ll**2*c2b**2
      lHpHmcc(2) = (g**2+gp**2)*c2b**2/8d0 + ll**2/4d0*s2b**2  

      end

