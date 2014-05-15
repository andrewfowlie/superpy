* This file contains
*
* the main routine FULLHIG for the full 1- and 2-Loop corrections
*   from Degrassi/Slavich,
*   ``On the radiative corrections to the neutral Higgs boson masses in
*    the NMSSM,'' Nucl. Phys.  B {\bf 825} (2010) 119
*   [arXiv:0907.4682 [hep-ph]].
*
* subroutine mysort(msd,ZS)
* subroutine getdVB(g,gp,mzpole,mwpole,Q,dVB)
* subroutine effpot(lp,mt,mg,T1,T2,st,ct,q2,tanb,vv,l,xx,as,
*     .     DMS,DMP)
* subroutine makefuncs(mt,mg,T1,T2,s2t,c2t,q,tanb,At,mu,
*     .     F1t,F2t,F3t,Ft,FA)
* subroutine makederiv(mt,mg,T1,T2,s2t,c2t,q,
*     .     DT1,DT1T1,DT1t,DT1c2t,DT1T2,Dtt,Dc2t,Dc2tc2t,Dtc2t,Dcptmptt)
* subroutine diagonalize(mq,mw,mz,msql,msqr,Aq,mu,tb,iq,
*     .     msq2,sth,cth)
* subroutine squarks(mtsm,mbsm,mz,mw,mqs,mts,mbs,At,Ab,mu,tb,mg,as,
*     .     A0,vv,q2,mtmssm,mbmssm,mstop2,msbot2,cst,sst,csb,ssb,
*     .     errsqua,mtpole,asma,ast,v2)
* subroutine getdmt(mt,mstop2,sst,cst,mb,msbot2,ssb,csb,
*     .      vv,mg,A0,mu,tb,as,q2,dmt)
* subroutine getdmb(mt,mstop2,sst,cst,mb,msbot2,ssb,csb,
*     .      vv,mg,A0,mu,tb,as,q2,dd,eb)
* subroutine swap12(M)
* subroutine jacobi(a,np,d,v)
* subroutine HEigensystem(n, A,ldA, d, U,ldU, sort)
*****************************************************************
* The following subroutines are in "fullhig1.f":
*
* subroutine gettadS(g,gp,ht,hb,htau,v1,v2,Q,tadS)
* subroutine getPiSS(g,gp,ht,hb,htau,v1,v2,p,Q,piSS)
* subroutine getPiPP(g,gp,ht,hb,htau,v1,v2,p,Q,piPP)
* subroutine getPiZZ(g,gp,ht,hb,htau,v1,v2,p,Q,piZZ)
* subroutine getPiWW(g,gp,ht,hb,htau,v1,v2,p,Q,piWW)
* subroutine treemasses(g,gp,ll,kk,ht,hb,htau,v1,v2,xx,M1,M2,
*     .     Ak,Al,At,Ab,Atau,mQ3,mtr,mbr,mQ,mur,mdr,mL3,mtaur,mL,mer,
*     .     Q,errmass)
*     defining the
*   common/TREE_MASSES/mhh,maa,mhc,RS,RP,RC,mch,UU,VV,mne,NN,
*     .     mstop,msbot,mstau,Rt,Rb,Rtau,msup,msdown,msel,msnutau,msnue
* subroutine tree_charginos(g,ll,v1,v2,M2,xx,xmc,u,v)
* subroutine tree_neutralinos(g,gp,ll,kk,v1,v2,xx,M1,M2,xmn,Z)
* subroutine tree_higgses(g,gp,ll,kk,v1,v2,xx,Ak,Al,
*     .     mss,maa,mhc,RS,RP,RC,errhiggs)
* subroutine tree_sfermions(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
*     .     mQ3,mtr,mbr,mQ,mur,mdr,mL3,mtaur,mL,mer,mstop,msbot,
*     .     mstau,Rt,Rb,Rtau,msup,msdown,msel,msnutau,msnue,errsfer)
* subroutine diagsfe(n,hf,g,gp,ll,v1,v2,xx,Af,mL,mR,mass,R,error)
* subroutine scalarcouplings(g,gp,ll,kk,ht,hb,htau,v1,v2,xx,
*     .    Al,Ak,At,Ab,Atau)
* subroutine coupl_s_sf(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
*     .     Rt,Rb,Rtau,lsstt,lssbb,lsstata,lssntnt,lssuu,lssdd,
*     .     lssee,lssnn,lstt,lsbb,lstata,lsntnt,lsuu,lsdd,lsee,lsnn)
* subroutine coupl_s_hh(g,gp,ll,kk,v1,v2,xx,Al,Ak,RS,RP,
*     .     lsshh,lssaa,lshh,lsaa,lsscc,lscc)
* subroutine coupl_s_ino(g,gp,ll,kk,NN,UU,VV,lsnene,lschch)
* subroutine coupl_p_sf(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
*     .     Rt,Rb,Rtau,lpptt,lppbb,lpptata,lppntnt,lppuu,lppdd,
*     .     lppee,lppnn,lptt,lpbb,lptata,lpntnt,lpuu,lpdd,lpee,lpnn)
* subroutine coupl_p_hh(g,gp,ll,kk,v1,v2,xx,Al,Ak,RS,RP,
*     .     lpphh,lppaa,lpah,lppcc,lpcc)
* subroutine coupl_p_ino(g,gp,ll,kk,NN,UU,VV,lpnene,lpchch)
* subroutine coupl_Z_ino(g,gp,NN,UU,VV,lznene,azchch,bzchch)
* subroutine coupl_W_ino(g,NN,UU,VV,awnech,bwnech)
*****************************************************************
* The following subroutines are in "fullhig2.f":
*
* subroutine twlpyuk(mt,mb,A0,T1,T2,B1,B2,st,ct,sb,cb,q,l,xx,tanb,
*     .     vv,DMS,DMP)
* subroutine makefuncstb(t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,
*     .     q,mu,vv,tanb,F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,
*     .     Ft,Fb,FA)
* subroutine makederivtb(t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,
*     .     q,mu,vv,tanb)
*****************************************************************


      subroutine FULLHIG(mom,l,k,mu,tanb,mq3,mtr,mbr,
     .     mQ,mur,mdr,mL3,mtaur,mL,mer,At,Ab,Atau,Al,Ak,M1,M2,mg,Q,
     .     mss,maa,OS,OP,mhc,err)

      implicit none

      logical mom,errsqua
      integer i,j,kk
      integer err,errmass,IL
      double precision l,k,mu,tanb,mq3,mtr,mbr,mQ,mur,mdr,mL3,mtaur,mL,
     .     mer,At,Ab,Atau,Al,Ak,M1,M2,mg,q,mss(3),maa(3),OS(3,3),
     .     OP(3,3),mhc
      double precision g,gp,gb2,vv,cb,sb,v1,v2,slasf,asq,mtsm,mbsm,
     .     runt,runb,pi,mtmssm,mbmssm,mstop2(2),msbot2(2),cst,sst,
     .     csb,ssb,As,MStr(3,3),MPtr(3,3),DMS_2lt(3,3),
     .     DMP_2lt(3,3),DMS_2lb(3,3),DMP_2lb(3,3),MS(3,3),MP(3,3),
     .     ms2(3),ma2(3),xx,tadS(3),piSS(3,3),piPP(3,3),
     .     gold,ms2cor(3),ma2cor(3),
     .     MSfull(3,3),MPfull(3,3),pa,ps,vev(3),mz,mw,
     .     piZZ_MZ,piWW_MW,piWW_0,dVB,sq2,T1,T2,B1,B2,
     .     A0,DMS_2ly(3,3),DMP_2ly(3,3)
      double precision OS0(3,3),OP0(3,3),asma,ast
      double precision mZpole,mWpole,mtpole,mbmb,mtau,asmz,GF
* For charged Higgs:
      double precision mhcsqtree, mhcsqcorr,phc,piHpHm

      common/SMINPUTS/mZpole,mWpole,mtpole,mbmb,mtau,asmz,GF

*      write(*,*) 'inputs_1',mom,l,k,mu,tanb,mq3,mtr,mbr,
*     .     mQ,mur,mdr,mL3,mtaur,mL,mer,At,Ab,Atau,Al,Ak,M1,M2,mg,Q
*      write(*,*) 'inputs_2',mZpole,mWpole,mtpole,mbmb,mtau,asmz,GF

      err=0

*     some preliminary quantities

      pi = 4d0*atan(1d0)
      sq2 = sqrt(2d0)
      cb = 1d0/sqrt(1d0+tanb**2)
      sb = tanb*cb

      asq = slasf(q,mtpole,asmz,mzpole)

      mtsm = runt(q,mtpole,asmz,mzpole) ! SM, DRbar running masses at q
      mbsm = runb(q,mbmb,mtpole,asmz,mzpole)

*     determine the running couplings and quark/squark masses

      vv  = 1d0/sqrt(2d0*sq2*GF)    ! note: v ~ 174
      g = sqrt(2d0*mwpole**2/vv**2)
      gp = sqrt(2d0*(mzpole**2-mwpole**2)/vv**2)
      mz = mzpole
      mw = mwpole

      v1 = vv*cb
      v2 = vv*sb
      xx = mu/l

      A0 = l*xx*(Al+k*xx)/sb/cb      ! the would-be mA^2 in the MSSM limit

* mod. by UE:
	IF(A0.le.1.d4) A0=1.d4
* asma=Alpha_s at the scale A0:
      asma = slasf(dsqrt(A0),mtpole,asmz,mzpole)
* ast=Alpha_s at the scale mtsm:
      ast = slasf(mtsm,mtpole,asmz,mzpole)
* (Used in subroutine squarks)
* End mod. by UE      
      
      
      do i=1,5
 
      CALL squarks(mtsm,mbsm,mz,mw,mq3,mtr,mbr,At,Ab,l*xx,tanb,mg,
     .      asq,A0,vv,q**2,mtmssm,mbmssm,mstop2,msbot2,cst,sst,
     .      csb,ssb,errsqua,mtpole,asma,ast,v2)

       IF(errsqua) THEN
         err=3
         RETURN
       ENDIF

 
      CALL treemasses(g,gp,l,k,mtmssm/v2,mbmssm/v1,mtau/v1,v1,v2,xx,
     .      M1,M2,Ak,Al,At,Ab,Atau,mQ3,mtr,mbr,mQ,mur,mdr,mL3,
     .      mtaur,mL,mer,Q,errmass)

       IF(errmass.ne.0) THEN
         err=errmass
         RETURN
       ENDIF

 
      CALL getPiZZ(g,gp,mtmssm/v2,mbmssm/v1,mtau/v1,
     .      v1,v2,mzpole,Q,piZZ_MZ)

 
      CALL getPiWW(g,gp,mtmssm/v2,mbmssm/v1,mtau/v1,
     .      v1,v2,mwpole,Q,piWW_MW)

 
      CALL getPiWW(g,gp,mtmssm/v2,mbmssm/v1,mtau/v1,
     .      v1,v2,0d0,Q,piWW_0)

 
      CALL getdVB(g,gp,mzpole,mwpole,Q,dVB)

       vv = 1d0/sqrt(2d0*sq2*GF)/sqrt(1d0-piWW_0/mwpole**2-dVB)
       v1 = vv*cb
       v2 = vv*sb

       g = sqrt(2d0*mwpole**2/vv**2*(1d0+piWW_MW/mwpole**2))
       gp = sqrt(2d0*mzpole**2/vv**2*(1d0+piZZ_MZ/mzpole**2)-g**2)

       mz = sqrt(g**2+gp**2)/sq2*vv
       mw = g/sq2*vv

      enddo

*     tree-level mass matrices

      As = Al+k*xx            ! shortcuts
      gb2 = (g**2+gp**2)/2d0

      MStr(1,1) = gb2*v1**2+l*xx*v2/v1*As
      MStr(1,2) = (2d0*l**2-gb2)*v1*v2-l*xx*As
      MStr(1,3) = 2d0*l**2*v1*xx-l*v2*(As+k*xx)
      MStr(2,2) = gb2*v2**2+l*xx*v1/v2*As
      MStr(2,3) = 2d0*l**2*v2*xx-l*v1*(As+k*xx)
      MStr(3,3) = l*Al*v1*v2/xx+k*xx*(Ak+4d0*k*xx)
      MStr(2,1) = MStr(1,2)
      MStr(3,1) = MStr(1,3)
      MStr(3,2) = MStr(2,3)

      gold = mz**2        ! gauge-fixing mass

      MPtr(1,1) = l*xx*v2/v1*As+cb**2*gold
      MPtr(1,2) = l*xx*As-cb*sb*gold
      MPtr(1,3) = l*v2*(As-3d0*k*xx)
      MPtr(2,2) = l*xx*v1/v2*As+sb**2*gold
      MPtr(2,3) = l*v1*(As-3d0*k*xx)
      MPtr(3,3) = 4d0*l*k*v1*v2+l*Al*v1*v2/xx-3d0*k*Ak*xx
      MPtr(2,1) = MPtr(1,2)
      MPtr(3,1) = MPtr(1,3)
      MPtr(3,2) = MPtr(2,3)

*     compute the top corrections at zero momentum

      T1 = mstop2(1)
      T2 = mstop2(2)


      CALL effpot(2,mtmssm,mg,T1,T2,sst,cst,q**2,
     .     tanb,vv,l,xx,asq,DMS_2lt,DMP_2lt)

*     compute the bottom corrections at zero momentum

      B1 = msbot2(1)
      B2 = msbot2(2)


      CALL effpot(2,mbmssm,mg,B1,B2,ssb,csb,q**2,
     .     1d0/tanb,vv,l,xx,asq,DMS_2lb,DMP_2lb)


      CALL swap12(DMS_2lb)

      CALL swap12(DMP_2lb)

*     compute the two-loop top-bot Yukawa corrections

 
      CALL twlpyuk(mtmssm,mbmssm,A0,T1,T2,B1,B2,sst,cst,ssb,csb,q**2,
     .     l,xx,tanb,vv,DMS_2ly,DMP_2ly)

*     compute the full one-loop at zero momentum


      CALL treemasses(g,gp,l,k,mtmssm/v2,mbmssm/v1,mtau/v1,v1,
     .     v2,xx,M1,M2,Ak,Al,At,Ab,Atau,mQ3,mtr,mbr,mQ,mur,mdr,mL3,
     .     mtaur,mL,mer,Q,errmass)

       IF(errmass.ne.0) THEN
         err=errmass
         RETURN
       ENDIF


      CALL scalarcouplings(g,gp,l,k,mtmssm/v2,mbmssm/v1,
     .     mtau/v1,v1,v2,xx,Al,Ak,At,Ab,Atau)


      CALL gettadS(g,gp,mtmssm/v2,mbmssm/v1,
     .     mtau/v1,v1,v2,Q,tadS)


      CALL getPiSS(g,gp,mtmssm/v2,mbmssm/v1,
     .     mtau/v1,v1,v2,0d0,Q,piSS)


      CALL getPiPP(g,gp,mtmssm/v2,mbmssm/v1,
     .     mtau/v1,v1,v2,0d0,Q,piPP)

*     put all together (except the momentum-dependent corrections)

      vev(1) = v1
      vev(2) = v2
      vev(3) = xx

      do i = 1,3
       do j = 1,3

          MS(i,j) = MStr(i,j) ! start with tree level
          MP(i,j) = MPtr(i,j)

          ! add the full two loop

             MS(i,j) = MS(i,j)+DMS_2lt(i,j)
             MP(i,j) = MP(i,j)+DMP_2lt(i,j)
             MS(i,j) = MS(i,j)+DMS_2lb(i,j)
             MP(i,j) = MP(i,j)+DMP_2lb(i,j)
             MS(i,j) = MS(i,j)+DMS_2ly(i,j)
             MP(i,j) = MP(i,j)+DMP_2ly(i,j)

          ! add zero-mom 1-loop
          ! (If mom=true: Mfull matrices are recalculated)

            MSfull(i,j) = MS(i,j)-piSS(i,j)
            MPfull(i,j) = MP(i,j)-piPP(i,j)

             if(i.eq.j) then

              MSfull(i,j) = MSfull(i,j)+tadS(i)/sqrt(2d0)/vev(i)
              MPfull(i,j) = MPfull(i,j)+tadS(i)/sqrt(2d0)/vev(i)

             endif

       enddo
      enddo

*     simple diagonalization
*    (check of pos. masses^2, and for the mixing matrices)

      CALL jacobi(MSfull,3,ms2,OS)

      CALL jacobi(MPfull,3,ma2,OP)

      CALL mysort(ms2,OS)

      CALL mysort(ma2,OP)

      IF(ms2(1).le.0d0) THEN
       mss(1)=ms2(1)
       err=1
       RETURN
      ENDIF

      IF(ma2(1).le.0d0) THEN
       maa(1)=ma2(1)
       err=2
       RETURN
      ENDIF

      if(mom) then    ! iterative procedure for the external mom

       do kk = 1,3          ! one round for each eigenstate

        ps = sqrt(ms2(kk))

        if(abs(ma2(kk)).le.1d-6) then
         pa = 0d0      ! otherwise the PV functions freak out
        else
         pa = sqrt(ma2(kk))
        endif
      
        IL=0
 500    continue
        IL=IL+1

        CALL getPiSS(g,gp,mtmssm/v2,mbmssm/v1,
     .         mtau/v1,v1,v2,ps,Q,piSS)
    
        CALL getPiPP(g,gp,mtmssm/v2,mbmssm/v1,
     .         mtau/v1,v1,v2,pa,Q,piPP)

        do i=1,3
         do j = 1,3

          MSfull(i,j) = MS(i,j)-piSS(i,j)
          MPfull(i,j) = MP(i,j)-piPP(i,j)

          if(i.eq.j) then

           MSfull(i,j) = MSfull(i,j)+tadS(i)/sqrt(2d0)/vev(i)
           MPfull(i,j) = MPfull(i,j)+tadS(i)/sqrt(2d0)/vev(i)

          endif
         enddo
        enddo

        CALL jacobi(MSfull,3,ms2cor,OS0)
    
        CALL jacobi(MPfull,3,ma2cor,OP0)

        CALL mysort(ms2cor,OS0)
    
        CALL mysort(ma2cor,OP0)

        if(abs(ps**2-abs(ms2cor(kk)))/ps**2
     .         +abs(pa**2-abs(ma2cor(kk)))/pa**2.gt.1d-4) then

         if(IL.ge.10) then
          err=4
          return
         endif

         ps = sqrt(abs(ms2cor(kk)))
         pa = sqrt(abs(ma2cor(kk)))
         goto 500

        else

         ms2(kk) = ms2cor(kk)
         ma2(kk) = ma2cor(kk)

        endif

       enddo
      endif

*     take the square roots and check that it's all right

      do i = 1,3
       if(ms2(i).ge.0d0) then
          mss(i) = sqrt(ms2(i))
       else
          mss(i) = -sqrt(abs(ms2(i)))
          err=1
          return
       endif
      enddo

      do i = 1,3
       if(ma2(i).ge.0d0) then
          maa(i) = sqrt(ma2(i))
       else
          maa(i) = -sqrt(abs(ma2(i)))
          err=2
          return
       endif
      enddo

* Charged Higgs mass (by courtesy of P. Slavich and K.H. Phan)

      mhcsqtree = (l*xx*As-l**2*v1*v2)/sb/cb + mw**2

* UE: Estimate 2-loop corrs. borrowed from the CP-odd Higgs:

      mhcsqtree = mhcsqtree +
     .     (DMP_2lt(1,2)+DMP_2lb(1,2)+DMP_2ly(1,2))/sb/cb

      mhcsqcorr = mhcsqtree

      phc = sqrt(dabs(mhcsqcorr))

      call getPiHpHm(g,gp,l,k,mtmssm/v2,mbmssm/v1,
     .        mtau/v1,v1,v2,xx,Al,At,Ab,Atau,phc,Q,piHpHm)

      mhcsqcorr = mhcsqtree - piHpHm 
     .        + sb**2*tadS(1)/v1/sqrt(2d0) + cb**2*tadS(2)/v2/sqrt(2d0)
     
      phc = sqrt(mhcsqcorr)
     
      if(mom) then

       IL=0
 600   continue
       IL=IL+1

       call getPiHpHm(g,gp,l,k,mtmssm/v2,mbmssm/v1,
     .           mtau/v1,v1,v2,xx,Al,At,Ab,Atau,phc,Q,piHpHm)

       mhcsqcorr = mhcsqtree - piHpHm
     .        + sb**2*tadS(1)/v1/sqrt(2d0) + cb**2*tadS(2)/v2/sqrt(2d0)

        if(abs(phc**2-abs(mhcsqcorr))/phc**2.gt.1d-4)then

         if(IL.ge.10) then
          err=4
          return
         endif

         phc = sqrt(abs(mhcsqcorr))
         goto 600

        endif
       endif

*     now take the square root (if allowed)

       if(mhcsqcorr.gt.0d0) then
          mhc = sqrt(mhcsqcorr)
       else
          mhc=-dsqrt(dabs(mhc))
          err=5
       endif

      end

*
***********************************************************************
*

      subroutine mysort(msd,ZS)

      implicit none

      double precision msd(3),ZS(3,3),mss(3),OS(3,3)
      integer i,j,kk

*     order the disorder

      mss(1) = dmin1(msd(1),msd(2),msd(3))
      mss(3) = dmax1(msd(1),msd(2),msd(3))

      do i=1,3
       if(msd(i).gt.mss(1).and.msd(i).lt.mss(3)) then
          mss(2) = msd(i)
       endif
      enddo

      do i = 1,3
       do j = 1,3
          if(mss(i).eq.msd(j)) then
             do kk = 1,3
              OS(i,kk) = ZS(kk,j)
             enddo
          endif
       enddo
      enddo

      do i=1,3
       msd(i) = mss(i)
       do j=1,3
          ZS(i,j)=OS(i,j)
       enddo
      enddo

      end

*
***********************************************************************
*

      subroutine getdVB(g,gp,mzpole,mwpole,Q,dVB)

      implicit none

      double precision g,gp,mzpole,mwpole,Q,dVB
      double precision ch2,sh2,c2,s2,pi,rho

      pi = 4d0*atan(1d0)

      ch2 = g**2/(g**2+gp**2)
      sh2 = 1d0-ch2

      c2 = (mwpole/mzpole)**2
      s2 = 1d0-c2

      rho = c2/ch2

      dVB = 6d0+log(c2)/s2*(3.5d0-2.5d0*s2-sh2*(5d0-1.5d0*c2/ch2))
      dVB = dVB-4d0*log(mzpole**2/Q**2)
      dVB = g**2/16d0/pi**2*dVB
*     dVB = dVB*rho

      end

*
***********************************************************************
*
      subroutine effpot(lp,mt,mg,T1,T2,st,ct,q2,tanb,vv,l,xx,as,
     .     DMS,DMP)

      implicit none

      integer i,j,lp
      double precision mt,mg,T1,T2,st,ct,q2,tanb,vv,l,xx,as,
     .     DMS(3,3),DMP(3,3)
      double precision c2t,s2t,At,mu,Xt,ht,sbe,pi,k
      double precision F1t,F2t,F3t,Ft,FA

      pi = 4d0*atan(1d0)

      if(lp.eq.1) then
       k = 3d0/(16d0*Pi**2)   ! one-loop factor
      elseif(lp.eq.2) then
       k = as/(16d0*Pi**3)    ! two-loop factor
      else
       k = 0d0
      endif

      s2t = 2d0*ct*st
      c2t = ct**2-st**2

      mu = l*xx
      Xt = (T1-T2)*s2t/2d0/mt
      At = Xt+mu/tanb

      sbe = dsin(datan(tanb))

      ht = mt/vv/sbe          ! v ~ 174

      if(lp.eq.1) then        !the usual one-loop functions

       Ft = T1*(log(T1/q2)-1d0)-T2*(log(T2/q2)-1d0)

       F1t = log(T1*T2/mt**4)

       F2t = log(T1/T2)

       F3t = 2d0-(T1+T2)/(T1-T2)*log(T1/T2)

       FA = Ft

      elseif(lp.eq.2) then

 
      CALL makefuncs(mt,mg,T1,T2,s2t,c2t,q2,tanb,At,mu,
     .      F1t,F2t,F3t,Ft,FA)

      endif

*     now build up the results

      DMS(1,1) = .5d0*ht**2*mu**2*s2t**2*F3t
     .    +ht**2*tanb*mu*At/(T1-T2)*Ft

      DMS(1,2) =-ht**2*mu*mt*s2t*F2t-.5d0*ht**2*At*mu*s2t**2*F3t
     .    -ht**2*mu*At/(T1-T2)*Ft

      DMS(1,3) = .5d0*ht*l*mu*mt*s2t**2/tanb*F3t
     .    -ht*l*mt*(At-2d0*mu/tanb)/(T1-T2)*Ft

      DMS(2,2) = 2d0*ht**2*mt**2*F1t+2d0*ht**2*At*mt*s2t*F2t
     .    +.5d0*ht**2*At**2*s2t**2*F3t
     .    +ht**2/tanb*mu*At/(T1-T2)*Ft

      DMS(2,3) = -.5d0*ht*l*At*mt*s2t**2/tanb*F3t
     .    -ht*l*mt**2*s2t/tanb*F2t-ht*l*mt*At/(T1-T2)/tanb*Ft

      DMS(3,3) = .5d0*l**2*s2t**2*mt**2/tanb**2*F3t
     .     +l**2*mt**2/tanb*At/mu/(T1-T2)*Ft

      DMS(2,1) = DMS(1,2)
      DMS(3,1) = DMS(1,3)
      DMS(3,2) = DMS(2,3)

      DMP(1,1) = ht**2*mu*At/(T1-T2)*FA*tanb

      DMP(1,2) = ht**2*mu*At/(T1-T2)*FA

      DMP(1,3) = l*ht*mt*At/(T1-T2)*FA

      DMP(2,2) = ht**2*mu*At/(T1-T2)*FA/tanb

      DMP(2,3) = l*ht*mt*At/(T1-T2)*FA/tanb

      DMP(3,3) = l**2*mt**2*At/mu/(T1-T2)*FA/tanb

      DMP(2,1) = DMP(1,2)
      DMP(3,1) = DMP(1,3)
      DMP(3,2) = DMP(2,3)

      do i=1,3
       do j=1,3
          DMS(i,j) = k*DMS(i,j)
          DMP(i,j) = k*DMP(i,j)
       enddo
      enddo

      end

*
***********************************************************************
*

      subroutine makefuncs(mt,mg,T1,T2,s2t,c2t,q,tanb,At,mu,
     .     F1t,F2t,F3t,Ft,FA)

      implicit none

      double precision mt,mg,T1,T2,s2t,c2t,q,tanb,At,mu,
     .     F1t,F2t,F3t,Ft,FA

      double precision DT1,DT2,Dc2t,DT1T1,DT2T2,Dtt,Dc2tc2t,
     .     DT1t,DT2t,DT1T2,Dtc2t,DT1c2t,DT2c2t,Dcptmptt,
     .     Dtt_1,Dc2t_1,Dc2tc2t_1,Dtc2t_1,Dcptmptt_1,
     .     Dtt_2,Dc2t_2,Dc2tc2t_2,Dtc2t_2,Dcptmptt_2


      CALL makederiv(mt,mg,T1,T2,s2t,c2t,q,
     .     DT1,DT1T1,DT1t,DT1c2t,DT1T2,
     .     Dtt_1,Dc2t_1,Dc2tc2t_1,Dtc2t_1,Dcptmptt_1)


      CALL makederiv(mt,mg,T2,T1,-s2t,c2t,q,
     .     DT2,DT2T2,DT2t,DT2c2t,DT1T2,
     .     Dtt_2,Dc2t_2,Dc2tc2t_2,Dtc2t_2,Dcptmptt_2)

      Dtt = Dtt_1+Dtt_2
      Dc2t = Dc2t_1+Dc2t_2
      Dc2tc2t = Dc2tc2t_1+Dc2tc2t_2
      Dtc2t = Dtc2t_1+Dtc2t_2
      Dcptmptt = Dcptmptt_1+Dcptmptt_2

      F1t = Dtt+DT1T1+DT2T2+2d0*(DT1t+DT2t+DT1T2)

      F2t = DT1T1-DT2T2+DT1t-DT2t
     .     -4d0*c2t**2/(T1-T2)*(Dtc2t+DT1c2t+DT2c2t)

      F3t = DT1T1+DT2T2-2d0*DT1T2
     .    -2d0/(T1-T2)*(DT1-DT2)
     .    +16d0*c2t**2/(T1-T2)**2*(c2t**2*Dc2tc2t+2d0*Dc2t)
     .     -8d0*c2t**2/(T1-T2)*(DT1c2t-DT2c2t)

      Ft = DT1-DT2-4d0*c2t**2/(T1-T2)*Dc2t

      FA = Ft-2d0*mu/tanb/At/(T1-T2)/s2t**2*Dcptmptt

      end

*
***********************************************************************
*

      subroutine makederiv(mt,mg,T1,T2,s2t,c2t,q,
     .     DT1,DT1T1,DT1t,DT1c2t,DT1T2,Dtt,Dc2t,Dc2tc2t,Dtc2t,Dcptmptt)

      implicit none

      double precision mt,mg,T1,T2,s2t,c2t,q,
     .     DT1,DT1T1,DT1t,DT1c2t,DT1T2,Dtt,Dc2t,Dc2tc2t,Dtc2t,Dcptmptt
      double precision delt,phi,II,JJ

      double precision t,g,Logt,Logg,LogT1,LogT2,pphi,del,III

      t = mt**2
      g = mg**2

      Logt = Log(t/q)
      Logg = Log(g/q)
      LogT1 = Log(T1/q)
      LogT2 = Log(T2/q)
      pphi = phi(T1,g,t)
      del = delt(T1,g,t)
      III = II(q,T1,g,t)

      Dc2t = .5d0*JJ(q,T1,T1)-.5d0*JJ(q,T1,T2)+2d0*mg*mt/s2t*III

      Dc2tc2t = mg*mt/s2t**3*III

      Dcptmptt = -4d0*mg*mt*s2t*III

      DT1= -6d0*T1+2d0*mg*mt*s2t+4d0*t*(1d0-logt+logT1)
     .     +4d0*g*(1d0-logg+logT1)+((5d0-c2t**2)*T1-s2t**2*T2
     .     -4d0*mg*mt*s2t)*logT1+(-3d0+c2t**2)*T1*logT1**2
     .     +s2t**2*T2*logT1*logT2-(2d0*(g+t-T1)-2d0*mg*mt*s2t)
     .     *(logt*(logT1-logg)+logT1*logg)+(2d0/t*(del+2d0*g*t)
     .     -2d0*mg/mt*s2t*(g+t-T1))*pphi

      DT1T1= -(1d0+c2t**2)+4d0/T1*(g+t-mg*mt*s2t)-s2t**2*T2/T1
     .     *(1d0-logT2)+(3d0+c2t**2+8d0*g*t/del-4d0*mg*mt*s2t/del
     .     *(g+t-T1))*logT1-4d0*t/del/T1*(del-g*(g-t-T1)+mg*mt*s2t
     .     *(g-t+T1))*logt-4d0*g/del/T1*(del+t*(g-t+T1)-mg*mt*s2t
     .     *(g-t-T1))*logg+(-3d0+c2t**2)*logT1**2
     .     +2d0*(logt*(logT1-logg)+logT1*logg)-2d0/t/del*((g+t-T1)
     .     *(del-2d0*g*t)+4d0*mg**3*mt**3*s2t)*pphi

      DT1c2t= (T2*(1d0-logT2)-T1*(1d0-logT1))*LogT1
     .     -mg*mt/s2t*(1d0-2d0*logT1+logt*(logT1-logg)
     .     +logT1*logg-(g+t-T1)/t*pphi)

      DT1t= mg/mt*s2t+4d0*g/del*(T1-g-t+2d0*mg*mt*s2t)*logg
     .     +4d0/del*(2d0*g*t-mg*mt*s2t*(g+t-T1))*logt+2d0/del
     .     *(2d0*g*(g-t-T1)-mg/mt*s2t*(del-2d0*t*(t-g-T1)))*logT1
     .     +(-2d0+mg/mt*s2t)*(logt*(logT1-logg)+logT1*logg)
     .     +1d0/del/t*(mg/mt*s2t*(del*(T1-g-3d0*t)
     .     +2d0*t*((t-T1)**2-g**2))
     .     +2d0*(g-T1)**3+2d0*t*(del+(2d0*T1-t)*(g+T1)))*pphi

      DT1T2= s2t**2*logT1*logT2

      Dtt= -2d0-5d0/2d0*mg/mt**3*s2t*T1+6d0*logt**2
     .     +4d0*g/del*(g-t-T1+mt/mg*s2t*(t-g-T1))*logt
     .     -4d0*g/del*(g-t+T1+mg/mt*s2t*(t-g+T1))*logg
     .     +(8d0*g*T1/del+2d0*mg/mt**3*s2t*T1
     .      *(1d0-2d0*t/del*(t+g-T1)))*logT1
     .     -(2d0-mg*s2t*T1/2d0/mt**3)*logg*(logt-logT1)
     .     -(2d0+mg*s2t*T1/2d0/mt**3)*logt*logT1
     .     +mg*s2t/2d0/mt**3*(g+3d0*t)*logT1*(logt-logg)
     .     -2d0/del/t*(mg/mt**3*s2t*(del**2/4d0+t*(g-2d0*t+T1)*del
     .     +t**2*(g-t+T1)**2)-T1*(del+(g+t)*(2d0*t-T1))-(g-t)**3)*pphi

      Dtc2t = -mg/mt/s2t/2d0*(5d0*T1
     .     -4d0*T1*logT1+(g-3d0*t)*logT1*(logg-logt)
     .     +T1*(logt*(logT1-logg)+logg*logT1)+(del/t-2d0*(g-t+T1))*pphi)

      end

*
***********************************************************************
*

      double precision function myA0(m,q)
      double precision m,q

      if(m.ne.0d0) then
       myA0 = m*(1d0-Log(m/q))
      else
       myA0 = 0d0
      endif

      end

*
***********************************************************************
*

      double precision function myB0(q,m1,m2,mu2)

*     from Degrassi and Sirlin, Phys. Rev. D46 (1992) 3104.

      double precision q,m1,m2,Omega,mu2

      if(q.eq.0d0) then

       if(m1.eq.0d0.and.m2.ne.0d0) then
          myB0 = 1d0-Log(m2/mu2)
       elseif(m1.ne.0d0.and.m2.eq.0d0) then
          myB0 = 1d0-Log(m1/mu2)
       elseif(abs(m1-m2).le.1d-8) then
          myB0 = -Log(m1/mu2)
       else
          myB0 = 1d0-Log(m2/mu2)+m1/(m1-m2)*Log(m2/m1)
       endif

      else

       if(m1.eq.0d0.and.m2.ne.0d0) then

          if(m2.ne.q) then
             myB0 = -(Log(m2/mu2)-2d0-(m2/q-1d0)*Log(abs(1d0-q/m2)))
          else
             myB0 = -(Log(m2/mu2)-2d0)
          endif

       elseif(m2.eq.0d0.and.m1.ne.0d0) then

          if(m1.ne.q) then
             myB0 = -(Log(m1/mu2)-2d0-(m1/q-1d0)*Log(abs(1d0-q/m1)))
          else
             myB0 = -(Log(m1/mu2)-2d0)
          endif

       elseif(m2.eq.0d0.and.m1.eq.0d0) then

          myB0 = -(Log(q/mu2)-2d0) ! cut the imaginary part (I Pi)

       else

          myB0 = -( log(q/mu2)-2d0+
     .         1d0/2d0*( 1d0+(m1/q-m2/q))*log(m1/q) +
     .         1d0/2d0*( 1d0-(m1/q-m2/q))*log(m2/q) +
     .         2d0*Omega(m1/q,m2/q))

       endif

      endif

      end

*     function Omega(a,b) contained in myB0
       double precision function Omega(a,b)
       double precision a,b,cbig,sqCbig
       cbig = 0.5d0*(a+b)-0.25d0*(a-b)*(a-b)-0.25d0
       if(Cbig.gt.0d0) then
        sqCbig = sqrt(Cbig)
        Omega = sqCbig*
     .       (atan((1d0+a-b)/(2d0*sqCbig)) +
     .      atan((1d0-a+b)/(2d0*sqCbig)) )
       elseif(Cbig.lt.0d0) then
        sqCbig = sqrt(-Cbig)
	if((a+b-1d0+2d0*sqCbig).eq.0d0)then
	 Omega = 0.5d0*log(2d0/(a+b))
	else
         Omega = 0.5d0*sqCbig*
     .       log((a+b-1d0-2d0*sqCbig)/
     .       (a+b-1d0+2d0*sqCbig))
        endif
       else
        Omega = 0d0
       endif

       end

*
**********************************************************************
*

      double precision function myB1(p,m1,m2,q)

      implicit none

      double precision p,m1,m2,q
      double precision myA0,myB0

      if(p.eq.0d0) then
       if(abs(m1-m2).le.1d-8) then
          myB1 = -Log(m1/q)/2d0
       else
          if(m1.eq.0d0) then
             myB1 = (1d0-2d0*Log(m2/q))/4d0
          elseif(m2.eq.0d0) then
             myB1 = (3d0-2d0*Log(m1/q))/4d0
          else
             myB1 = (1d0-Log(m2/q)+m1**2/(m1-m2)**2*Log(m2/m1)
     .            +(m1+m2)/(m1-m2)/2d0)/2d0
          endif
       endif
      else
       myB1 = (myA0(m2,q)-myA0(m1,q)+(p+m1-m2)*myB0(p,m1,m2,q))/2d0/p
      endif

      end

*
**********************************************************************
*

      double precision function myF(q,m1,m2,mu2)

      implicit none
      double precision q,m1,m2,mu2,myA0,myB0

      myF = myA0(m1,mu2)-2d0*myA0(m2,mu2)
     .     -(2d0*q+2d0*m1-m2)*myB0(q,m1,m2,mu2)

      end

*
***********************************************************************
*

      double precision function myG(q,m1,m2,mu2)

      implicit none
      double precision q,m1,m2,mu2,myA0,myB0

      if(q.eq.0d0.and.m1.eq.0d0.and.m2.eq.0d0) then
       myG = 0d0
      else
       myG = (q-m1-m2)*myB0(q,m1,m2,mu2)-myA0(m1,mu2)-myA0(m2,mu2)
      endif

      end

*
***********************************************************************
*

      double precision function myB22(q,m1,m2,mu2)

      implicit none
      double precision q,m1,m2,mu2,myA0,myB0,myB1

      if(q.eq.0d0.and.m1.eq.0d0.and.m2.eq.0d0) then
       myB22 = 0d0
      else
       myB22 = ((myA0(m1,mu2)+myA0(m2,mu2))/2d0
     .      +(m1+m2-q/2d0)*myB0(q,m1,m2,mu2)
     .      +(m2-m1)*(myB1(q,m1,m2,mu2)-myB0(q,m1,m2,mu2)/2d0)
     .      +m1+m2-q/3d0)/6d0
      endif

      end

*
***********************************************************************
*

      double precision function myB22T(q,m1,m2,mu2)

      implicit none
      double precision q,m1,m2,mu2,myA0,myB22

      myB22T = myB22(q,m1,m2,mu2)-myA0(m1,mu2)/4d0-myA0(m2,mu2)/4d0

      end

*
**********************************************************************
*

      double precision function myH(q,m1,m2,mu2)

      implicit none
      double precision q,m1,m2,mu2,myG,myB22

      myH = 4d0*myB22(q,m1,m2,mu2)+myG(q,m1,m2,mu2)

      end

*
**********************************************************************
*

      double precision function JJ(q,m1,m2)

      implicit none
      double precision q,m1,m2

      JJ = m1*m2*(Log(m1/q)-1d0)*(Log(m2/q)-1d0)

      end

*
**********************************************************************
*

      double precision function II(q,m1,m2,m3)

      implicit none
      double precision q,m1,m2,m3,delt,phi

      II = (m1-m2-m3)/2d0*Log(m2/q)*Log(m3/q)
     .     +(m2-m1-m3)/2d0*Log(m1/q)*Log(m3/q)
     .     +(m3-m1-m2)/2d0*Log(m1/q)*Log(m2/q)
     .     +2d0*m1*log(m1/q)+2d0*m2*Log(m2/q)+2d0*m3*Log(m3/q)
     .     -2.5d0*(m1+m2+m3)-delt(m1,m2,m3)/2d0/m3*phi(m1,m2,m3)

      end

*
**********************************************************************
*

      function phi(x,y,z)

*     from Davydychev and Tausk, Nucl. Phys. B397 (1993) 23

      implicit none
      double precision x,y,z,phi,pphi,myphi

      if(x.le.z.and.y.le.z) then
       pphi = myphi(x,y,z)
      elseif(z.le.x.and.y.le.x) then
       pphi = z/x*myphi(z,y,x)
      elseif(z.le.y.and.x.le.y) then
       pphi = z/y*myphi(z,x,y)
      endif

      phi = pphi

      end

      function myphi(x,y,z)

      implicit none

      double precision x,y,z,myphi
      double precision u,v
      double precision Pi,SLLi2
      double complex clam,cxp,cxm,SLCLI2,ccphi

      parameter (pi = 3.1415926535897932384626433832795029d0)

*     auxiliary variables

      u = x/z
      v = y/z

      if(u.le.1d-8) then

       if(v.ne.1d0) then
          myphi = (log(u)*log(v)+2d0*SLLi2(1d0-v))/(1d0-v)
       else
          myphi = 2d0-log(u)
       endif

      elseif(v.le.1d-8) then

       if(u.ne.1d0) then
          myphi = (log(v)*log(u)+2d0*SLLi2(1d0-u))/(1d0-u)
       else
          myphi = 2d0-log(v)
       endif

      else

       if((1d0-u-v)**2.ge.4d0*u*v) then
          clam = DCMPLX(sqrt((1d0-u-v)**2-4d0*u*v),0d0)
       else
          clam = DCMPLX(0d0,sqrt(4d0*u*v-(1d0-u-v)**2))
       endif

       cxp = (1d0+(u-v)-clam)/2d0
       cxm = (1d0-(u-v)-clam)/2d0

*     phi function from eq. (A4)

       ccphi = (2d0*log(cxp)*log(cxm)-log(u)*log(v)-
     .      2d0*(SLCLI2(cxp)+SLCLI2(cxm))+Pi**2/3d0)/clam
       myphi = DBLE(ccphi)

      endif

      end

*
***********************************************************************
*

      function SLLi2(x)

      implicit none

      double complex SLCLI2,z
      double precision x,SLLi2

      z = DCMPLX(x,0d0)
      SLLi2 = DBLE(SLCLI2(z))

      end

*
***********************************************************************
*

      DOUBLE COMPLEX FUNCTION SLCLI2(Z)

*     just CALL the Dilog routine

      DOUBLE COMPLEX Z,Dilog

      SLCLI2 = Dilog(Z)

      end

*
**********************************************************************
*
* Dilog.F
* complex dilogarithm
* this file is part of FeynHiggs
* last modified 20 Oct 05 th

      double complex function Dilog(z)
      implicit none
      double complex z

      double complex Dilogsum
      external Dilogsum

      double precision absz, abs1z
      double complex t, mlogz

      double precision pi, zeta2
      parameter (pi = 3.1415926535897932384626433832795029d0)
      parameter (zeta2 = pi*pi/6d0)

      absz = abs(z)
      if( absz .lt. 1D-20 ) then
       Dilog = -log(1d0-z)
       return
      endif

      abs1z = abs(1d0-z)
      if( abs1z .lt. 1D-20 ) then
         Dilog = zeta2
         return
      endif

      if( DBLE(z) .gt. .5d0 ) then
         mlogz = -log(z)
         t = zeta2+mlogz*log(1d0-z)
         if( abs1z .gt. 1 ) then
            Dilog = Dilogsum(log(1d0-1d0/z))+zeta2 +
     .           .5d0*log(z-1d0)**2+t
         else
          Dilog = -Dilogsum(mlogz)+t
       endif
      else
       if( absz .gt. 1 ) then
          Dilog = -Dilogsum(-log(1d0-1d0/z))-zeta2-.5d0*log(-z)**2
       else
          Dilog = Dilogsum(-log(1d0-z))
       endif
      endif
      end

************************************************************************

      double complex function Dilogsum(w)
      implicit none
      double complex w

      double complex u, t
      integer k

      double precision b2, b4, b6, b8, b10, b12, b14
      double precision b16, b18, b20, b22
      parameter (b2 = 1d0/6d0)
      parameter (b4 = -1d0/30d0)
      parameter (b6 = 1d0/42d0)
      parameter (b8 = -1d0/30d0)
      parameter (b10 = 5d0/66d0)
      parameter (b12 = -691d0/2730d0)
      parameter (b14 = 7d0/6d0)
      parameter (b16 = -3617d0/510d0)
      parameter (b18 = 43867d0/798d0)
      parameter (b20 = -174611d0/330d0)
      parameter (b22 = 854513d0/138d0)

      double precision bernoulliB(11)
      data bernoulliB /b2, b4, b6, b8, b10, b12, b14,
     .     b16, b18, b20, b22/

      Dilogsum = w*(1d0-.25d0*w)
      if( abs(w) .lt. 1d-10 ) return

      u = w
      do k = 1, 11
       u = u*w**2/DBLE(2d0*k*(2d0*k+1d0))
       t = u*bernoulliB(k)
       Dilogsum = Dilogsum+t
       if( abs(t) .lt. 1d-16*abs(Dilogsum) ) return
      enddo

      end

      function delt(x,y,z)
      double precision delt,x,y,z

      delt = x**2+y**2+z**2-2d0*(x*y+x*z+y*z)

      end

***********************************************************************
*     FUNCTION IN THE HALL-RATTAZZI-SARID TERM
***********************************************************************

      function SLH2(x,y)

      implicit none

      double precision x,y,SLH2,eps

      eps = 1d-8

      if(abs(x-y).ge.eps) then
       if(abs(x-1d0).ge.eps.and.abs(y-1d0).ge.eps) then
          SLH2 = x*log(x)/(1d0-x)/(x-y)+y*log(y)/(1d0-y)/(y-x)
       elseif(abs(x-1d0).ge.eps.and.abs(y-1d0).lt.eps) then
          SLH2 = (-1d0+x-x*log(x))/(x-1d0)**2
       elseif(abs(x-1d0).lt.eps.and.abs(y-1d0).ge.eps) then
          SLH2 = (-1d0+y-y*log(y))/(y-1d0)**2
       else
          SLH2 = -.5d0
       endif
      else
       if(abs(x-1d0).ge.eps) then
          SLH2 = (1d0-x+log(x))/(x-1d0)**2
       else
          SLH2 = -.5d0
       endif
      endif

      end

*
***********************************************************************
*

      subroutine diagonalize(mq,mw,mz,msql,msqr,Aq,mu,tb,iq,
     .     msq2,sth,cth)

      implicit none

      integer iq
      double precision mq,mw,mz,msql,msqr,Aq,mu,tb,msq2(2),sth,cth

      double precision mq2,mz2,mw2,c2b,xq,yq,zq,tth,dx,dy

      mq2 = mq**2
      mz2 = mz**2
      mw2 = mw**2
      c2b = (1d0-tb**2)/(1d0+tb**2)

      if(iq.eq.1) then
       zq = mq*(Aq-mu/tb)
       dx = 1d0/4d0*mz2*c2b
       dy = 1d0/12d0*(8d0*mw2-5d0*mz2)*c2b
      elseif(iq.eq.2) then
       zq = mq*(Aq-mu*tb)
       dx = -1d0/4d0*mz2*c2b
       dy = -1d0/12d0*(4d0*mw2-mz2)*c2b
      else
       write(*,*) 'ERROR: iq out of range'
      endif

      xq = mq2+1d0/2d0*(msql**2+msqr**2)+dx
      yq = 1d0/2d0*(msql**2-msqr**2)+dy

      msq2(1) = xq+sqrt(yq**2+zq**2)
      msq2(2) = xq-sqrt(yq**2+zq**2)

      if(zq.eq.0.and.yq.ge.0)then
       cth=1d0
       sth=0d0
      elseif(zq.eq.0.and.yq.lt.0)then
       cth=0d0
       sth=1d0
      else
       tth = 1d0/zq*(sqrt(yq**2+zq**2)-yq)
       cth = 1d0/sqrt(1d0+tth**2)
       sth = cth*tth
      endif

      end

*
***********************************************************************
*

      subroutine squarks(mtsm,mbsm,mz,mw,mqs,mts,mbs,At,Ab,mu,tb,mg,as,
     .     A0,vv,q2,mtmssm,mbmssm,mstop2,msbot2,cst,sst,csb,ssb,
     .     errsqua,mtpole,asma,ast,v2)

*     compute the quark running masses and the squark masses and mixing

      implicit none

      double precision mtsm,mbsm,mz,mw,mqs,mts,mbs,At,Ab,mu,tb,mg,as,
     .     A0,vv,q2,mtmssm,mbmssm,mstop2(2),msbot2(2),cst,sst,csb,ssb,
     .     mtpole,asma,ast,v2,htmt,pi,dla,sb2,fact,mtma,htma,dlqa

      logical errsqua

      double precision dmt,eb,dd,K

      integer i

      pi=4d0*atan(1d0)

      errsqua = .false.

*     first computation of the squarks using the SM running quark masses

      mtmssm = mtsm
      mbmssm = mbsm


      CALL diagonalize(mtmssm,mw,mz,mqs,mts,At,mu,tb,1,
     .     mstop2,sst,cst)


      CALL diagonalize(mbmssm,mw,mz,mqs,mbs,Ab,mu,tb,2,
     .     msbot2,ssb,csb)
      
      do i =1,5

*     compute the threshold corrections to the quark masses

 
      CALL getdmt(mtmssm,mstop2,sst,cst,mbmssm,msbot2,ssb,csb,
     .      vv,mg,A0,mu,tb,as,q2,dmt)
       
 
      CALL getdmb(mtmssm,mstop2,sst,cst,mbmssm,msbot2,ssb,csb,
     .      vv,mg,A0,mu,tb,as,q2,dd,eb)

*     recompute the quark masses

       mtmssm = mtsm+dmt
      
* mod. by UE: resum large logs ~ht^2*ln(Q^2/mt^2) (not included in dmt):

       htmt=mtpole/v2/(1d0+4d0*ast/(3d0*pi)+10.9d0*(ast/pi)**2)
       dla=(asma/ast)**(1d0/7d0)
       sb2=tb**2/(1d0+tb**2)
       fact=1d0-9d0*sb2*htmt**2/(8d0*pi*ast)*(1d0-dla)
      
       mtma=mtmssm*fact**(-1d0/6d0)
       htma=htmt*dla**4*fact**(-0.5d0)
       dlqa=(as/asma)**(1d0/7d0)
       mtmssm=mtma*(1d0-9d0*htma**2/(8d0*pi*asma)*
     .    (1d0-dlqa))**(-1d0/6d0)

* end mod. by UE

       K = (1d0-dd)/(1d0+eb*tb) ! "resummation"

       mbmssm = K*mbsm

*     now recompute the squark masses and mixing

 
      CALL diagonalize(mtmssm,mw,mz,mqs,mts,At,mu,tb,1,
     .      mstop2,sst,cst)

 
      CALL diagonalize(mbmssm,mw,mz,mqs,mbs,Ab,mu,tb,2,
     .      msbot2,ssb,csb)

      if(.not.mstop2(2).ge.0d0) then
       errsqua = .true.
       return
      endif

      if(.not.msbot2(2).ge.0d0) then
       errsqua = .true.
       return
      endif

      enddo

      end

*
***********************************************************************
*

      subroutine getdmt(mt,mstop2,sst,cst,mb,msbot2,ssb,csb,
     .      vv,mg,A0,mu,tb,as,q2,dmt)

      implicit none

      double precision mt,mstop2(2),sst,cst,mb,msbot2(2),ssb,csb,
     .     vv,mg,A0,mu,tb,as,q2,dmt

      double precision pi,cbe,sbe,ht,hb,mt2,mb2,mu2,mg2,T1,T2,B1,B2,
     .     dmts,dmty,myB0,myB1,cbe2,sbe2,ht2,hb2

      pi=4d0*atan(1d0)

      cbe = 1d0/sqrt(1d0+tb**2)
      sbe = tb*cbe

      cbe2 = cbe**2
      sbe2 = sbe**2

      ht = mt/vv/sbe
      hb = mb/vv/cbe

      ht2 = ht**2
      hb2 = hb**2

      T1 = mstop2(1)
      T2 = mstop2(2)

      B1 = msbot2(1)
      B2 = msbot2(2)

      mt2 = mt**2
      mb2 = mb**2
      mu2 = mu**2
      mg2 = mg**2

      dmts = as/3d0/pi*mt*(
     .     myB1(mt2,mg2,T1,q2)+myB1(mt2,mg2,T2,q2)
     .     -2d0*sst*cst*mg/mt*(myB0(mt2,mg2,T1,q2)-myB0(mt2,mg2,T2,q2)))

      
      IF(A0.ge.0d0) THEN
       dmty = mt/32d0/pi**2*(
* mod. by UE:
* Scale q of myB1 replaced; large logs are resummed in subroutine
* squarks
     .     ht2*(2d0*myB1(mt2,mt2,0d0,mt2)+myB1(mt2,mb2,0d0,mt2)
     .       +2d0*cbe2*(myB1(mt2,mt2,A0,A0)-myB1(mt2,mt2,0d0,mt2))
     .       +cbe2*(myB1(mt2,mb2,A0,A0)-myB1(mt2,mb2,0d0,mt2))
* end mod. by UE
     .       +myB1(mt2,mu2,T1,q2)+myB1(mt2,mu2,T2,q2))
     .     +hb2*(cbe2*myB1(mt2,mb2,0d0,q2)+sbe2*myB1(mt2,mb2,A0,q2))
     .     -2d0*ht*hb*sbe*cbe*mb/mt*
     .     (myB0(mt2,mb2,0d0,q2)-myB0(mt2,mb2,A0,q2))
     .     +(ht2*csb**2+hb2*ssb**2)*myB1(mt2,mu2,B1,q2)
     .     +(ht2*ssb**2+hb2*csb**2)*myB1(mt2,mu2,B2,q2)
     .     +2d0*ht*hb*ssb*csb*mu/mt*
     .     (myB0(mt2,mu2,B1,q2)-myB0(mt2,mu2,B2,q2)))
      ELSE
       dmty=0d0
      endif
      dmt = dmts+dmty

      end

*
***********************************************************************
*

      subroutine getdmb(mt,mstop2,sst,cst,mb,msbot2,ssb,csb,
     .      vv,mg,A0,mu,tb,as,q2,dd,eb)

      implicit none

      double precision mt,mstop2(2),sst,cst,mb,msbot2(2),ssb,csb,
     .     vv,mg,A0,mu,tb,as,q2,dd,eb

      double precision pi,cbe,sbe,ht,hb,mt2,mu2,mg2,T1,T2,B1,B2,
     .     myB0,myB1,ht2,hb2,At,Ab,dds,ebs,ddy,eby,
     .     lhiggs,hhiggs,higgsino

      pi=4d0*atan(1d0)

      cbe = 1d0/sqrt(1d0+tb**2)
      sbe = tb*cbe

      ht = mt/vv/sbe
      hb = mb/vv/cbe

      ht2 = ht**2
      hb2 = hb**2

      T1 = mstop2(1)
      T2 = mstop2(2)

      B1 = msbot2(1)
      B2 = msbot2(2)

      mt2 = mt**2
      mu2 = mu**2
      mg2 = mg**2

      At = sst*cst*(T1-T2)/mt+mu/tb
      Ab = ssb*csb*(B1-B2)/mb+mu*tb

*     the strong part

      dds = as/3d0/pi*(
     .    2d0*Ab*mg/(B1-B2)*(myB0(0d0,mg2,B1,q2)-myB0(0d0,mg2,B2,q2))
     .    -myB1(0d0,mg2,B1,q2)-myB1(0d0,mg2,B2,q2))

      ebs = -as/3d0/pi*(
     .    2d0*mu*mg/(B1-B2)*(myB0(0d0,mg2,B1,q2)-myB0(0d0,mg2,B2,q2)))

*     the Yukawa part (in the large-tanB limit)

      lhiggs = -ht2/4d0*(5d0-6d0*Log(mt2/q2))

      IF(A0.gt.0d0) THEN
       hhiggs = 2d0*ht2*(mt2-A0+A0*Log(A0/q2)-mt2*Log(mt2/q2))/(mt2-A0)
     .    +hb2*mt2/2d0/(mt2-A0)**2*(mt2-A0+mt2*Log(A0/mt2))
     .    +3d0*hb**2/4d0*(1d0-2d0*Log(A0/q2))
      ELSE
       hhiggs=0d0
      ENDIF

      higgsino = hb**2*(myB1(0d0,mu2,B1,q2)+myB1(0d0,mu2,B2,q2))
     .    +(hb2*cst**2+ht2*sst**2)*myB1(0d0,mu2,T1,q2)
     .    +(hb2*sst**2+ht2*cst**2)*myB1(0d0,mu2,T2,q2)
     .    -2d0*ht2*mu2/(T1-T2)*(myB0(0d0,mu2,T1,q2)-myB0(0d0,mu2,T2,q2))

      ddy = -(lhiggs+hhiggs+higgsino)/32d0/pi**2

      eby = ht**2/16d0/pi**2*At*mu/(T1-T2)*
     .    (T1/(T1-mu2)*Log(T1/mu2)-T2/(T2-mu2)*Log(T2/mu2))

*     all together

      dd = dds+ddy
      eb = ebs+eby

      end

*
***********************************************************************
*

      double precision function runt(x,tc,asc,zm)

*     compute the running (SM, MSbar) top mass

      implicit double precision (a-h,o-z)

      pi=4d0*atan(1d0)

      if(x.le.tc)then
       fn=5d0
      else
       fn=6d0
      endif

      b0=11d0-2d0*fn/3d0
      b1=102d0-38d0*fn/3d0
      g0=8d0
      g1=404d0/3d0-40d0*fn/9d0

      asx=slasf(x,tc,asc,zm)
      ast=slasf(tc,tc,asc,zm)
      rrr=tc*(asx/ast)**(g0/(2d0*b0))

*     this is the relation between mpole/mtrun(mtrun)
*      pol1=1d0+4d0*ast/(3d0*pi)+8.243d0*(ast/pi)**2

*     this is the relation between mpole/mt(mpole)
      pol2= 1d0+4d0*ast/(3d0*pi) +10.9d0*(ast/pi)**2

      corr=1d0+ast*g0/(4d0*pi*2d0*b0)*(-b1/b0+g1/g0)*(asx/ast-1d0)

      runt=rrr*corr/pol2

      runt = runt*(1d0-asx/pi/3d0-(asx/pi)**2*29d0/72d0) ! shift to DRbar

      end

*
***********************************************************************
*

      double precision function slasf(x,tc,asc,zm)

*     compute the running (SM, MSbar) alpha_s

      implicit double precision (a-h,o-z)

      pi = 4d0*atan(1d0)

      fn=5d0

      b0=11d0-2d0*fn/3d0
      b1=102d0-38d0*fn/3d0
      vvv=1d0-b0*asc/(2d0*pi)*log(zm/x)

      if(x.le.tc) then        ! five flavors

       slasf=asc/vvv*(1d0-b1/b0*asc/(4d0*pi*vvv)*log(vvv))

      else

       vvv=1d0-b0*asc/(2d0*pi)*log(zm/tc) ! first evolve up to q=mt

       ast=asc/vvv*(1d0-b1/b0*asc/(4d0*pi*vvv)*log(vvv))

       b0t=b0-2d0/3d0       ! six flavours
       b1t=b1-38d0/3d0
       vvv=1d0-b0t*ast/(2d0*pi)*log(tc/x) !     now evolve up to the scale >mt

       slasf=ast/vvv*(1d0-b1t/b0t*ast/(4d0*pi*vvv)*log(vvv))

      endif

      end

*
************************************************************************
*

      double precision function runb(x,mbmb,mt,asmz,mz)

*     compute the running (SM, DRbar) bottom mass

      implicit none

      double precision x,mbmb,mt,asmz,mz
      double precision pi,b0,b1,g0,g1,asx,asb,slasf,ast,rrr,corr,mbmt
      integer nf

      pi=4d0*atan(1d0)

      if(x.le.mt) then

       nf = 5             ! evolve from mb to x with nf=5
       b0 = 11d0-2d0*nf/3d0
       b1 = 102d0-38d0*nf/3d0

       g0 = 8d0
       g1 = 404d0/3d0-40d0*nf/9d0

       asx=slasf(x,mt,asmz,mz)
       asb=slasf(mbmb,mt,asmz,mz)

       rrr=mbmb*(asx/asb)**(g0/(2d0*b0))
       corr=1d0+asb*g0/(4d0*pi*2d0*b0)*(-b1/b0+g1/g0)*(asx/asb-1d0)

       runb=rrr*corr

      else

       nf = 5             ! first evolve from mb to mt with nf=5
       b0 = 11d0-2d0*nf/3d0
       b1 = 102d0-38d0*nf/3d0

       g0 = 8d0
       g1 = 404d0/3d0-40d0*nf/9d0

       ast=slasf(mt,mt,asmz,mz)
       asb=slasf(mbmb,mt,asmz,mz)

       rrr=mbmb*(ast/asb)**(g0/(2d0*b0))
       corr=1d0+asb*g0/(4d0*pi*2d0*b0)*(-b1/b0+g1/g0)*(ast/asb-1d0)

       mbmt=rrr*corr

       nf = 6             ! then evolve from mt to x with nf=6
       b0 = 11d0-2d0*nf/3d0
       b1 = 102d0-38d0*nf/3d0

       g0 = 8d0
       g1 = 404d0/3d0-40d0*nf/9d0

       asx=slasf(x,mt,asmz,mz)

       rrr=mbmt*(asx/ast)**(g0/(2d0*b0))
       corr=1d0+ast*g0/(4d0*pi*2d0*b0)*(-b1/b0+g1/g0)*(asx/ast-1d0)

       runb=rrr*corr

      endif

      runb = runb*(1d0-asx/pi/3d0-(asx/pi)**2*29d0/72d0) ! shift to DRbar

      end

*
***********************************************************************
*

      subroutine swap12(M)

      implicit none
      double precision M(3,3),temp

      temp = M(1,1)
      M(1,1) = M(2,2)
      M(2,2) = temp

      temp = M(1,3)
      M(1,3) = M(2,3)
      M(2,3) = temp

      temp = M(3,1)
      M(3,1) = M(3,2)
      M(3,2) = temp

      end

*
***********************************************************************
*

      SUBROUTINE jacobi(a,np,d,v)

*     just CALLs Tomas Hahn's diagonalization routine

      INTEGER np,i,j
      double precision a(np,np),d(np),v(np,np)
      double complex M(np,np),U(np,np)

      do i=1,np              ! turn to complex
       do j=1,np
          M(i,j) = DCMPLX(a(i,j))
       enddo
      enddo


      CALL HEigensystem(np, M, np, d, U, np, 0)

      do i=1,np              ! back to real
       do j=1,np
          v(i,j) = DBLE(U(j,i)) ! the other jacobi had this convention
       enddo
      enddo

      end

* diagonalization of a Hermitian n-by-n matrix using the Jacobi algorithm
* code adapted from the "Handbook" routines for complex A
* (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
* this file is part of the Diag library
* last modified 27 Sep 07 th
************************************************************************
** HEigensystem diagonalizes a Hermitian n-by-n matrix.
** Input: n, A = n-by-n matrix, Hermitian
** (only the upper triangle of A needs to be filled).
** Output: d = vector of eigenvalues, U = transformation matrix
** these fulfill diag(d) = U A U^+ = U A U^-1 with U unitary.

      subroutine HEigensystem(n, A,ldA, d, U,ldU, sort)
      implicit none
      integer n, ldA, ldU, sort
      double complex A(ldA,*), U(ldU,*)
      double precision d(*)

      integer p, q, j
      double precision red, off, thresh
      double precision delta, t, invc, s
      double complex x, y, Apq
      double precision ev(2,16)

      integer sweep
      common /nsweeps/ sweep

      double precision sq
      double complex c
      sq(c) = DBLE(c*DCONJG(c))

      if( n .gt. 16 ) then
        print *, "Dimension too large"
        d(1) = -999d0
        return
      endif

      do p = 1, n
        ev(1,p) = 0d0
        ev(2,p) = DBLE(A(p,p))
        d(p) = ev(2,p)
      enddo

      do p = 1, n
        do q = 1, n
          U(q,p) = 0
        enddo
        U(p,p) = 1
      enddo

      red = .04d0/n**4

      do sweep = 1, 50
        off = 0
        do q = 2, n
          do p = 1, q-1
            off = off+sq(A(p,q))
          enddo
        enddo
        if( off .lt. 2d0**(-103) ) goto 1

        thresh = 0
        if( sweep .lt. 4 ) thresh = off*red

        do q = 2, n
          do p = 1, q-1
            off = sq(A(p,q))
            if( sweep .gt. 4 .and. off .lt.
     .            2d0**(-103)*max(ev(2,p)**2, ev(2,q)**2) ) then
            A(p,q) = 0
            else
            if( off .gt. thresh ) then
              t = .5d0*(ev(2,p)-ev(2,q))
              t = 1d0/(t+sign(sqrt(t**2+off), t))

              delta = t*off
              ev(1,p) = ev(1,p)+delta
              ev(2,p) = d(p)+ev(1,p)
              ev(1,q) = ev(1,q)-delta
              ev(2,q) = d(q)+ev(1,q)

              invc = sqrt(delta*t+1d0)
              s = t/invc
              t = delta/(invc+1d0)

              Apq = A(p,q)

              do j = 1, p-1
                x = A(j,p)
                y = A(j,q)
                A(j,p) = x+s*(DCONJG(Apq)*y-t*x)
                A(j,q) = y-s*(Apq*x+t*y)
              enddo

              do j = p+1, q-1
                x = A(p,j)
                y = A(j,q)
                A(p,j) = x+s*(Apq*DCONJG(y)-t*x)
                A(j,q) = y-s*(Apq*DCONJG(x)+t*y)
              enddo

              do j = q+1, n
                x = A(p,j)
                y = A(q,j)
                A(p,j) = x+s*(Apq*y-t*x)
                A(q,j) = y-s*(DCONJG(Apq)*x+t*y)
              enddo

              A(p,q) = 0

              do j = 1, n
                x = U(p,j)
                y = U(q,j)
                U(p,j) = x+s*(Apq*y-t*x)
                U(q,j) = y-s*(DCONJG(Apq)*x+t*y)
              enddo
            endif
            endif
          enddo
        enddo

        do p = 1, n
          ev(1,p) = 0
          d(p) = ev(2,p)
        enddo
      enddo

      print *, "Bad convergence in HEigensystem"

1     if( sort .eq. 0 ) return

* sort the eigenvalues

      do p = 1, n-1
        j = p
        t = d(p)
        do q = p+1, n
          if( sort*(t-d(q)) .gt. 0 ) then
            j = q
            t = d(q)
          endif
        enddo

        if( j .ne. p ) then
          d(j) = d(p)
          d(p) = t
          do q = 1, n
            x = U(p,q)
            U(p,q) = U(j,q)
            U(j,q) = x
          enddo
        endif
      enddo
      end
