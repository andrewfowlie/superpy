c ==================================================================== c
c                           top 2-body decays                          c
c ==================================================================== c
c Incl. dominant rad. corrs.:
c top --> bottom + W+/-: 
c   1-loop QCD as in Li, Oakes, Yuan, PRD43 (1991) 3759 for mb -> 0
c   no Susy (<~ 4%)
c top --> bottom + H+/-:
c   1-loop QCD as in Czarnecki, Davidson, hep-ph/9301237, for mb -> 0
c      (BUT: with terms ~ mb*tanb/mt)
c   Susy: ~delta_mb/mb only, from Guasch et al., hep-ph/9507461
c      and hep-ph/0003109
c ==================================================================== c

        SUBROUTINE TDECAY(PAR)

        IMPLICIT NONE 
        INTEGER I,J

        DOUBLE PRECISION PAR(*)
        DOUBLE PRECISION topneutrstop(5,2),brtopneutrstop(5,2)
        DOUBLE PRECISION atopr(2,5),btopr(2,5)
        DOUBLE PRECISION topbw,topbh,toptot,brtopbw,brtopbh
        DOUBLE PRECISION gmst(2),PI,SQR2,lamb_funct,rmt,rmb,integ
        DOUBLE PRECISION tanb,mhc,alsmt,runmb,sp,xmb,xmw,xmh
        DOUBLE PRECISION dtbwqcd,gp,gm,nb,nllr,delb
        DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     C     MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     C     CST,CSB,CSL,MSMU1,MSMU2,MSMUNT
        DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
        DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
        DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
        DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),N(5,5)
        DOUBLE PRECISION sst,sw,cw,tw,scalt

        COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
        COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     C     MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     C     CST,CSB,CSL,MSMU1,MSMU2,MSMUNT
        COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        COMMON/BR_top2body/brtopbw,brtopbh,brtopneutrstop
        COMMON/topwidth/toptot
        COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
        COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
        COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,N

        PI=4d0*DATAN(1d0)
        sqr2=dsqrt(2d0)
        tanb=PAR(3)
        mhc=cmass

* Parameters at M_top:

        ALSMT=ALSMZ/(1d0+23d0/(12d0*PI)*ALSMZ*DLOG(MT**2/MZ**2))
        rmt=MT/(1d0+4d0*ALSMT/(3d0*PI)+11d0*(ALSMT/PI)**2)
        rmb=RUNMB(MT)

        xmb=(rmb/rmt)**2
        xmw=(mw/rmt)**2
        xmh=(mhc/mt)**2

c -------------------------------------------------------------------- c
c top --> bottom + W+/-

      if(mt.gt.(mb+mw)) then
c   QCD correction (mb -> 0):
        dtbwqcd=2d0*alsmt/(3d0*pi)*(
     .    2d0*xmw*(2d0*xmw-1d0)*(1d0+xmw)*dlog(xmw)/
     .           ((1d0-xmw)**2*(1d0+2d0*xmw))
     .    -(5d0+4d0*xmw)/(1d0+2d0*xmw)*dlog(1-xmw)
     .    +2d0*(sp(1d0-xmw)-sp(xmw))-pi**2
     .    +(5d0+9d0*xmw+6d0*xmw**2)/(2d0*(1d0-xmw)*(1d0+2d0*xmw)))
         topbw = GF/(8d0*pi*sqr2)*mt**3*lamb_funct(mb/mt,mw/mt)*
     .     ((1-xmb)**2+(1+xmb)*xmw-2d0*xmw**2)*(1d0+dtbwqcd)
      else
         topbw = 0.D0
      endif
c -------------------------------------------------------------------- c
c top --> bottom + H+/-

      if(mt.gt.(mb+mhc)) then
c  QCD corrections (mb -> 0):
        gp=(1-xmh)*(sp(1-xmh)+9d0/8d0-pi**2/3d0
     .    +dlog(xmh)*dlog(1-xmh)/2d0+xmh*dlog(xmh)/(xmh-1)/2d0
     .    +(1/xmh/2d0-5d0/4d0)*dlog(1-xmh)+3d0/8d0*dlog(xmb))
        gm=-3d0/4d0*(1-xmh)*dlog(xmb)
c For tree level + Susy corrections:
        nb=(1+xmb-xmh)*(1/tanb**2+xmb*tanb**2)+4d0*xmb
        nllr=(1d0+xmb-xmh)*xmb*tanb**2+2d0*xmb
c Susy corrections:
        delb=2d0*alsmt/(3d0*pi)*MGL*(par(13)-par(4)*tanb)*
     .      integ(msb1,msb2,mgl)
     .     -rmt**2*GF*(1+tanb**2)*par(4)/(4d0*pi**2*sqr2*tanb)
     .      *(par(12)-par(4)/tanb)*integ(mst1,mst2,par(4))
c
         topbh = GF/(8d0*pi*sqr2)*mt**3*lamb_funct(mb/mt,mhc/mt)*(nb
     .     +2d0*nllr*delb
     .     +4d0*alsmt/(3d0*pi)*(2d0*(1/tanb**2+xmb*tanb**2)*gp
     .     +(1/tanb**2-xmb*tanb**2)*gm))
      else
         topbh = 0.D0
      endif

c -------------------------------------------------------------------- c
c top --> stop(j) + chi^0_i from Sdecay:

      do i=1,5,1
         do j=1,2,1
               topneutrstop(i,j) = 0.D0
         end do
      end do

       IF(dabs(mneu(1))+MST1.lt.mt) THEN

         sst=DSQRT(1.D0-CST**2)
         cw=DSQRT(G2/(G1+G2))
         sw=DSQRT(G1/(G1s+G2))
         tw=sw/cw
         scalt=rmt*dsqrt(2d0*sqr2*GF*(1+tanb**2))/(tanb*dsqrt(g2))

        gmst(1) = Mst1
        gmst(2) = Mst2

      do i=1,5,1
         atopr(1,i)=cst*sqr2*(-2d0*(N(i,1)*cw+N(i,2)*sw)*sw/3.d0
     .        +(-0.5d0+2d0/3d0*sw**2)*(-N(i,1)*sw+N(i,2)*cw)/cw)
     .        -sst*scalt*N(i,4)
         atopr(2,i)=-sst*sqr2*(-2d0*(N(i,1)*cw+N(i,2)*sw)*sw/3d0+(-.5d0
     .        +2d0/3d0*sw**2)*(-N(i,1)*sw+N(i,2)*cw)/cw)
     .        -cst*scalt*N(i,4)
         btopr(1,i)=-2.d0*sst*sqr2*sw*((-N(i,1)*sw+N(i,2)*cw)*tw
     .   -(N(i,1)*cw+N(i,2)*sw))/3.d0
     .        -cst*scalt*N(i,4)
         btopr(2,i)=-2.d0*cst*sqr2*sw*((N(i,1)*cw+N(i,2)*sw)*tw
     .   -(N(i,1)*cw+N(i,2)*sw))/3.d0
     .         +sst*scalt*N(i,4)
         do j=1,2,1
            if(mt.gt.(dabs(mneu(i))+gmst(j))) then
               topneutrstop(i,j) = g2/32.D0/pi/mt*(
     .              atopr(j,i)*btopr(j,i)*4.D0*mt*mneu(i) +
     .              (atopr(j,i)**2+btopr(j,i)**2)*
     .              (mt**2-gmst(j)**2+mneu(i)**2) )*
     .              lamb_funct(dabs(mneu(i))/mt,gmst(j)/mt)
            else
               topneutrstop(i,j) = 0.D0
            endif
         end do
      end do
      
      ENDIF
c -------------------------------------------------------------------- c

         toptot = topbw+topbh
         do i=1,5,1
            do j=1,2,1
               toptot = toptot + topneutrstop(i,j)
            end do
         end do
         
         brtopbw = topbw/toptot
         brtopbh = topbh/toptot

         do i=1,5,1
            do j=1,2,1
               brtopneutrstop(i,j) = topneutrstop(i,j)/toptot
            end do
         end do

      end
c -------------------------------------------------------------------- c

      double precision function lamb_funct(x,y)      

      implicit double precision (a-h,k-z)

      lamb_funct=dsqrt((1.D0-x**2-y**2)**2-4.D0*x**2*y**2) 

      return
      end
