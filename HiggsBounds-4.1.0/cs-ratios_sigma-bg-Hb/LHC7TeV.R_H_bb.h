!******************************************************
        real*8 function lhc7_rH_bb(x)
!******************************************************
!* ggH process : consensus numbers used (PDF: MSTW 2008 NNLO)
!*		 see: arXiv:1101.0593 [hep-ph] 
!* 		 numbers taken from: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt7TeV
!* bbH process : numbers in range [70:1020] generated with bbh@nnlo v. 1.3 (PDF: MSTW 2008 NNLO) 
!* 		 with input parameters according to arXiv:1101.0593 [hep-ph]
!******************************************************
!* x : Higgs mass in GeV
!* fit: strictly valid in range [90:1000], deviations from data below 0.9%
!*      slight extrapolation to range [70:1020] is allowed.
!* 27/4/2011, Oliver Brein (revised: 29/4/2011)
!******************************
        implicit none
        real*8 x,log_lhc7_rH_bb
        real*8 a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1
        real*8 a2,b2,d2,f2,h2,i2,k2
        real*8 a3,b3,d3,f3,h3,i3,k3
        real*8 a4,b4,c4,d4,e4,f4,g4,h4,i4,j4,k4

        i1 = 85204695.5338915d0
        j1 = -6912182.64018204d0
        k1 = 243962.07649863d0
        a1 = -4928.02789274431d0
        b1 = 62.909070611551d0
        c1 = -0.531357197253766d0
        d1 = 0.00300387026858175d0
        e1 = -1.12321432811065d-05
        f1 = 2.66224722543116d-08
        g1 = -3.61842683083223d-11
        h1 = 2.14600362136921d-14

        i2 = 216.998619080247d0
        k2 = -10455.1340076185d0
        a2 = 126.132581948455d0
        b2 = -0.486714890764453d0
        d2 = 2.66710181968047d-06
        f2 = -1.15648133355064d-11
        h2 = 2.12621914385691d-17

        i3 = 382.000632359183d0
        k3 = 20351.3699811471d0
        a3 = -151.281993851824d0
        b3 = 0.30675527455342d0
        d3 = -5.95717005189322d-07
        f3 = 9.30998363879783d-13
        h3 = -6.45177335713779d-19

        i4 = -84200046211.4068d0
        j4 = 1300614298.56323d0
        k4 = -8896220.2344996d0
        a4 = 35502.2157503435d0
        b4 = -91.6156640905987d0
        c4 = 0.159636867717387d0
        d4 = -0.000190242945344468d0
        e4 = 1.53100237628858d-07
        f4 = -7.96212154113225d-11
        g4 = 2.41607356489374d-14
        h4 = -3.24794747858385d-18

        log_lhc7_rH_bb=0d0

        if(x .lt. 70d0) then
        write(*,*)'function lhc7_rH_bb might not be a good fit (m_H < 70 GeV)'
        endif
        if((x .ge. 70d0) .and. (x .le. 250d0)) then
        log_lhc7_rH_bb=i1/x**3+j1/x**2+k1/x+a1+b1*x+c1*x**2+d1*x**3+e1*x**4+f1*x**5+g1*x**6+h1*x**7
        endif
        if((x .gt. 250d0) .and. (x .le. 340d0)) then
        log_lhc7_rH_bb=i2/x**3+k2/x+a2+b2*x+d2*x**3+f2*x**5+h2*x**7

        endif
        if((x .gt. 340d0) .and. (x .le. 400d0)) then
        log_lhc7_rH_bb=i3/x**3+k3/x+a3+b3*x+d3*x**3+f3*x**5+h3*x**7
        endif
        if((x .gt. 400d0) .and. (x .le. 1020d0)) then
        log_lhc7_rH_bb=i4/x**3+j4/x**2+k4/x+a4+b4*x+c4*x**2+d4*x**3+e4*x**4+f4*x**5+g4*x**6+h4*x**7
        endif
        if(x .gt. 1020d0) then
        write(*,*)'function lhc7_rH_bb might not be a good fit (m_H > 1020 GeV)'
        endif
        
        lhc7_rH_bb=exp(log_lhc7_rH_bb)
        
        end function
        
