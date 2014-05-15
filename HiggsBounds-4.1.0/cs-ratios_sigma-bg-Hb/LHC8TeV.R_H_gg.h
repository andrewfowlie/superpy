!******************************************************
        real*8 function lhc8_rH_gg(x)
!******************************************************
!* ggH process : consensus numbers used (PDF: MSTW 2008 NNLO)
!*		 see: arXiv:1101.0593 [hep-ph] 
!* 		 numbers taken from: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt7TeV
!* bbH process : numbers in range [70:1020] generated with bbh@nnlo v. 1.3 (PDF: MSTW 2008 NNLO) 
!* 		 with input parameters according to arXiv:1101.0593 [hep-ph]
!******************************************************
!* x : Higgs mass in GeV
!* fit: strictly valid in range [90:1000], deviations from data below 0.012%
!*      slight extrapolation to range [70:1020] is allowed.
!* 5/9/2012, Oliver Brein (implemented by OS)
!******************************
        implicit none
        real*8 a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1,x
        real*8 a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,k2

        i1 = 964113.343439237D0
        j1 = -63148.8079879927D0
        k1 = 1796.32834256942D0
        a1 = -28.2111518227145D0
        b1 = 0.299213570678883D0
        c1 = -0.00201616549445028D0
        d1 = 9.05199533157006D-06
        e1 = -2.67657331056638D-08
        f1 = 4.99699199024971D-11
        g1 = -5.3303608557019D-14
        h1 = 2.47284262943397D-17
        
        i2 = 45311966.2179785D0
        j2 = -756076.101453613D0
        k2 = 5503.92080433926D0
        a2 = -22.0725813195331D0
        b2 = 0.0617412845820624D0
        c2 = -0.000110242205100537D0
        d2 = 1.32982976036479D-07
        e2 = -1.06924880989388D-10
        f2 = 5.47593956447315D-14
        g2 = -1.60873930927401D-17
        h2 = 2.04985949040854D-21

        lhc8_rH_gg=0d0

        if(x .lt. 70d0) then
        write(*,*)'function lhc8_rH_gg might not be a good fit (m_H < 70 GeV)'
        endif
        if((x .ge. 70d0) .and. (x .le. 350d0)) then
        lhc8_rH_gg=i1/x**3+j1/x**2+k1/x+a1+b1*x+c1*x**2+d1*x**3+e1*x**4+f1*x**5+g1*x**6+h1*x**7
        endif
        if((x .ge. 350d0) .and. (x .le. 1020d0)) then
        lhc8_rH_gg=i2/x**3+j2/x**2+k2/x+a2+b2*x+c2*x**2+d2*x**3+e2*x**4+f2*x**5+g2*x**6+h2*x**7
        endif
        if(x .gt. 1020d0) then
        write(*,*)'function lhc8_rH_gg might not be a good fit (m_H > 1020 GeV)'
        endif
                                                                                                                                        
        end function
        
