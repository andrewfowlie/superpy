!******************************************************
        real*8 function lhc7_rHZ_gg(x)
!******************************************************
!* Used LO calculation and NNLO QCD and EW K-factors.
!* This reproduces the numbers of arXiv:1101.0593 [hep-ph] 
!* to better than 1%. (PDF: MSTW 2008 NNLO)
!******************************************************
!* x : Higgs mass in GeV
!* fit: strictly valid in range [90:300], deviations from data below 0.4%
!*      slight extrapolation to range [80:310] is allowed.
!* /!\: For the fit it is important that the 
!*      coefficients are given explicitly with a double precision "d".
!*      Otherwise the fit result would be off by a few percent!
!*      It seems, there are strong cancellations occuring in this
!*      expressions.
!* 27/4/2011, Oliver Brein (revised: 29/4/2011)
!******************************
        implicit none
        real*8 x
        real*8 a1,b1,c1,d1,e1,f1,g1,h1,k1
        real*8 a2,b2,c2,d2,e2,f2,g2,h2,k2
        real*8 a0,b0,c0,d0
        real*8 a4,b4,c4,d4
        
        a0 = 0.0191763343414287d0
        b0 = -0.00061646800365048d0
        c0 = 8.86573978886801d-06
        d0 = -2.36430682752164d-08

        a4 = -0.123939763255212d0
        b4 = 0.00282168359675478d0
        c4 = -1.27833948672529d-05
        d4 = 1.67990251927367d-08

        k1 = 8829.3243692852d0
        a1 = -598.834333600709d0
        b1 = 17.6794228692499d0
        c1 = -0.296767828181659d0
        d1 = 0.00309822063906165d0
        e1 = -2.06007944611655d-05
        f1 = 8.52084806171817d-08
        g1 = -2.00468391301571d-10
        h1 = 2.0542042537335d-13

        k2 = 58641.408114179d0
        a2 = -2112.32446415142d0
        b2 = 33.0718507552994d0
        c2 = -0.293959260853248d0
        d2 = 0.00162247924882443d0
        e2 = -5.69398492892907d-06
        f2 = 1.24073925811521d-08
        g2 = -1.53480930354464d-11
        h2 = 8.25235966239085d-15

        lhc7_rHZ_gg=0d0
        
        if(x .lt. 80d0) then
        write(*,*)'function lhc7_rHZ_gg might not be a good fit (m_H < 80 GeV)'
        endif
        if((x .ge. 80d0) .and. (x .lt. 90d0)) then
        lhc7_rHZ_gg=a0+b0*x+c0*x**2+d0*x**3
        endif
        if((x .ge. 90d0) .and. (x .lt. 160d0)) then
        lhc7_rHZ_gg=k1/x+a1+b1*x+c1*x**2+d1*x**3+e1*x**4+f1*x**5+g1*x**6+h1*x**7
        endif
        if((x .ge. 160d0) .and. (x .le. 300d0)) then
        lhc7_rHZ_gg=k2/x+a2+b2*x+c2*x**2+d2*x**3+e2*x**4+f2*x**5+g2*x**6+h2*x**7
        endif
        if((x .gt. 300d0) .and. (x .le. 310d0)) then
        lhc7_rHZ_gg=a4+b4*x+c4*x**2+d4*x**3
        endif
        if(x .gt. 310d0) then
        write(*,*)'function lhc7_rHZ_gg might not be a good fit (m_H > 310 GeV)'
        endif
                                                                                                                                                                                                                                                                
        end function
        
