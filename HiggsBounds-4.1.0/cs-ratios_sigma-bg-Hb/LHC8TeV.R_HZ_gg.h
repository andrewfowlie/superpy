!******************************************************
        real*8 function lhc8_rHZ_gg(x)
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
!* 24/9/2012, Oliver Brein / Oscar Stål
!******************************
        implicit none
        real*8 x
        real*8 a1,b1,c1,d1,e1,f1,g1,h1,k1
        real*8 a2,b2,c2,d2,e2,f2,g2,h2,k2
        real*8 a0,b0,c0,d0
        real*8 a4,b4,c4,d4
        
        a0 = 0.0307601672019924D0
        b0 = -0.000941689922892115D0
        c0 = 1.23163054836611d-05
        d0 = -3.36444019453493d-08

        a4 = -0.00838046716107623D0
        b4 = 0.00170826630713703D0
        c4 = -8.89904234127779D-06
        d4 = 1.20856588898756D-08

        k1 = 55518.4120509902D0
        a1 = -3539.01145059983D0
        b1 = 97.9427820915078D0
        c1 = -1.53707397369298D0
        d1 = 0.0149615537230227D0
        e1 = -9.24970920201235D-05
        f1 = 3.54704863046333D-07
        g1 = -7.71434487263541D-10
        h1 = 7.28563152809408D-13

        k2 = 89760.9888063226D0
        a2 = -3247.56824671049D0
        b2 = 51.063336784092D0
        c2 = -0.455752015226749D0
        d2 = 0.00252551125920739D0
        e2 = -8.89744038247416D-06
        f2 = 1.94613303763967D-08
        g2 = -2.4163848330695D-11
        h2 = 1.30404687186049D-14

        lhc8_rHZ_gg=0d0
        
        if(x .lt. 80d0) then
        write(*,*)'function lhc8_rHZ_gg might not be a good fit (m_H < 80 GeV)'
        endif
        if((x .ge. 80d0) .and. (x .lt. 90d0)) then
        lhc8_rHZ_gg=a0+b0*x+c0*x**2+d0*x**3
        endif
        if((x .ge. 90d0) .and. (x .lt. 160d0)) then
        lhc8_rHZ_gg=k1/x+a1+b1*x+c1*x**2+d1*x**3+e1*x**4+f1*x**5+g1*x**6+h1*x**7
        endif
        if((x .ge. 160d0) .and. (x .le. 300d0)) then
        lhc8_rHZ_gg=k2/x+a2+b2*x+c2*x**2+d2*x**3+e2*x**4+f2*x**5+g2*x**6+h2*x**7
        endif
        if((x .gt. 300d0) .and. (x .le. 310d0)) then
        lhc8_rHZ_gg=a4+b4*x+c4*x**2+d4*x**3
        endif
        if(x .gt. 310d0) then
        write(*,*)'function lhc8_rHZ_gg might not be a good fit (m_H > 310 GeV)'
        endif
                                                                                                                                                                                                                                                                
        end function
        
