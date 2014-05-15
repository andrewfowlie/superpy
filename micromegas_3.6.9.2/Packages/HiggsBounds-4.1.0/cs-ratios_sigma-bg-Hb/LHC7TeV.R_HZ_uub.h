!******************************************************
        real*8 function lhc7_rHZ_uub(x)
!******************************************************
!* Used LO calculation and NNLO QCD and EW K-factors.
!* This reproduces the numbers of arXiv:1101.0593 [hep-ph] 
!* to better than 1%. (PDF: MSTW 2008 NNLO)
!******************************************************
!* x : Higgs mass in GeV
!* fit: strictly valid in range [90:300], deviations from data below 0.5%
!*      slight extrapolation to range [80:310] is allowed.
!* 27/4/2011, Oliver Brein (revised: 29/4/2011)
!******************************
        implicit none
        real*8 a,b,c,d,x

        a = 0.486454541874656
        b = -9.78294961659717e-05
        c = 3.63077314490425e-07
        d = 1.62568809034006e-09

        lhc7_rHZ_uub=0d0
                                                                        
        if(x .lt. 80d0) then
        write(*,*)'function lhc7_rHZ_uub might not be a good fit (m_H < 80 GeV)'
        endif
        if((x .ge. 80d0) .and. (x .le. 310d0)) then
        lhc7_rHZ_uub=a+b*x+c*x**2+d*x**3
        endif
        if(x .gt. 310d0) then
        write(*,*)'function lhc7_rHZ_uub might not be a good fit (m_H > 310 GeV)'
        endif
                                                                                                                                        
        end function
        
