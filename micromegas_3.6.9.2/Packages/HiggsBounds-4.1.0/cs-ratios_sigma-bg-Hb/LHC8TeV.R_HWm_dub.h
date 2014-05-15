!******************************************************
        real*8 function lhc8_rHWm_dub(x)
!******************************************************
!* Used LO calculation and NNLO QCD and EW K-factors.
!* This reproduces the numbers of arXiv:1101.0593 [hep-ph] 
!* to better than 1%. (PDF: MSTW 2008 NNLO)
!******************************************************
!* x : Higgs mass in GeV
!* fit: strictly valid in range [90:300], deviations from data below 0.3%
!*      slight extrapolation to range [80:310] is allowed.
!* 27/4/2011, Oliver Brein (revised: 29/4/2011)
!******************************
        implicit none
        real*8 a1,b1,c1,d1,x,log_lhc8_rHWm_dub

        a1 = -1.087125097642D0
        b1 = -0.000284899649D0
        c1 = -8.692156405009D-07
        d1 = 9.6651154337596D-10
                                                                        
        lhc8_rHWm_dub=0d0
        
        if(x .lt. 80d0) then
        write(*,*)'function lhc8_rHWm_dub might not be a good fit (m_H < 80 GeV)'
        endif
        if((x .ge. 80d0) .and. (x .le. 310d0)) then
        log_lhc8_rHWm_dub=a1+b1*x+c1*x**2+d1*x**3
        lhc8_rHWm_dub=exp(log_lhc8_rHWm_dub)
        endif
        if(x .gt. 310d0) then
        write(*,*)'function lhc8_rHWm_dub might not be a good fit (m_H > 310 GeV)'
        endif
                                                                                                                        
                                                                                                                                        
        end function
        
