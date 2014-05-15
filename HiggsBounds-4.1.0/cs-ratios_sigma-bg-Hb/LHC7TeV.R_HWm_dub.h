!******************************************************
        real*8 function lhc7_rHWm_dub(x)
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
        real*8 a1,b1,c1,d1,x,log_lhc7_rHWm_dub

        a1 = -1.08712509764222
        b1 = -0.000284899649133387
        c1 = -8.69215640500941e-07
        d1 = 9.6651154337596e-10
                                                                        
        lhc7_rHWm_dub=0d0
        
        if(x .lt. 80d0) then
        write(*,*)'function lhc7_rHWm_dub might not be a good fit (m_H < 80 GeV)'
        endif
        if((x .ge. 80d0) .and. (x .le. 310d0)) then
        log_lhc7_rHWm_dub=a1+b1*x+c1*x**2+d1*x**3
        lhc7_rHWm_dub=exp(log_lhc7_rHWm_dub)
        endif
        if(x .gt. 310d0) then
        write(*,*)'function lhc7_rHWm_dub might not be a good fit (m_H > 310 GeV)'
        endif
                                                                                                                        
                                                                                                                                        
        end function
        
