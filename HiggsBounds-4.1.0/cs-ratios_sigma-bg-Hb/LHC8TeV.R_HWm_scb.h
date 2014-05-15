!******************************************************
        real*8 function lhc8_rHWm_scb(x)
!******************************************************
!* Used LO calculation and NNLO QCD and EW K-factors.
!* This reproduces the numbers of arXiv:1101.0593 [hep-ph] 
!* to better than 1%. (PDF: MSTW 2008 NNLO)
!******************************************************
!* x : Higgs mass in GeV
!* fit: strictly valid in range [90:300], deviations from data below 0.3%
!*      slight extrapolation to range [80:310] is allowed.
!* 24/9/2012, Oliver Brein / Oscar Stål
!******************************
        implicit none
        real*8 a1,b1,c1,d1,e1,x,log_lhc8_rHWm_scb

        a1 = -2.620114254755D0
        b1 = -0.004760506544D0
        c1 = 6.9165415989983D-07
        d1 = 1.2915247823349D-08
        e1 = -2.264120499263D-11

        lhc8_rHWm_scb=0d0

        if(x .lt. 80d0) then
        write(*,*)'function lhc8_rHWm_scb might not be a good fit (m_H < 80 GeV)'
        endif
        if((x .ge. 80d0) .and. (x .le. 310d0)) then
        log_lhc8_rHWm_scb=a1+b1*x+c1*x**2+d1*x**3+e1*x**4
        lhc8_rHWm_scb=exp(log_lhc8_rHWm_scb)
        endif
        if(x .gt. 310d0) then
        write(*,*)'function lhc8_rHWm_scb might not be a good fit (m_H > 310 GeV)'
        endif

        end function
        
