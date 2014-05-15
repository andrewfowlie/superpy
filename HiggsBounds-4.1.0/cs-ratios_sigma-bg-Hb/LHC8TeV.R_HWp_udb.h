!******************************************************
        real*8 function lhc8_rHWp_udb(x)
!******************************************************
!* Used LO calculation and NNLO QCD and EW K-factors.
!* This reproduces the numbers of arXiv:1101.0593 [hep-ph] 
!* to better than 1%. (PDF: MSTW 2008 NNLO)
!******************************************************
!* x : Higgs mass in GeV
!* fit: strictly valid in range [90:300], deviations from data below 0.2%
!*      slight extrapolation to range [80:310] is allowed.
!* 24/9/2012, Oliver Brein / Oscar Stål
!******************************
        implicit none
        real*8 a1,b1,c1,d1,x,log_lhc8_rHWp_udb

        a1 = -0.6484690466396D0
        b1 = 0.00136216142736D0
        c1 = -2.6547629079855D-06
        d1 = 2.57155565269944D-09

        lhc8_rHWp_udb=0d0

        if(x .lt. 80d0) then
        write(*,*)'function lhc8_rHWp_udb might not be a good fit (m_H < 80 GeV)'
        endif
        if((x .ge. 80d0) .and. (x .le. 310d0)) then
        log_lhc8_rHWp_udb=a1+b1*x+c1*x**2+d1*x**3
        lhc8_rHWp_udb=exp(log_lhc8_rHWp_udb)
        endif
        if(x .gt. 310d0) then
        write(*,*)'function lhc8_rHWp_udb might not be a good fit (m_H > 310 GeV)'
        endif
                                                                                                                                        
        end function
        
