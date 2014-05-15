!******************************************************
        real*8 function tev_rHWpm_scb(x)
!******************************************************
!* PDF used: MRST NNLO 2006
!* K-factors used: drop out
!******************************************************
!c x : Higgs mass
!c fit: valid in range [50:400], deviations from data below 0.4%
        implicit none
        real*8 a1,b1,c1,d1,e1,f1,g1,h1,k,x,log_tev_rHWpm_scb

        k = 1.37596672687619
        a1 = -3.73090688631112
        b1 = -0.0145906687520353
        c1 = 6.08966130970616e-06
        d1 = 2.07499977017497e-07
        e1 = -1.31386656808655e-09
        f1 = 4.04569655822192e-12
        g1 = -6.3892305976534e-15
        h1 = 4.09548889836751e-18
!        FIT_CONVERGED = 1
!        FIT_NDF = 62
!        FIT_STDFIT = 0.00143662357896548
!        FIT_WSSR = 0.000127961013073654

        if(x .lt. 50d0) then
        write(*,*)'function tev_rHWpm_scb might not be a good fit (m_H < 50 GeV)'
        endif
        if((x .ge. 50d0) .and. (x .le. 400d0)) then
        log_tev_rHWpm_scb=k/x+a1+b1*x+c1*x**2+d1*x**3+e1*x**4+f1*x**5+g1*x**6+h1*x**7
        tev_rHWpm_scb=exp(log_tev_rHWpm_scb)
        endif
        if(x .gt. 400d0) then
        write(*,*)'function tev_rHWpm_scb might not be a good fit (m_H > 400 GeV)'
        endif

                                                                                                                        
                                                                                                                                        
        end function
