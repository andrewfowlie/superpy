!******************************************************
        real*8 function tev_rHZ_ccb(x)
!******************************************************
!* PDF used: MRST NNLO 2006
!* K-factors used: drop out
!******************************************************
!c x : Higgs macc
!c fit: valid in range [50:400], deviations from data below 0.4%
        implicit none
        real*8 a1,b1,c1,d1,e1,f1,g1,h1,k,x,log_tev_rHZ_ccb
        
        k = -5.90932178929537
        a1 = -4.04714305390897
        b1 = -0.027748553051398
        c1 = 0.00017806395867611
        d1 = -1.14479567497352e-06
        e1 = 5.00189692600156e-09
        f1 = -1.31722787138353e-11
        g1 = 1.88924774486347e-14
        h1 = -1.13267553145036e-17
!        FIT_CONVERGED = 1
!        FIT_NDF = 62
!        FIT_STDFIT = 0.00134597095517168
!        FIT_WSSR = 0.000112321544354277

        if(x .lt. 50d0) then
        write(*,*)'function tev_rHZ_ccb might not be a good fit (m_H < 50 GeV)'
        endif
        if((x .ge. 50d0) .and. (x .le. 400d0)) then
        log_tev_rHZ_ccb=k/x+a1+b1*x+c1*x**2+d1*x**3+e1*x**4+f1*x**5+g1*x**6+h1*x**7
        tev_rHZ_ccb=exp(log_tev_rHZ_ccb)
        endif
        if(x .gt. 400d0) then
        write(*,*)'function tev_rHZ_ccb might not be a good fit (m_H > 400 GeV)'
        endif        

                                                                                                                        
                                                                                                                                        
        end function

