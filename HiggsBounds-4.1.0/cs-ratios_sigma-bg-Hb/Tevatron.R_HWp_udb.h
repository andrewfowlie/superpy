!******************************************************
        real*8 function tev_rHWpm_udb(x)
!******************************************************
!* PDF used: MRST NNLO 2006
!* K-factors used: drop out
!******************************************************
!c x : Higgs mass
!c fit: valid in range [50:400], deviations from data below 0.2%
        implicit none
        real*8 a1,b1,c1,d1,e1,f1,g1,h1,k,x,log_tev_rHWpm_udb

        k = -3.04939673975276
        a1 = -0.556924228377434
        b1 = -0.00378141534901606
        c1 = 5.14697730322126e-05
        d1 = -3.90998743435728e-07
        e1 = 1.74369386764978e-09
        f1 = -4.53051110686384e-12
        g1 = 6.34012832550451e-15
        h1 = -3.68905631076079e-18
!        FIT_CONVERGED = 1
!        FIT_NDF = 62
!        FIT_STDFIT = 0.000694595363424964
!        FIT_WSSR = 2.99126885712703e-05

        if(x .lt. 50d0) then
        write(*,*)'function tev_rHWpm_udb might not be a good fit (m_H < 50 GeV)'
        endif
        if((x .ge. 50d0) .and. (x .le. 400d0)) then
        log_tev_rHWpm_udb=k/x+a1+b1*x+c1*x**2+d1*x**3+e1*x**4+f1*x**5+g1*x**6+h1*x**7
        tev_rHWpm_udb=exp(log_tev_rHWpm_udb)
        endif
        if(x .gt. 400d0) then
        write(*,*)'function tev_rHWpm_udb might not be a good fit (m_H > 400 GeV)'
        endif


                                                                                                                        
                                                                                                                                        
        end function

