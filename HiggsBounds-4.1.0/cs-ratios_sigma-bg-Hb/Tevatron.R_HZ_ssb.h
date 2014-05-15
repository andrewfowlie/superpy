!******************************************************
        real*8 function tev_rHZ_ssb(x)
!******************************************************
!* PDF used: MRST NNLO 2006
!* K-factors used: drop out
!******************************************************
!c x : Higgs mass
!c fit: valid in range [50:400], deviations from data below 0.4%
        implicit none
        real*8 a1,b1,c1,d1,e1,f1,g1,h1,k,x,log_tev_rHZ_ssb
        
        k = -12.6606626703327
        a1 = -2.44637299116304
        b1 = -0.0337220374701087
        c1 = 0.00024875850507805
        d1 = -1.59530506283027e-06
        e1 = 6.52910331412345e-09
        f1 = -1.59380826107913e-11
        g1 = 2.11709366853816e-14
        h1 = -1.17774825576315e-17
!        FIT_CONVERGED = 1
!        FIT_NDF = 62
!        FIT_STDFIT = 0.00108021511022444
!        FIT_WSSR = 7.2345610430146e-05

        if(x .lt. 50d0) then
        write(*,*)'function tev_rHZ_ssb might not be a good fit (m_H < 50 GeV)'
        endif
        if((x .ge. 50d0) .and. (x .le. 400d0)) then
        log_tev_rHZ_ssb=k/x+a1+b1*x+c1*x**2+d1*x**3+e1*x**4+f1*x**5+g1*x**6+h1*x**7
        tev_rHZ_ssb=exp(log_tev_rHZ_ssb)
        endif
        if(x .gt. 400d0) then
        write(*,*)'function tev_rHZ_ssb might not be a good fit (m_H > 400 GeV)'
        endif


                                                                                                                        
                                                                                                                                        
        end function
