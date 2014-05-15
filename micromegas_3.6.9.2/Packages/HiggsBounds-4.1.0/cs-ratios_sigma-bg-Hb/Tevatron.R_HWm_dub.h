!******************************************************
        real*8 function tev_rHWpm_dub(x)
!******************************************************
!* PDF used: MRST NNLO 2006
!* K-factors used: drop out
!******************************************************
!c x : Higgs mass
!c fit: valid in range [50:400], deviations from data below 0.2%
        implicit none
        real*8 a1,b1,c1,d1,e1,f1,g1,h1,k,x,log_tev_rHWpm_dub

        k = 2.75168749742592
        a1 = -0.916249712981603
        b1 = 0.00501661388296234
        c1 = -6.09312384221491e-05
        d1 = 4.36703298238695e-07
        e1 = -1.88582796595321e-09
        f1 = 4.80568676735866e-12
        g1 = -6.64120139865788e-15
        h1 = 3.83101236728375e-18
!        FIT_CONVERGED = 1
!        FIT_NDF = 62
!        FIT_STDFIT = 0.000694750227272757
!        FIT_WSSR = 2.99260284543239e-05


        if(x .lt. 50d0) then
        write(*,*)'function tev_rHWpm_dub might not be a good fit (m_H < 50 GeV)'
        endif
        if((x .ge. 50d0) .and. (x .le. 400d0)) then
        log_tev_rHWpm_dub=k/x+a1+b1*x+c1*x**2+d1*x**3+e1*x**4+f1*x**5+g1*x**6+h1*x**7
        tev_rHWpm_dub=exp(log_tev_rHWpm_dub)
        endif
        if(x .gt. 400d0) then
        write(*,*)'function tev_rHWpm_dub might not be a good fit (m_H > 400 GeV)'
        endif
                                                                                                                        
                                                                                                                                        
        end function
