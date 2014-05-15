!******************************************************
        real*8 function tev_rHWpm_csb(x)
!******************************************************
!* PDF used: MRST NNLO 2006
!* K-factors used: drop out
!******************************************************
!c x : Higgs mass
!c fit: valid in range [50:400], deviations from data below 0.3%
        implicit none
        real*8 a1,b1,c1,d1,e1,f1,g1,h1,k,x,log_tev_rHWpm_csb

        k = 2.69004047335856
        a1 = -3.79140986380442
        b1 = -0.0134166875946903
        c1 = -6.06503772384786e-06
        d1 = 2.7368462602839e-07
        e1 = -1.48020344245938e-09
        f1 = 4.11062191957571e-12
        g1 = -5.94166617139696e-15
        h1 = 3.51944113090225e-18
!        FIT_CONVERGED = 1
!        FIT_NDF = 62
!        FIT_STDFIT = 0.00122967200049379
!        FIT_WSSR = 9.37497801855009e-05

        if(x .lt. 50d0) then
        write(*,*)'function tev_rHWpm_csb might not be a good fit (m_H < 50 GeV)'
        endif
        if((x .ge. 50d0) .and. (x .le. 400d0)) then
        log_tev_rHWpm_csb=k/x+a1+b1*x+c1*x**2+d1*x**3+e1*x**4+f1*x**5+g1*x**6+h1*x**7
        tev_rHWpm_csb=exp(log_tev_rHWpm_csb)
        endif
        if(x .gt. 400d0) then
        write(*,*)'function tev_rHWpm_csb might not be a good fit (m_H > 400 GeV)'
        endif

                                                                                                                        
                                                                                                                                        
        end function
