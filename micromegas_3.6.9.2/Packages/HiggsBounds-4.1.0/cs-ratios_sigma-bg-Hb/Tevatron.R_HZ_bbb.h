!******************************************************
        real*8 function tev_rHZ_bbb(x)
!******************************************************
!* PDF used: MRST NNLO 2006
!* K-factors used: drop out
!******************************************************
!c x : Higgs mabb
!c fit: valid in range [50:400], deviations from data below 0.3%
        implicit none
        real*8 a1,b1,c1,d1,e1,f1,g1,h1,k,x,log_tev_rHZ_bbb

        k = -4.17175352576335
        a1 = -4.96425924657662
        b1 = -0.0235199724168584
        c1 = 0.000104972287882997
        d1 = -5.24036826273742e-07
        e1 = 1.94151528615174e-09
        f1 = -4.41213051593084e-12
        g1 = 5.50964825163976e-15
        h1 = -2.89849178633677e-18
!        FIT_CONVERGED = 1
!        FIT_NDF = 62
!        FIT_STDFIT = 0.00118575488110775
!        FIT_WSSR = 8.71729075603935e-05
        
        if(x .lt. 50d0) then
        write(*,*)'function tev_rHZ_bbb might not be a good fit (m_H < 50 GeV)'
        endif
        if((x .ge. 50d0) .and. (x .le. 400d0)) then
        log_tev_rHZ_bbb=k/x+a1+b1*x+c1*x**2+d1*x**3+e1*x**4+f1*x**5+g1*x**6+h1*x**7
        tev_rHZ_bbb=exp(log_tev_rHZ_bbb)
        endif
        if(x .gt. 400d0) then
        write(*,*)'function tev_rHZ_bbb might not be a good fit (m_H > 400 GeV)'
        endif
        

                                                                                                                        
                                                                                                                                        
        end function

