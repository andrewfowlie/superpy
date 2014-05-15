!******************************************************
        real*8 function tev_rHZ_uub(x)
!******************************************************
!* PDF used: MRST NNLO 2006
!* K-factors used: drop out
!******************************************************
!c x : Higgs mass
!c fit: valid in range [50:400], deviations from data below 0.08%
        implicit none
        real*8 a,b,c,d,e,f,g,h,i,j,k,x

        i = -71574.5518426517
        j = 5160.62795328258
        k = -154.217142570514
        a = 3.08939823563244
        b = -0.022276581937085
        c = 0.000136990393425251
        d = -5.05349956146734e-07
        e = 1.06940967044471e-09
        f = -1.06545623236754e-12
        g = 3.40159712939385e-17
        h = 5.32813313631463e-19
!        FIT_CONVERGED = 1
!        FIT_NDF = 60
!        FIT_STDFIT = 0.000224402040673476
!        FIT_WSSR = 3.02137655150523e-06

        if(x .lt. 50d0) then
        write(*,*)'function tev_rHZ_uub might not be a good fit (m_H < 50 GeV)'
        endif
        if((x .ge. 50d0) .and. (x .le. 400d0)) then
        tev_rHZ_uub=i/x**3+j/x**2+k/x+a+b*x+c*x**2+d*x**3+e*x**4+f*x**5+g*x**6+h*x**7
        endif
        if(x .gt. 400d0) then
        write(*,*)'function tev_rHZ_uub might not be a good fit (m_H > 400 GeV)'
        endif


                                                                                                                                        
        end function

