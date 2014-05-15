!******************************************************
        real*8 function tev_rH_gg(x)
!******************************************************
!* PDF used: MRST NNLO 2006
!* K-factors used: K_ggH_TEV = 3.6
!******************************************************
!c x : Higgs mass
!c fit: valid in range [50:400], deviations from data below 0.04%
        implicit none
        real*8 a,b,c,d,e,f,g,h,i,j,k,x

       i = 6643.16353482064
       j = -472.482085531536
       k = 15.3118835634636
       a = 0.64668414135986
       b = 0.00450601858533845
       c = -3.58931308628694e-05
       d = 1.84446994484542e-07
       e = -6.06967613406461e-10
       f = 1.22838296555632e-12
       g = -1.38453592806246e-15
       h = 6.61589556779331e-19
!       FIT_CONVERGED = 1
!       FIT_NDF = 340
!       FIT_STDFIT = 9.30844138157386e-06
!       FIT_WSSR = 2.94600075244269e-08

        if(x .lt. 50d0) then
        write(*,*)'function tev_rH_gg might not be a good fit (m_H < 50 GeV)'
        endif
        if((x .ge. 50d0) .and. (x .le. 400d0)) then
        tev_rH_gg=i/x**3+j/x**2+k/x+a+b*x+c*x**2+d*x**3+e*x**4+f*x**5+g*x**6+h*x**7
        endif
        if(x .gt. 400d0) then
        write(*,*)'function tev_rH_gg might not be a good fit (m_H > 400 GeV)'
        endif
                                                                                                                                        
        end function
