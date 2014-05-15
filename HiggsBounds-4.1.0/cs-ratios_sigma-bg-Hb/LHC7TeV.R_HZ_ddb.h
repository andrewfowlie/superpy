!******************************************************
        real*8 function lhc7_rHZ_ddb(x)
!******************************************************
!* Used LO calculation and NNLO QCD and EW K-factors.
!* This reproduces the numbers of arXiv:1101.0593 [hep-ph] 
!* to better than 1%. (PDF: MSTW 2008 NNLO)
!******************************************************
!* x : Higgs mass in GeV
!* fit: strictly valid in range [90:300], deviations from data below 0.2%
!*      slight extrapolation to range [80:310] is allowed.
!* 27/4/2011, Oliver Brein (revised: 29/4/2011)
!******************************
        implicit none
        real*8 a,b,c,d,e,x

        a = 0.375742808544361
        b = 1.20956191165674e-05
        c = 1.87943156511064e-06
        d = -1.38392054025039e-09
        e = -7.1277797851365e-11

        lhc7_rHZ_ddb=0d0

        if(x .lt. 80d0) then
        write(*,*)'function lhc7_rHZ_ddb might not be a good fit (m_H < 80 GeV)'
        endif
        if((x .ge. 80d0) .and. (x .le. 310d0)) then
        lhc7_rHZ_ddb=a+b*(x-200d0)+c*(x-200d0)**2+d*(x-200d0)**3+e*(x-200d0)**4
        endif
        if(x .gt. 310d0) then
        write(*,*)'function lhc7_rHZ_ddb might not be a good fit (m_H > 310 GeV)'
        endif

                                                                                                                                        
        end function
        
