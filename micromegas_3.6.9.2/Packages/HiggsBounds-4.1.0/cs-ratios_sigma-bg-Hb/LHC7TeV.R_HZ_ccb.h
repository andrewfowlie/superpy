!******************************************************
        real*8 function lhc7_rHZ_ccb(x)
!******************************************************
!* Used LO calculation and NNLO QCD and EW K-factors.
!* This reproduces the numbers of arXiv:1101.0593 [hep-ph] 
!* to better than 1%. (PDF: MSTW 2008 NNLO)
!******************************************************
!* x : Higgs mass in GeV
!* fit: strictly valid in range [90:300], deviations from data below 0.4%
!*      slight extrapolation to range [80:310] is allowed.
!* 27/4/2011, Oliver Brein (revised: 29/4/2011)
!******************************
        implicit none
        real*8 a1,b1,c1,d1,e1,x,log_lhc7_rHZ_ccb
        
        a1 = -3.24950342921778
        b1 = -0.00221218122605409
        c1 = -3.03932963115695e-05
        d1 = 1.36905436364982e-07
        e1 = -1.80512708031642e-10

        lhc7_rHZ_ccb=0d0
                                                                                
        if(x .lt. 80d0) then
        write(*,*)'function lhc7_rHZ_ccb might not be a good fit (m_H < 80 GeV)'
        endif
        if((x .ge. 80d0) .and. (x .le. 310d0)) then
        log_lhc7_rHZ_ccb=a1+b1*x+c1*x**2+d1*x**3+e1*x**4
        lhc7_rHZ_ccb=exp(log_lhc7_rHZ_ccb)
        endif
        if(x .gt. 310d0) then
        write(*,*)'function lhc7_rHZ_ccb might not be a good fit (m_H > 310 GeV)'
        endif        

                                                                                                                        
                                                                                                                                        
        end function
        
