!******************************************************
        real*8 function tev_cs_bg_Hb_c1_SM(x)
!******************************************************
!* Fit for sigma(pp-bar -> b g -> H b) in pikobarn as a function of Higgs mass in GeV
!* However, sigma(pp-bar -> b g -> H b)+sigma(pp-bar -> b-bar g -> H b-bar) is 
!* needed here. Hence, the factor of 2 below.
!*
!* PDF used: MRST NNLO 2006
!* K-factors used: none. 
!* However, used: a) scale choice MU_R=MU_F=M_h/4
!*                b) 1-loop running bottom mass instead of pole mass
!* This resembles pretty closely the NLO result of Maltoni et al. (hep-ph/0204093)
!* cuts used: ptmin=15 GeV, etamax=2 (c1)
!******************************************************
!c x : Higgs mass
!c fit: valid in range [10:400], deviations from data below 0.6%
        implicit none
        real*8 a,b,c,d,e,f,g,h,i,j,k,x,log_tev_cs_bg_Hb_c1

       i = 979.232120836288
       j = -264.696290260115
       k = 23.7789096288902
       a = -0.976837806470803
       b = -0.0694778939102754
       c = 0.000300026595336677
       d = -1.31248432594797e-06
       e = 4.11789388677158e-09
       f = -8.31929589817916e-12
       g = 9.55612036163061e-15
       h = -4.71205719524022e-18
!       FIT_CONVERGED = 1
!       FIT_NDF = 68
!       FIT_STDFIT = 0.00192144957958659
!       FIT_WSSR = 0.000251053857108756

       log_tev_cs_bg_Hb_c1=i/x**3+j/x**2+k/x+a+b*x+c*x**2+d*x**3+e*x**4+f*x**5+g*x**6+h*x**7
       tev_cs_bg_Hb_c1_SM=2d0*exp(log_tev_cs_bg_Hb_c1)

       end function

