!******************************************************
        real*8 function tev_cs_bg_Hb_SM(x)
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
!******************************************************
! x : Higgs mass in GeV.
! fit: valid in range [10:400], deviations from data below 0.3%
        implicit none
        real*8 x,a,b,c,d,e,f,g,h,i,j,k,log_tev_cs_bg_Hb

        i = 1123.83787799968
        j = -356.707392839726
        k = 45.9128876730515
        a = 1.6356919815081
        b = -0.0936888262673351
        c = 0.000517589926136613
        d = -2.60428788311317e-06
        e = 9.05347214042015e-09
        f = -1.99410338384936e-11
        g = 2.48409692520096e-14
        h = -1.32760959238599e-17
!        FIT_CONVERGED = 1
!        FIT_NDF = 68
!        FIT_STDFIT = 0.000936864903217008
!        FIT_WSSR = 5.96846775878273e-05

        log_tev_cs_bg_Hb=i/x**3+j/x**2+k/x+a+b*x+c*x**2+d*x**3+e*x**4+f*x**5+g*x**6+h*x**7
        tev_cs_bg_Hb_SM=2d0*exp(log_tev_cs_bg_Hb)
                                                                                                                        
                                                                                                                                        
        end function

