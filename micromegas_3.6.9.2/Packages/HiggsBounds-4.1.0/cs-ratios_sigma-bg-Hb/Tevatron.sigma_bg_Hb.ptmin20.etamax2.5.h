!******************************************************
        real*8 function tev_cs_bg_Hb_c2_SM(x)
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
!* cuts used: ptmin=20 GeV, etamax=2.5 (c2)
!******************************************************
!c x : Higgs mass
!c fit: valid in range [10:400], deviations from data below 0.5%
        implicit none
        real*8 a,b,c,d,e,f,g,h,i,j,k,x,log_tev_cs_bg_Hb_c2

        i = 834.839818583352
        j = -191.225997471594
        k = 11.7585430990651
        a = -1.19870317148684
        b = -0.0704359029295396
        c = 0.000329995529227224
        d = -1.56164857859755e-06
        e = 5.21024154218793e-09
        f = -1.10321498760019e-11
        g = 1.31618393769414e-14
        h = -6.7078954437179e-18
!        FIT_CONVERGED = 1
!        FIT_NDF = 68
!        FIT_STDFIT = 0.00191125120085312
!        FIT_WSSR = 0.00024839591838785

        log_tev_cs_bg_Hb_c2=i/x**3+j/x**2+k/x+a+b*x+c*x**2+d*x**3+e*x**4+f*x**5+g*x**6+h*x**7
        tev_cs_bg_Hb_c2_SM=2d0*exp(log_tev_cs_bg_Hb_c2)
                                                                                                                        
                                                                                                                                        
        end function
