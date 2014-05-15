
	if( m .lt. 110.D0 )
     &    Warning('Extrapolating GammaSM(H0VV('//Digit(h)//',3)) in MHiggs')

	if( m.lt.161.D0 ) then

        GammaSM(H0VV(h,3)) = 
     &  exp(-65.8900129090698D0 + 
     &    m*(0.9688132244049187D0 + 
     &       m*(-0.005607462162886842D0 + 0.000011893449991338892D0*m))
     &    )

	else

        GammaSM(H0VV(h,3)) = 
     &  exp((3.735823586939495D10 + 
     &      m*(-2.2752344484718385D9 + 
     &         m*(2.853110882358539D7 + 
     &            m*(-132662.46162488303D0 + 210.65166534520054D0*m))))
     &     /(-5.753715195825776D10 + 
     &      m*(9.98748114009375D8 + 
     &         m*(-5.084653579394836D6 + 
     &            m*(3473.7525214318503D0 + 21.288377319851392D0*m)))))

	endif

#ifdef DETAILED_DEBUG
	DPROD 'GammaSM(H0VV(h,3)) =', GammaSM(H0VV(h,3)) ENDL
#endif


	if( m .lt. 110.D0 )
     &    Warning('Extrapolating GammaSM(H0VV('//Digit(h)//',4)) in MHiggs')

	if( m.lt.161.D0 ) then

        GammaSM(H0VV(h,4)) = 
     &  exp(-112.31137784841D0 + 
     &    m*(2.1967220647228136D0 + 
     &       m*(-0.015783982754797143D0 + 0.00003955767499872661D0*m)))

	else

        GammaSM(H0VV(h,4)) = 
     &  exp((2.243397065250441D9 + 
     &      m*(-3.008339624796666D7 + 
     &         m*(124805.53107759365D0 + 
     &            m*(-184.6192174620866D0 + 0.1617913046889949D0*m))))/
     &    (-2.095596108520445D8 + 
     &      m*(6409.693296513176D0 + 
     &         m*(10859.356293770943D0 + 
     &            m*(-18.11269530580481D0 + 0.020676545445741868D0*m)))
     &      ))

	endif

#ifdef DETAILED_DEBUG
	DPROD 'GammaSM(H0VV(h,4)) =', GammaSM(H0VV(h,4)) ENDL
#endif

