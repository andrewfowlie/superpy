
	if( sqrtm .lt. 8.366600265340756D0 )
     &    Warning('Extrapolating bbhSM('//Digit(h)//') in MHiggs')

        bbhSM(h) = exp(12.737840556011133D0 + 
     &    sqrtm*(-1.1385193320022062D0 + 
     &       sqrtm*(0.01974128710664379D0 - 
     &          0.00037385660449064606D0*sqrtm)))

#ifdef DETAILED_DEBUG
	DPROD 'bbhSM(h) =', bbhSM(h) ENDL
#endif


	if( sqrtm .lt. 10.D0 )
     &    Warning('Extrapolating btagbhSM('//Digit(h)//') in MHiggs')

        btagbhSM(h) = exp(9.877524541618678D0 + 
     &    sqrtm*(-0.8219518503005108D0 + 0.0027226543828609293D0*sqrtm)
     &    )

#ifdef DETAILED_DEBUG
	DPROD 'btagbhSM(h) =', btagbhSM(h) ENDL
#endif


	if( sqrtm .lt. 10.D0 )
     &    Warning('Extrapolating gghSM('//Digit(h)//') in MHiggs')

	if( sqrtm.lt.18.61719635176038D0 ) then

        gghSM(h) = exp(-46.27548863223735D0 + 
     &    sqrtm*(21.05153412243787D0 + 
     &       sqrtm*(-3.114426854810264D0 + 
     &          sqrtm*(0.22259698908928452D0 + 
     &             sqrtm*(-0.007906674212839828D0 + 
     &                0.00011214075131305213D0*sqrtm)))))

	else

        gghSM(h) = exp(950745.2707612333D0 + 
     &    sqrtm*(-322501.50118803204D0 + 
     &       sqrtm*(46836.078415459095D0 + 
     &          sqrtm*(-3775.078459976367D0 + 
     &             sqrtm*(182.39054792755945D0 + 
     &                sqrtm*
     &                 (-5.282227206486968D0 + 
     &                   sqrtm*
     &                    (0.08490843739771639D0 - 
     &                      0.0005843878365243091D0*sqrtm)))))))

	endif

#ifdef DETAILED_DEBUG
	DPROD 'gghSM(h) =', gghSM(h) ENDL
#endif


	if( sqrtm .lt. 10.D0 )
     &    Warning('Extrapolating gghSMNLO('//Digit(h)//') in MHiggs')

	if( sqrtm.lt.18.61719635176038D0 ) then

        gghSMNLO(h) = exp(-62.1528748368274D0 + 
     &    sqrtm*(28.022572761357544D0 + 
     &       sqrtm*(-4.315416736512078D0 + 
     &          sqrtm*(0.3217952269508727D0 + 
     &             sqrtm*(-0.011859299816808468D0 + 
     &                0.0001732451930086779D0*sqrtm)))))

	else if( sqrtm.lt.28.284271247461902D0 ) then

        gghSMNLO(h) = exp(-18664.02649145578D0 + 
     &    sqrtm*(5478.080272797016D0 + 
     &       sqrtm*(-687.7399519226109D0 + 
     &          sqrtm*(47.88565030357969D0 + 
     &             sqrtm*(-1.9969217914097743D0 + 
     &                sqrtm*
     &                 (0.04986840045209054D0 + 
     &                   sqrtm*
     &                    (-0.0006904451025769653D0 + 
     &                      4.088186386751036D-6*sqrtm)))))))

	else

        gghSMNLO(h) = exp(0.3675877978198016D0 + 
     &    sqrtm*(1.8152620147303542D0 + 
     &       sqrtm*(-0.14515247129992012D0 + 
     &          sqrtm*(0.0038809028465699442D0 - 
     &             0.00004282737400642601D0*sqrtm))))

	endif

#ifdef DETAILED_DEBUG
	DPROD 'gghSMNLO(h) =', gghSMNLO(h) ENDL
#endif


	if( sqrtm .lt. 10.D0 )
     &    Warning('Extrapolating tthSM('//Digit(h)//') in MHiggs')

        tthSM(h) = exp(6.2156010658600405D0 + 
     &    sqrtm*(-0.29064474105268673D0 - 0.01219609615555369D0*sqrtm))

#ifdef DETAILED_DEBUG
	DPROD 'tthSM(h) =', tthSM(h) ENDL
#endif


	if( sqrtm .lt. 10.D0 )
     &    Warning('Extrapolating qqhSM('//Digit(h)//') in MHiggs')

        qqhSM(h) = exp(7.002201651011047D0 + 
     &    sqrtm*(-0.14443770718417903D0 - 0.009571109288237433D0*sqrtm)
     &    )

#ifdef DETAILED_DEBUG
	DPROD 'qqhSM(h) =', qqhSM(h) ENDL
#endif


	if( sqrtm .lt. 10.D0 )
     &    Warning('Extrapolating WhSM('//Digit(h)//') in MHiggs')

        WhSM(h) = exp(12.381567520024452D0 + 
     &    sqrtm*(-0.6881926914681327D0 + 0.0015468871671186748D0*sqrtm)
     &    )

#ifdef DETAILED_DEBUG
	DPROD 'WhSM(h) =', WhSM(h) ENDL
#endif


	if( sqrtm .lt. 10.D0 )
     &    Warning('Extrapolating ZhSM('//Digit(h)//') in MHiggs')

        ZhSM(h) = exp(11.328287585810497D0 + 
     &    sqrtm*(-0.6319216224730427D0 + 0.001036383720161046D0*sqrtm))

#ifdef DETAILED_DEBUG
	DPROD 'ZhSM(h) =', ZhSM(h) ENDL
#endif

