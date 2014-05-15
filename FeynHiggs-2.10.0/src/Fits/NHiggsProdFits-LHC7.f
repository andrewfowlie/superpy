
	if( sqrtm .lt. 9.486832980505138D0 )
     &    Warning('Extrapolating tthSM('//Digit(h)//') in MHiggs')

        tthSM(h) = exp(12.150198749224147D0 + 
     &     sqrtm*(-0.8247169745199131D0 + 0.012192698761097822D0*sqrtm)
     &     )

#ifdef DETAILED_DEBUG
	DPROD "tthSM(h) =", tthSM(h) ENDL
#endif


	if( sqrtm .lt. 9.486832980505138D0 )
     &    Warning('Extrapolating qqhSM('//Digit(h)//') in MHiggs')

        qqhSM(h) = exp(9.556644671514436D0 + 
     &     sqrtm*(-0.21542683374250993D0 - 
     &        0.0003657746642486446D0*sqrtm))

#ifdef DETAILED_DEBUG
	DPROD "qqhSM(h) =", qqhSM(h) ENDL
#endif


	if( sqrtm .lt. 9.486832980505138D0 )
     &    Warning('Extrapolating WhSM('//Digit(h)//') in MHiggs')

        WhSM(h) = exp(14.613444804667251D0 + 
     &     sqrtm*(-0.8671330612724978D0 + 0.01141699023461514D0*sqrtm))

#ifdef DETAILED_DEBUG
	DPROD "WhSM(h) =", WhSM(h) ENDL
#endif


	if( sqrtm .lt. 9.486832980505138D0 )
     &    Warning('Extrapolating ZhSM('//Digit(h)//') in MHiggs')

        ZhSM(h) = exp(13.03878791914205D0 + 
     &     sqrtm*(-0.7269995240810023D0 + 0.006743670763495621D0*sqrtm)
     &     )

#ifdef DETAILED_DEBUG
	DPROD "ZhSM(h) =", ZhSM(h) ENDL
#endif

