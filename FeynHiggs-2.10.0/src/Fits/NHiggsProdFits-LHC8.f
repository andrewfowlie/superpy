
	if( sqrtm .lt. 8.94427190999916D0 )
     &    Warning('Extrapolating tthSM('//Digit(h)//') in MHiggs')

        tthSM(h) = exp(12.29188828965328D0 + 
     &     sqrtm*(-0.7960204182731393D0 + 0.011714090379440354D0*sqrtm)
     &     )

#ifdef DETAILED_DEBUG
	DPROD "tthSM(h) =", tthSM(h) ENDL
#endif


	if( sqrtm .lt. 8.94427190999916D0 )
     &    Warning('Extrapolating qqhSM('//Digit(h)//') in MHiggs')

	if( sqrtm.lt.17.029386365926403D0 ) then

        qqhSM(h) = exp(9.349914150012411D0 + 
     &     sqrtm*(-0.15568199440962405D0 - 0.00206798159522837D0*sqrtm)
     &     )

	else

        qqhSM(h) = exp(-149000.7180863036D0 + 
     &     sqrtm*(61264.66246891521D0 + 
     &        sqrtm*(-11265.229724678064D0 + 
     &           sqrtm*(1220.1402864472052D0 + 
     &              sqrtm*(-86.21877187999989D0 + 
     &                 sqrtm*
     &                  (4.153941351062099D0 + 
     &                    sqrtm*
     &                     (-0.13821348179032897D0 + 
     &                       sqrtm*
     &                       (0.003136517963886148D0 + 
     &                       sqrtm*
     &                       (-0.00004646729673954824D0 + 
     &                       sqrtm*
     &                       (4.058853519457965D-7 - 
     &                       1.5875864564765808D-9*sqrtm))))))))))

	endif

#ifdef DETAILED_DEBUG
	DPROD "qqhSM(h) =", qqhSM(h) ENDL
#endif


	if( sqrtm .lt. 8.94427190999916D0 )
     &    Warning('Extrapolating WhSM('//Digit(h)//') in MHiggs')

        WhSM(h) = exp(8.638638290870523D0 + 
     &     sqrtm*(1.0675625132840283D0 + 
     &        sqrtm*(-0.21443231616049974D0 + 
     &           sqrtm*(0.011604258534538489D0 - 
     &              0.00022048468153196195D0*sqrtm))))

#ifdef DETAILED_DEBUG
	DPROD "WhSM(h) =", WhSM(h) ENDL
#endif


	if( sqrtm .lt. 8.94427190999916D0 )
     &    Warning('Extrapolating ZhSM('//Digit(h)//') in MHiggs')

        ZhSM(h) = exp(11.257012045186782D0 + 
     &     sqrtm*(-0.037183545587492574D0 + 
     &        sqrtm*(-0.08158474392414189D0 + 
     &           sqrtm*(0.004946796102005691D0 - 
     &              0.00010120919299465597D0*sqrtm))))

#ifdef DETAILED_DEBUG
	DPROD "ZhSM(h) =", ZhSM(h) ENDL
#endif

