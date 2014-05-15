#include <iostream>
#include <vector>
#include "TRandom3.h"
#include "TMath.h"
//#include <stdio.h>

void smearErrorsBR () {

  const bool addPUandTHUlin = true;
  const bool deriveTHUsf = false; 
  const int  THUshape    = 1; // 0: gaus, 1: box
  const double THUBoxErrorScaleFactor = 1.;
  const bool THU100percCorr = false;
  const bool THU100percCorrToPU = false;
  const int nPU = 4;
  const int nBR  = 9;
  const int nToys = 1000000;
  vector< vector<double> > relUncPU;
  vector<double> relUncTHU;
  vector<double> referenceBR;
  vector<double> centralValueGamma; 
  vector< vector<double> > toyBR;
  vector< vector<double> > toyGammas;
  vector< vector<double> > covarianceMatrixBR; 
  vector< vector<double> > covarianceMatrixGamma;
  const double mtsf = 1.0;
  const double mqsf = 1.0;
  const double assf = 1.0;
  const double thsf = 1.0;   

  TRandom3 random;

  cout << "#----------------------------------------------------#" << endl;
  cout << "#                                                    #" << endl;
  cout << "#      Covariance matrix calculator for the          #" << endl;
  cout << "#    Higgs boson branching ratio uncertainties       #" << endl;
  cout << "#                                                    #" << endl;
  cout << "#          This ROOT macro is part of the            #" << endl;
  cout << "#              HiggsSignals package.                 #" << endl;
  cout << "#                                                    #" << endl;
  cout << "#      (latest change on 25/09/2013, PB & TS)        #" << endl;
  cout << "#----------------------------------------------------#" << endl;
  cout << endl;
  cout << "Filling input relative parameteric and theoretical uncertainties..." << endl;
  
  // fill inputs
  vector<double> relUncPUbb;
  relUncPUbb.push_back(assf*(-0.023));// bb alphas
  relUncPUbb.push_back(mqsf*0.033);// bb mb
  relUncPUbb.push_back(mqsf*0.0  );// bb mc
  relUncPUbb.push_back(mtsf*0.0  );// bb mt

  vector<double> relUncPUtautau;
  relUncPUtautau.push_back(assf*0.0  );// tautau alphas
  relUncPUtautau.push_back(mqsf*0.0  );// tautau mb
  relUncPUtautau.push_back(mqsf*0.0  );// tautau mc
  relUncPUtautau.push_back(mtsf*0.001);// tautau mt

  vector<double> relUncPUmumu;
  relUncPUmumu.push_back(assf*0.0  );// mumu alphas
  relUncPUmumu.push_back(mqsf*0.0  );// mumu mb
  relUncPUmumu.push_back(mqsf*(-0.001));// mumu mc
  relUncPUmumu.push_back(mtsf*0.001);// mumu mt

  vector<double> relUncPUcc;
  relUncPUcc.push_back(assf*(-0.071));// cc alphas
  relUncPUcc.push_back(mqsf*(-0.001));// cc mb
  relUncPUcc.push_back(mqsf*0.062);// cc mc
  relUncPUcc.push_back(mtsf*0.001);// cc mt

  vector<double> relUncPUgg;
  relUncPUgg.push_back(assf*0.042);// gg alphas
  relUncPUgg.push_back(mqsf*(-0.001));// gg mb
  relUncPUgg.push_back(mqsf*0.000);// gg mc
  relUncPUgg.push_back(mtsf*(-0.002));// gg mt

  vector<double> relUncPUgaga;
  relUncPUgaga.push_back(assf*0.0  );// gaga alphas
  relUncPUgaga.push_back(mqsf*0.0  );// gaga mb
  relUncPUgaga.push_back(mqsf*0.0  );// gaga mc
  relUncPUgaga.push_back(mtsf*0.0  );// gaga mt

  vector<double> relUncPUZga;
  relUncPUZga.push_back(assf*0.0  );// Zga alphas
  relUncPUZga.push_back(mqsf*0.0  );// Zga mb
  relUncPUZga.push_back(mqsf*0.001);// Zga mc
  relUncPUZga.push_back(mtsf*0.001);// Zga mt

  vector<double> relUncPUWW;
  relUncPUWW.push_back(assf*0.0  );// WW alphas
  relUncPUWW.push_back(mqsf*0.0  );// WW mb
  relUncPUWW.push_back(mqsf*0.0  );// WW mc
  relUncPUWW.push_back(mtsf*0.0  );// WW mt

  vector<double> relUncPUZZ;
  relUncPUZZ.push_back(assf*0.0  );// ZZ alphas
  relUncPUZZ.push_back(mqsf*0.0  );// ZZ mb
  relUncPUZZ.push_back(mqsf*0.0  );// ZZ mc
  relUncPUZZ.push_back(mtsf*0.0  );// ZZ mt

// n.b.: Need the following ordering for HiggsSignals!

  relUncPU.push_back(relUncPUgaga);
  relUncPU.push_back(relUncPUWW);
  relUncPU.push_back(relUncPUZZ);
  relUncPU.push_back(relUncPUtautau);  
  relUncPU.push_back(relUncPUbb);
  relUncPU.push_back(relUncPUZga);
  relUncPU.push_back(relUncPUcc);
  relUncPU.push_back(relUncPUmumu);  
  relUncPU.push_back(relUncPUgg);

  relUncTHU.push_back(thsf*0.010);// gaga
  relUncTHU.push_back(0.005);// WW
  relUncTHU.push_back(0.005);// ZZ     
  relUncTHU.push_back(thsf*0.020);// tautau
  relUncTHU.push_back(thsf*0.020);// bb
  relUncTHU.push_back(thsf*0.050);// Zga
  relUncTHU.push_back(thsf*0.020);// cc
  relUncTHU.push_back(thsf*0.020);// mumu
  relUncTHU.push_back(thsf*0.030);// gg

  cout << "Filling reference branching ratios..." << endl;
  
// Please enter your reference BR values here:

  referenceBR.push_back(0.00228);// gaga
  referenceBR.push_back(0.231);// WW
  referenceBR.push_back(0.0289);// ZZ  
  referenceBR.push_back(0.0616);// tautau    
  referenceBR.push_back(0.561);// bb
  referenceBR.push_back(0.00162);// Zga  
  referenceBR.push_back(0.0283);// cc
  referenceBR.push_back(0.000214);// mumu  
  referenceBR.push_back(0.0848);// gg

  cout << "Filling input central values for partial widths..." << endl;
  centralValueGamma.push_back(0.00000959 );//gaga   
  centralValueGamma.push_back(0.000973   );//WW	    
  centralValueGamma.push_back(0.000122   );//ZZ     
  centralValueGamma.push_back(0.000259   );//tautau 
  centralValueGamma.push_back(0.00236    );//bb
  centralValueGamma.push_back(0.00000684 );//Zga    
  centralValueGamma.push_back(0.000119   );//cc	    
  centralValueGamma.push_back(0.000000899);//mumu   
  centralValueGamma.push_back(0.000357   );//gg	    
  
  cout << "Start generating toys..." << endl;

  // make toys
  for (int iToy = 0; iToy < nToys; iToy++) {
    //cout << "toy " << iToy << "\n";
    // first calculate all gammas
    vector<double> thisPUsmearing;
    double THUscalefactor = 0.;
    for (int iPU = 0; iPU < nPU; iPU++) {
      thisPUsmearing.push_back(random.Gaus(0,1));
      THUscalefactor += thisPUsmearing[iPU]/2.;
    }
    double THURandomNumber = 0;
    if (THUshape == 0) {
      THURandomNumber = random.Gaus(0,1);
    } else if (THUshape == 1) {
      THURandomNumber = random.Uniform(-THUBoxErrorScaleFactor,THUBoxErrorScaleFactor);
    } else {
      cout << "nothing else implemented" << endl;
      return;
    }
    vector<double> theseGammas;
    for (int iBR = 0; iBR < nBR; iBR++) {
      theseGammas.push_back(centralValueGamma[iBR]);
      for (int iPU = 0; iPU < nPU; iPU++) {
	theseGammas[iBR] += 
	  centralValueGamma[iBR]*thisPUsmearing[iPU]*relUncPU[iBR][iPU];
      }
      if (!addPUandTHUlin) {
	double x;
	if (THUshape == 0) {
	  x = random.Gaus(0,1);
	} else if (THUshape == 1) {
	  x = random.Uniform(-THUBoxErrorScaleFactor,THUBoxErrorScaleFactor);
	} else {
	  cout << "nothing else implemented" << endl;
	  return;
	}
	theseGammas[iBR] += 
	  centralValueGamma[iBR]*x*relUncTHU[iBR];
      } else {
	double x;
	if(deriveTHUsf){
	  x = TMath::Abs(THUscalefactor);
	  if (theseGammas[iBR]>centralValueGamma[iBR]) {
	    theseGammas[iBR] += 
	      centralValueGamma[iBR]*x*relUncTHU[iBR];
	  } else {
	    theseGammas[iBR] += 
	      -centralValueGamma[iBR]*x*relUncTHU[iBR];
	  }
	} else {
	  if (THU100percCorr) {
	    x = THURandomNumber;
	    theseGammas[iBR] += 
	      centralValueGamma[iBR]*x*relUncTHU[iBR];
	  } else {
	    if (THUshape == 0) {
	      x = random.Gaus(0,1);
	    } else if (THUshape == 1) {
	      x = random.Uniform(-THUBoxErrorScaleFactor,THUBoxErrorScaleFactor);
	    } else {
	      cout << "nothing else implemented" << endl;
	      return;
	    }
	    if (THU100percCorrToPU) {
	      x = TMath::Abs(x);
	      if (theseGammas[iBR]>centralValueGamma[iBR]) {
		theseGammas[iBR] += 
		  centralValueGamma[iBR]*x*relUncTHU[iBR];
	      } else {
		theseGammas[iBR] += 
		  -centralValueGamma[iBR]*x*relUncTHU[iBR];
	      }	
	    } else {
	      theseGammas[iBR] += 
		centralValueGamma[iBR]*x*relUncTHU[iBR];
	    }
	  }
	}
      }
    } 
    // then calculate GammaTot
    double thisGammaTot = 0;
    for (int iBR = 0; iBR < nBR; iBR++) {
      thisGammaTot += theseGammas[iBR];
    }    
    // then calculate BR's
    vector<double> theseBRs;
    for (int iBR = 0; iBR < nBR; iBR++) {
      theseBRs.push_back(theseGammas[iBR]/thisGammaTot);
    }    
    // then store result
    toyBR.push_back(theseBRs);
    toyGammas.push_back(theseGammas);
  }
  cout << "Calculated the toys, now calculate the covariance matrix..." << endl;
  
// calculate covariances for the partial widths
  for (int iBR = 0; iBR < nBR; iBR++) {
    // calculate E[iBR]
    double meanGammaI = 0;
    for (int iToy = 0; iToy < nToys; iToy++) {
      meanGammaI += toyGammas[iToy][iBR];
    }
    meanGammaI = meanGammaI/(double)nToys;
    vector<double> thisCovMatrixRow;
    for (int jBR = 0; jBR < nBR; jBR++) {
      // calculate E[jBR] and E[iBR*jBR]
      double meanGammaJ = 0;
      double meanGammaIJ = 0;
      for (int iToy = 0; iToy < nToys; iToy++) {
	meanGammaJ  += toyGammas[iToy][jBR];
	meanGammaIJ += toyGammas[iToy][jBR]*toyGammas[iToy][iBR];
      }
      meanGammaJ = meanGammaJ/(double)nToys;
      meanGammaIJ = meanGammaIJ/(double)nToys;
      thisCovMatrixRow.push_back(meanGammaIJ - meanGammaI*meanGammaJ);
    }      
    covarianceMatrixGamma.push_back(thisCovMatrixRow);
  }      
  
  // calculate covariances for the BRs
  for (int iBR = 0; iBR < nBR; iBR++) {
    // calculate E[iBR]
    double meanI = 0;
    for (int iToy = 0; iToy < nToys; iToy++) {
      meanI += toyBR[iToy][iBR];
    }
    meanI = meanI/(double)nToys;
    vector<double> thisCovMatrixRow;
    for (int jBR = 0; jBR < nBR; jBR++) {
      // calculate E[jBR] and E[iBR*jBR]
      double meanJ = 0;
      double meanIJ = 0;
      for (int iToy = 0; iToy < nToys; iToy++) {
	meanJ  += toyBR[iToy][jBR];
	meanIJ += toyBR[iToy][jBR]*toyBR[iToy][iBR];
      }
      meanJ = meanJ/(double)nToys;
      meanIJ = meanIJ/(double)nToys;
      thisCovMatrixRow.push_back(meanIJ - meanI*meanJ);
    }      
    covarianceMatrixBR.push_back(thisCovMatrixRow);
  }    
//   cout << "#---------------------------------------------------------------"<< endl;
//   cout << "Covariance matrix for the partial widths:\n" << endl;
//  
//   for (int iBR = 0; iBR < nBR; iBR++) {
//     for (int jBR = 0; jBR < nBR; jBR++) {
//       cout << " " << covarianceMatrixGamma[iBR][jBR];
//     }    
//     cout << endl;
//   }
//   cout << "\nCorrelation matrix for the partial widths:\n" << endl;
// 
//   for (int iBR = 0; iBR < nBR; iBR++) {
//     for (int jBR = 0; jBR < nBR; jBR++) {
//       cout << " " << covarianceMatrixGamma[iBR][jBR]/TMath::Sqrt(covarianceMatrixGamma[iBR][iBR]*covarianceMatrixGamma[jBR][jBR]);
//     }    
//     cout << endl;
//   }
// 
//   cout << "#---------------------------------------------------------------"<< endl;
// 
//   // print covariances
//   cout << "\nCovariance matrix for the BR with absolute errors:\n" << endl;
//  
//   for (int iBR = 0; iBR < nBR; iBR++) {
//     for (int jBR = 0; jBR < nBR; jBR++) {
//       cout << " " << covarianceMatrixBR[iBR][jBR];
//     }    
//     cout << endl;
//   }
  
  cout << "\nCovariance matrix for the BR with relative errors:" << endl;
  cout << "(to be used in HiggsSignals under the name BRcov.in or BRcovSM.in)\n" << endl;
  
  
  for (int iBR = 0; iBR < nBR; iBR++) {
    for (int jBR = 0; jBR < nBR; jBR++) {
      cout << " " << covarianceMatrixBR[iBR][jBR]/(referenceBR[iBR]*referenceBR[jBR]);
    }    
    cout << endl;
  }  

  cout << "\nCorrelation matrix (reference BR's cancel out!):\n" << endl;

  for (int iBR = 0; iBR < nBR; iBR++) {
    for (int jBR = 0; jBR < nBR; jBR++) {
      cout << " " << covarianceMatrixBR[iBR][jBR]/TMath::Sqrt(covarianceMatrixBR[iBR][iBR]*covarianceMatrixBR[jBR][jBR]);
    }    
    cout << endl;
  }
  
  
  
  cout << "#---------------------------------------------------------------"<< endl;
  cout << "FOR BETTER READING" << endl;
  cout << "#---------------------------------------------------------------"<< endl;
  const char * format = "%.1f\t";
  cout << "Correlation matrix (in %) for the partial widths (gaga, WW, ZZ, tautau, bb, Zga, cc, mumu, gg):" << endl;
  cout << endl;
  
  for (int iBR = 0; iBR < nBR; iBR++) {
    for (int jBR = 0; jBR < nBR; jBR++) {
      printf (format,100.*covarianceMatrixGamma[iBR][jBR]/TMath::Sqrt(covarianceMatrixGamma[iBR][iBR]*covarianceMatrixGamma[jBR][jBR]));
    }    
    cout << endl;
  }
  cout << "#---------------------------------------------------------------"<< endl;
  cout << "Correlation matrix (in %) for BRs (gaga, WW, ZZ, tautau, bb, Zga, cc, mumu, gg):" << endl;
  cout << endl;
  
  for (int iBR = 0; iBR < nBR; iBR++) {
    for (int jBR = 0; jBR < nBR; jBR++) {
      printf (format,100.*covarianceMatrixBR[iBR][jBR]/TMath::Sqrt(covarianceMatrixBR[iBR][iBR]*covarianceMatrixBR[jBR][jBR]));
    }    
    cout << endl;
  }
cout << "#---------------------------------------------------------------"<< endl;

 // gaga
 // WW
 // ZZ     
 // tautau
 // bb
 // Zga
 // cc
 // mumu
 // gg

 cout << "comparison quantities:" << endl;
 cout << "relative error on BR(h->gaga) =   " <<
   TMath::Sqrt(covarianceMatrixBR[0][0]/(referenceBR[0]*referenceBR[0])) << endl;
 cout << "relative error on BR(h->WW) =     " <<
   TMath::Sqrt(covarianceMatrixBR[1][1]/(referenceBR[1]*referenceBR[1])) << endl;
 cout << "relative error on BR(h->ZZ) =     " <<
   TMath::Sqrt(covarianceMatrixBR[2][2]/(referenceBR[2]*referenceBR[2])) << endl;
 cout << "relative error on BR(h->tautau) = " <<
   TMath::Sqrt(covarianceMatrixBR[3][3]/(referenceBR[3]*referenceBR[3])) << endl;
 cout << "relative error on BR(h->bb) =     " <<
   TMath::Sqrt(covarianceMatrixBR[4][4]/(referenceBR[4]*referenceBR[4])) << endl;
 cout << "relative error on BR(h->Zga) =    " <<
   TMath::Sqrt(covarianceMatrixBR[5][5]/(referenceBR[5]*referenceBR[5])) << endl;
 cout << "relative error on BR(h->cc) =     " <<
   TMath::Sqrt(covarianceMatrixBR[6][6]/(referenceBR[6]*referenceBR[6])) << endl;
 cout << "relative error on BR(h->mumu) =   " <<
   TMath::Sqrt(covarianceMatrixBR[7][7]/(referenceBR[7]*referenceBR[7])) << endl;
 cout << "relative error on BR(h->gg) =     " <<
   TMath::Sqrt(covarianceMatrixBR[8][8]/(referenceBR[8]*referenceBR[8])) << endl;

  
  return;

}
