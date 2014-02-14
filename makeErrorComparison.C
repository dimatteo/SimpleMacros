// This macro makes a comparison of errors to understand 
// what is the correct treatment of syst uncertainties

using namespace std;

void makeErrorComparison(){

    //--> first let's initialize all the relevant numbers
    
    //we have six pt bins 
    const int nBins = 6;

    //yields are taken from the datacard used for the limit computation 
    //(as shown in pre-approval talk)
	float yields[nBins] = {65.9227,106.341,102.289,53.9474,8.5374,0.42635};
	totYield = 0;
	for (int iBin=0; iBin < nBins; iBin++)
		totYield += yields[iBin];
    //errors and k-factors are taken from AN-2012-439 v7 
	float errors[nBins] = {0.108,0.137,0.178,0.194,0.319,0.911};
	float kfactors[nBins] = {1.320,1.430,1.499,1.528,1.642,1.880};
    //total k-factor error is taken from AN, page 43 
	float totRelError = 0.133/1.424;
    //convert the absolute errors on k-factors into relative errors 
	for (int iBin=0; iBin < nBins; iBin++)
		errors[iBin] = errors[iBin]/kfactors[iBin];
    //rhos are taken from J. Neveu mail. 
    //here the assumption is that any bin pair will follow the pattern   
    //found by Jeremy when analyzing the first bin against the others   
	float rhos[nBins] = {1.,0.97,0.96,0.96,0.93,0.74};

    //--> second let's compute the errors for the three different scenarios

	//this is the error assuming no bin correlations
	float totErrorCorr = 0;
	//this is the error assuming bins are fully correlated
	float totErrorUncorr = 0;
	//this is the error assuming bins are smoothly correlated
	float totErrorSmoothCorr = 0;
	//loop on pt bins
	for (int iBin=0; iBin < nBins; iBin++){
		//full correlations: sum errors linearly
		totErrorCorr += yields[iBin]*errors[iBin];
		//no correlations: sum errors in quadrature
		totErrorUncorr += yields[iBin]*errors[iBin] * yields[iBin]*errors[iBin];
		//smooth correlations: sum errors in quadrature and add covarince terms
		totErrorSmoothCorr += yields[iBin]*errors[iBin] * yields[iBin]*errors[iBin];
		//add covariances for bin less than the last
	    if (iBin < nBins-1) 
	    	for (int jBin = iBin+1; jBin < nBins; jBin++)
	    		totErrorSmoothCorr += 2*rhos[jBin-iBin]*yields[iBin]*errors[iBin]*yields[jBin]*errors[jBin];
	}
	//correct the errors (after quadrature sum you have to take the sqrt)
	totErrorUncorr = sqrt(totErrorUncorr);
	totErrorSmoothCorr = sqrt(totErrorSmoothCorr);

	cout << "total number of expected events                  : " << totYield << endl;

	cout << "total relative error with correlation is         : " << totErrorCorr/totYield << endl;
	cout << "total relative error with correct correlation is : " << totErrorSmoothCorr/totYield << endl;
	cout << "total relative error w/o  correlation is         : " << totErrorUncorr/totYield << endl;
	cout << "AN total relative error is                       : " << totRelError << endl;

    //exit
    return;

}