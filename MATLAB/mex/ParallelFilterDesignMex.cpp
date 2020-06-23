/*
* Copyright (c) 2020 Matouš Vrbík
*
* URL	https://github.com/TheMates/Diploma-thesis-public
* 06/2020
* matousvrbik[at]gmail.com
*
* Distributed under CC BY-NC-SA 4.0 licence see LICENCE for details.
*/


#include "ParallelFilterDesign.h"
#include "VectorExtensions.h"
#include "mex.h"

// prhs contents:
// [0] H - magnitude of frequency response
// [1] W - angular frequencies
// [2] Fs - sample rate
// [3] lambda1 - first warping parameter
// [4] lambda2 - second warping parameter
// [5] nPoles1 - number of poles in first band
// [6] nPoles2 - number of poles in second band
// [7] crossFreq - crossing frequency of two bands
// [8] crossLength - number of frequency bins where bands cross
// [9] NFIR - number of FIR coefficients
// [10] useNAKspline - bool
//
//Returns plhs:
// [0] Bm - numerator coefficients
// [1] Am - denominator coefficients
// [2] FIR - FIR part coefficients
//
//Notes :
//number of iterations of lsidfr = 5, NFIR = 5
void mexFunction(int nlhs, mxArray* plhs[],
        int nrhs, const mxArray* prhs[]) {
    
    /* Check for proper number of arguments */
    if(nrhs != 11) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin",
                "ParalellFilterDesignMex requires eleven input arguments.");
    }
    if(nlhs != 3) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout",
                "ParalellFilterDesignMex requires three output argument.");
    }
    
    // H - vector of double
    if(!mxIsDouble(prhs[0]) ||                                    // not double
            mxIsComplex(prhs[0]) ||                                   // or complex
            mxIsScalar(prhs[0])) {                                  // or not scalar
        mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin",
                "1. argument H has to be double vector.");
    }
    if(mxGetN(prhs[0]) != 1 ){                                  // column vec
        mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin",
                "1. argument H has to be column vector.");
    }
    
    // W - vector of double
    if(!mxIsDouble(prhs[1]) ||                                    // not double
            mxIsComplex(prhs[1]) ||                                   // or complex
            mxIsScalar(prhs[1])) {                                  // or not scalar
        mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin",
                "2. argument W has to be double vector.");
    }
    if(mxGetN(prhs[1]) != 1 ){                                  // or not scalar
        mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin",
                "2. argument W has to be column vector.");
    }
    // Fs - int
    if(//mxIsDouble(prhs[3]) ||
            mxIsComplex(prhs[2]) ||
            !mxIsScalar(prhs[2])){
        mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin",
                "3. argument Fs has to be int scalar.");
    }
    // lambda1 - double
    if(!mxIsDouble(prhs[3]) ||
            mxIsComplex(prhs[3]) ||
            !mxIsScalar(prhs[3])){
        mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin",
                "4. argument lambda1 has to be double scalar.");
    }
    // lambda2 - double
    if(!mxIsDouble(prhs[4]) ||
            mxIsComplex(prhs[4]) ||
            !mxIsScalar(prhs[4])){
        mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin",
                "5. argument lambda2 has to be double scalar.");
    }
    // nPoles1 - int
    if(//mxIsDouble(prhs[6]) ||
            mxIsComplex(prhs[5]) ||
            !mxIsScalar(prhs[5])){
        mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin",
                "6. argument nPoles1 has to be int scalar.");
    }
    
    // nPoles2 - int
    if(//mxIsDouble(prhs[7]) ||
            mxIsComplex(prhs[6]) ||
            !mxIsScalar(prhs[6])){
        mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin",
                "7. argument nPoles2 has to be int scalar.");
    }
    
    // crossFreq - int
    if(//mxIsDouble(prhs[8]) ||
            mxIsComplex(prhs[7]) ||
            !mxIsScalar(prhs[7])){
        mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin",
                "8. argument crossFreq has to be int scalar.");
    }
    
    // crossLength - int
    if(//mxIsDouble(prhs[9]) ||
            mxIsComplex(prhs[8]) ||
            !mxIsScalar(prhs[8])){
        mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin",
                "9. argument crossFreq has to be int scalar.");
    }
    // NFIR - int
    if(//mxIsDouble(prhs[10]) ||
            mxIsComplex(prhs[9]) ||
            !mxIsScalar(prhs[9])){
        mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin",
                "10. argument NFIR has to be int scalar.");
    }
    
    // useNAKspline - bool
    if(     !mxIsLogical(prhs[10]) ||
            !mxIsScalar(prhs[10])){
        mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin",
                "11. argument useNAKspilne has to be bool scalar.");
    }
    
//     mexPrintf("All input variables are correct.\n");
    try
    {
        double Fs = *mxGetPr(prhs[2]);
        double lambda1 = *mxGetPr(prhs[3]);
        double lambda2 = *mxGetPr(prhs[4]);
        unsigned nPoles1 = static_cast<unsigned>( *mxGetPr(prhs[5]));
        unsigned nPoles2 = static_cast<unsigned>( *mxGetPr(prhs[6]));
        double crossFreqHz = *mxGetPr(prhs[7]);
        unsigned crossLength = static_cast<unsigned>( *mxGetPr(prhs[8]));
        unsigned NFIR = static_cast<unsigned>( *mxGetPr(prhs[9]));
        bool useNAKspline = ( *mxGetLogicals (prhs[10]));
        
        Vectord H(mxGetPr(prhs[0]),mxGetPr(prhs[0])+mxGetM(prhs[0]) );
        Vectord W(mxGetPr(prhs[1]),mxGetPr(prhs[1])+mxGetM(prhs[1]) );
        cplxVectord minPhMagn(W.size());
        
        ParallelFilters::minPhaseN(H, W,minPhMagn,useNAKspline);
        
        cplxVectord X(W.size(),1.0);
        
        auto poles = ParallelFilters::dualWarpPolesFr(W, X, H,
                crossFreqHz, crossLength, lambda1, lambda2, nPoles1, nPoles2, Fs, 5, useNAKspline);
        
 
        auto coeffs = ParallelFilters::parFiltDesignFr(W, minPhMagn, poles,NFIR);
        
        auto Bm = std::get<0>(coeffs);
        auto Am = std::get<1>(coeffs);
        auto FIR = std::get<2>(coeffs);
        
        plhs[0] = mxCreateDoubleMatrix(2,Bm.size(), mxREAL);        //Bm
        plhs[1] = mxCreateDoubleMatrix(3,Am.size(), mxREAL);        //Am
        plhs[2] = mxCreateDoubleMatrix(FIR.size(),1, mxREAL);       //FIR
        
        auto outBm = mxGetPr(plhs[0]);
        auto outAm = mxGetPr(plhs[1]);
        auto outFIR = mxGetPr(plhs[2]);
        
        for (auto filter = 0; filter<Bm.size();++filter)
        {
            *(outBm + filter*2) = Bm[filter][0];
            *(outBm + filter*2 + 1) = Bm[filter][1];
            *(outAm + filter*3) = Am[filter][0];
            *(outAm + filter*3 + 1) = Am[filter][1];
            *(outAm + filter*3 + 2) = Am[filter][2];
        }
        for(auto fir = 0; fir< NFIR; ++fir)
            *(outFIR + fir) = FIR[fir];
    }
    catch(std::exception &e)
    {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin",
                e.what());
    }
    
}