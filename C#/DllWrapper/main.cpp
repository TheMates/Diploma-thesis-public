#pragma once
#include "C:\Users\Matouš\source\repos\Diploma-thesis-GitHub\CD_EXPORT\C++\ParallelFilters\ParalelFilters\ParallelFilters"

using namespace ParallelFilters;
using namespace std;

#define EXPORT_TYPE extern "C" __declspec(dllexport) 


/// <summary>
/// Loads reponse from csv file. First column is angular frequency vector, second column is the response magnitude (not in dB).
/// w and H must be preallocated.
/// Returns false if loaded data were longer than preallocated, than a new allocation must be done and called again.
/// </summary>
/// <param name="w">Angular frequency vector.</param>
/// <param name="H">Reponse magnitude vector.</param>
/// <param name="size">Preallocated size of the outputs. If loaded files have greater size, function returns false.</param>
/// <param name="path">Path to file.</param>
/// <returns>True if loading was successfull.</returns>
EXPORT_TYPE bool loadResponse(double w[], double H[], int * size, const char * path)
{
	// tady musím nìjakým zpùsobem pøijít k té odezvì - naèíst ze souboru
	Vectord2D data;
	try
	{
		Extensions::loadRealFromCsvVec2D(path, data);
	}
	catch(exception &e)
	{
		*size = -1;
		return false;
	}
	if (data.size() != 2)			// wrong file format - not 2 columns
	{
		*size = -1;
		return false;
	}
	if (data[0].size() > * size)	//loaded file size was greater than allocated
		return false;
	
	std::copy(data[0].begin(), data[0].end(), w);
	std::copy(data[1].begin(), data[1].end(), H);
	*size = data[0].size();
	return true;
}

/// <summary>
/// Computes filter response with dual warping algorithm. 
/// </summary>
/// <param name="response">Preallocated output array of length w.</param>
/// <param name="w">Angular frequencies.</param>
/// <param name="target">Target frequency response magnitude.</param>
/// <param name="size">Length of input arrays.</param>
/// <param name="nPoles1">Number of poles low.</param>
/// <param name="nPoles2">Number of poles high.</param>
/// <param name="crossFreq">Cross frequency in Hz.</param>
/// <param name="crossLength">Cross length in samples.</param>
/// <param name="lambda1">Lambda warpint parameter low.</param>
/// <param name="lambda2">Lambda warpint parameter high.</param>
/// <param name="sampleRate">Sample rate.</param>
/// <param name="NFIR">Number of FIR coefficients.</param>
/// <param name="useNAK">Use not-a-knot spline</param>
EXPORT_TYPE void computeResponse(double response[], double w[], double target[], int* size, int nPoles1, int nPoles2, double crossFreq, int crossLength, double lambda1, double lambda2, double sampleRate,int NFIR, bool useNAK )
{
	//spoèítat koeficienty, a rovnou aji magn odezvu a hodit do response
	
	Vectord W(w,w+*size);
	Vectord tf(target, target + *size);
	cplxVectord minPhMagn(W.size());
	ParallelFilters::minPhaseN(tf, W, minPhMagn, useNAK);
	cplxVectord X(W.size(), 1.0);
	auto poles = ParallelFilters::dualWarpPolesFr(W, X, tf, crossFreq, crossLength, lambda1, lambda2, nPoles1, nPoles2, sampleRate, 5, useNAK);
	auto coeffs = ParallelFilters::parFiltDesignFr(W, minPhMagn, poles,NFIR);

	auto Bm = get<0>(coeffs);
	auto Am = get<1>(coeffs);
	auto FIR = get<2>(coeffs);

	Vectord magn(*size);
	parFiltFreqResp(Bm, Am, FIR, W, magn);

	std::copy(magn.begin(), magn.end(), response);

}

/// <summary>
/// Computes parameters of parallel filter simulating the given response and saves it to .csv file.
/// </summary>
/// <param name="w">Angular frequencies.</param>
/// <param name="target">Target frequency response magnitude.</param>
/// <param name="size">Length of input arrays.</param>
/// <param name="nPoles1">Number of poles low.</param>
/// <param name="nPoles2">Number of poles high.</param>
/// <param name="crossFreq">Cross frequency in Hz.</param>
/// <param name="crossLength">Cross length in samples.</param>
/// <param name="lambda1">Lambda warpint parameter low.</param>
/// <param name="lambda2">Lambda warpint parameter high.</param>
/// <param name="sampleRate">Sample rate.</param>
/// <param name="NFIR">Number of FIR coefficients.</param>
/// <param name="useNAK">Use not-a-knot spline</param>
/// <param name="fname">Name of output file.</param>
EXPORT_TYPE void exportToCsv(double w[], double target[], int* size, int nPoles1, int nPoles2, double crossFreq, int crossLength, double lambda1, double lambda2, double sampleRate, int NFIR, bool useNAK, const char* fname)
{
	Vectord W(w, w + *size);
	Vectord tf(target, target + *size);
	cplxVectord minPhMagn(W.size());
	ParallelFilters::minPhaseN(tf, W, minPhMagn, useNAK);
	cplxVectord X(W.size(), 1.0);
	auto poles = ParallelFilters::dualWarpPolesFr(W, X, tf, crossFreq, crossLength, lambda1, lambda2, nPoles1, nPoles2, sampleRate, 5, useNAK);
	auto coeffs = ParallelFilters::parFiltDesignFr(W, minPhMagn, poles, NFIR);

	auto Bm = get<0>(coeffs);
	auto Am = get<1>(coeffs);
	auto FIR = get<2>(coeffs);
	string fileName(fname);
	Extensions::writeToCsv(Bm, fileName + "_Bm.csv");
	Extensions::writeToCsv(Am, fileName + "_Am.csv");
	Extensions::writeToCsv(FIR, fileName + "_FIR.csv");

}

/// <summary>
/// Private method, that returns the Bm, Am, FIR coefficients as 1D vectors.
/// </summary>
/// <param name="Am">Denominator coefficients.</param>
/// <param name="Bm">Nuemrator coefficients.</param>
/// <param name="FIR">FIR coefficients.</param>
/// <param name="w">Angular frequencies.</param>
/// <param name="target">Target frequency response magnitude.</param>
/// <param name="size">Length of input arrays.</param>
/// <param name="nPoles1">Number of poles low.</param>
/// <param name="nPoles2">Number of poles high.</param>
/// <param name="crossFreq">Cross frequency in Hz.</param>
/// <param name="crossLength">Cross length in samples.</param>
/// <param name="lambda1">Lambda warpint parameter low.</param>
/// <param name="lambda2">Lambda warpint parameter high.</param>
/// <param name="sampleRate">Sample rate.</param>
/// <param name="NFIR">Number of FIR coefficients.</param>
/// <param name="useNAK">Use not-a-knot spline</param>
EXPORT_TYPE void computeCoeffs(double Am[], double Bm[] ,double FIR[], double w[], double target[], int size, int nPoles1, int nPoles2, double crossFreq, int crossLength, double lambda1, double lambda2, double sampleRate, int NFIR, bool useNAK)
{
	Vectord W(w, w + size);
	Vectord tf(target, target + size);
	cplxVectord minPhMagn(W.size());
	ParallelFilters::minPhaseN(tf, W, minPhMagn, useNAK);
	cplxVectord X(W.size(), 1.0);
	auto poles = ParallelFilters::dualWarpPolesFr(W, X, tf, crossFreq, crossLength, lambda1, lambda2, nPoles1, nPoles2, sampleRate, 5, useNAK);
	auto coeffs = ParallelFilters::parFiltDesignFr(W, minPhMagn, poles, NFIR);

	Vectord2D Bmret = get<0>(coeffs);
	Vectord2D Amret = get<1>(coeffs);
	Vectord FIRret = get<2>(coeffs);
	
	Vectord BmCat(Bmret.size() * 2), AmCat(Amret.size() * 2);
	for(auto i =0;i<Bmret.size();++i)
	{
		Bm[i] = Bmret[i][0];
		Bm[i + Bmret.size()] = Bmret[i][1];
		Am[i] = Amret[i][1];
		Am[i + Amret.size()] = Amret[i][2];
	}
	std::copy(FIRret.begin(), FIRret.end(), FIR);
}
