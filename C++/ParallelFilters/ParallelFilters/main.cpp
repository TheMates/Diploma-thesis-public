// main.cpp : This file contains the 'main' function. Program execution begins and ends there.
//


#include <iostream>
#include "ParallelFilterDesign.h"
#include "rapidcsv.h"
#include <tuple>
#include "VectorExtensions.h"
#include <chrono>

#include "SplineNaK.h"

#define M_E 2.71828182845904523536


using namespace std;
using namespace ParallelFilters::Extensions;
using namespace VectorExtensions;

int Fs = 44100;		

using namespace Eigen;


int main()
{
	bool useNAK = true;
	string path = "";
	cout << "Fs is " << Fs << endl;

	while (path != "exit")
	{
		cout << "Name of var:" << endl;
		cin >> path;

		if(path == "Fs")
		{	
			std::string temp;
			cout << "Change Fs:" << endl;
			cin >> temp;
			Fs = stoi(temp);
			cout << endl << "Fs is " << Fs << endl;
			continue;
		}

		// Load data //
		//-----------//
		Vectord tf,W;
		Vectord2D data;
		try
		{
			loadRealFromCsvVec2D("input_data_csv/" + path + ".csv", data);
		}
		catch(exception &e)
		{
			cout << "Error loading file" << endl << endl;
			continue;
		}
		W = data[0];
		tf = data[1];


		// Get input data //
		//----------------//
		
		cplxVectord minPhMagn(W.size());
		ParallelFilters::minPhaseN(tf, W,minPhMagn, useNAK);

		// Parameters //
		//------------//

		cplxVectord X(W.size(), 1.0);
		
		
		/////////////////////////////
		// DUAL LOG   //

		auto nplog1 = 9, nplog2 = 11;		
		auto plog1 = ParallelFilters::Extensions::logspace(10, 500, nplog1);		
		auto plog2 = ParallelFilters::Extensions::logspace(780, 20000, nplog2);	
		plog1.insert(plog1.end(), plog2.begin(), plog2.end());		
		cplxVectord poles;
		ParallelFilters::freqpoles(plog1, poles,Fs);
		//
		/////////////////////////////
		auto coeffs = ParallelFilters::parFiltDesignFr(W, minPhMagn, poles);

		auto Bm = get<0>(coeffs);
		auto Am = get<1>(coeffs);
		auto FIR = get<2>(coeffs);

		writeToCsv(Bm, "coeffs_data/" + path + "_cpp_Bm_dl.csv");
		writeToCsv(Am, "coeffs_data/" + path + "_cpp_Am_dl.csv");
		writeToCsv(FIR, "coeffs_data/" + path + "_cpp_FIR_dl.csv");


		///////////////////////////////
		////   SINGLE WARP   //
		int pwp = 40;
		auto lambda = 0.92;

		poles = ParallelFilters::warpPolesFr(W, X, minPhMagn, pwp, lambda, 5);
		//
		////////////////////////////

		coeffs = ParallelFilters::parFiltDesignFr(W, minPhMagn, poles);

		Bm = get<0>(coeffs);
		Am = get<1>(coeffs);
		FIR = get<2>(coeffs);

		writeToCsv(Bm, "coeffs_data/" + path + "_cpp_Bm_w.csv");
		writeToCsv(Am, "coeffs_data/" + path + "_cpp_Am_w.csv");
		writeToCsv(FIR, "coeffs_data/" + path + "_cpp_FIR_w.csv");

		
		///////////////////////////
		//// DUAL WARP	//

		auto CrossFreq = 500.0;
		auto WindowCrossoverLength = 50;
		auto lambda1 = 0.986;
		auto lambda2 = 0.65;
		auto nPoles1 = 18;
		auto nPoles2 = 22;

		poles = ParallelFilters::dualWarpPolesFr(W, X, tf       , CrossFreq, WindowCrossoverLength, lambda1, lambda2, nPoles1, nPoles2, Fs, 5, useNAK);
		//////////////////////////////
		coeffs = ParallelFilters::parFiltDesignFr(W, minPhMagn, poles, 3); 

		Bm = get<0>(coeffs);
		Am = get<1>(coeffs);
		FIR = get<2>(coeffs);

		writeToCsv(Bm, "coeffs_data/" + path + "_cpp_Bm_dw.csv");
		writeToCsv(Am, "coeffs_data/" + path + "_cpp_Am_dw.csv");
		writeToCsv(FIR, "coeffs_data/" + path + "_cpp_FIR_dw.csv");

		cout << "Done" << endl<<endl;
	}
	std::cout << "Exit";

	return 0;
}

//benchmark testing
int main1()
{

 	string path = "1_1";

	// Load data //
	//-----------//
	Vectord tf, W;
	Vectord2D data;
	loadRealFromCsvVec2D("input_data_csv/" + path + ".csv", data);
	W = data[0];
	tf = data[1];

	// Get input data //
	//----------------//


	cplxVectord minPhMagn(W.size());
	ParallelFilters::minPhaseN(tf, W, minPhMagn,false);

	cplxVectord X(W.size(), 1.0);

	auto repetitions = 100;
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();		//tic

	for (auto rep = 0; rep < repetitions; ++rep)
	{

		///////////////////////////
		//// DUAL WARP	//

		auto CrossFreq = 500.0;
		auto WindowCrossoverLength = 50;
		auto lambda1 = 0.986;
		auto lambda2 = 0.65;
		auto nPoles1 = 18;
		auto nPoles2 = 22;

		auto poles = ParallelFilters::dualWarpPolesFr(W, X, tf,	CrossFreq, WindowCrossoverLength, lambda1, lambda2, nPoles1, nPoles2, Fs, 5);
		////
		//////////////////////////////

		// Par filter design //
		//-------------------//
		auto coeffs = ParallelFilters::parFiltDesignFr(W, minPhMagn, poles);

	}

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();	//toc

	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
	auto t = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
	cout << "Average time per computation: " << double(t) / repetitions << " [ms]" << endl;

	
	system("pause");

	std::cout << "Exit";
	return 0;

}