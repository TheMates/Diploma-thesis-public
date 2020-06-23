/*
* Copyright (c) 2020 Matouš Vrbík
*
* URL	https://github.com/TheMates/Diploma-thesis-public
* 06/2020
* matousvrbik[at]gmail.com
*
* Distributed under MIT licence see LICENCE for details.
*/

#pragma once
#include "ParallelFilterDesign.h"
#define _USE_MATH_DEFINES
#include <math.h>

#include <functional>
#include <algorithm>
#include <memory>
#include "Eigen/Eigen"
#include "rapidcsv.h"
#include "Extensions.h"
#include "VectorExtensions.h"


using namespace std;
using namespace ParallelFilters::Extensions;
using namespace VectorExtensions;

namespace ParallelFilters
{
	

	void freqpoles(const Vectord& fr, cplxVectord &output, const double Fs,const double Q)
	{
		auto pnum = fr.size();
		vector<double> wp = vector<double>(fr);

		transform(wp.begin(), wp.end(), wp.begin(),
			bind(multiplies<float>(), placeholders::_1, 2 * M_PI / Fs));	//changes the fr to normalized frequencies  

		std::unique_ptr<vector<double>> dwp = std::make_unique<vector<double>>(pnum);	//vector dwp, of size of number of frequencies
		if (Q <= 0.0)	//if Q not specified
		{
			for (auto i = 1; i < pnum - 1; i++)
				dwp->operator[](i) = (wp[i+1] - wp[i - 1]) / 2.0;
			dwp->operator[](0) = wp[1] - wp[0];
			dwp->operator[](pnum - 1) = wp[pnum - 1] - wp[pnum - 2];
		}
		else
		{
			for (auto i = 0; i < pnum; i++)
				dwp->operator[](i) = wp[i] / Q;
		}
		output.resize(pnum * 2);
		for (auto i = 0; i < pnum; i++)
			output[i] = exp(-(dwp->operator[](i) / 2.0)) * exp(complex<double>(0.0, wp[i]));
		for (auto i = pnum; i < pnum * 2; i++)
			output[i] = conj(output[i - pnum]);
	}

	void solve(Eigen::MatrixXd& A, Eigen::MatrixXd& b, Eigen::MatrixXd &x)
	{
		x = A.ldlt().solve(b);
	}


	tuple<vector<Vectord>, vector<Vectord>,Vectord> parFiltDesignFr(
		vector<double> W, cplxVectord H, cplxVectord poles, unsigned nFIR )
	{
		if (W.size() != H.size()) throw length_error("W and H length mismatch!");
		 
		//flip poles inside unit circle
		for (auto& pole : poles)
		{
			if (abs(pole) > 1)
				pole = 1.0 / conj(pole);
			if (abs(pole) > 0.9995 && imag(pole) == 0)
				pole = 0.995;
		}

		cplxpair(poles);

		int polesNumber = poles.size();
		int evenPolesNumber = int(2*floor( poles.size()/2));
		auto ODD = false;
		if (polesNumber > evenPolesNumber)
			ODD = true;
		
		Eigen::VectorXcd Z(W.size());

		for (int i = 0; i < W.size(); i++)
		{
			Z(i) = exp(complex<double>(0.0, -W[i]));
		}

		Eigen::VectorXcd z2 = Z.cwiseProduct(Z);

		// now making second-order sections //
		//----------------------------------//

		Eigen::MatrixXcd M(W.size(), polesNumber + nFIR);

		for (auto k = 0; k < polesNumber - 1; k += 2)
		{
			cplxVectord tempRoots(poles.begin()+k, poles.begin()+k+2);
			auto A = poly(tempRoots);		
			M.col(k) = (A[0] + (A[1] * Z + A[2]*z2).array()).cwiseInverse();
			M.col(k + 1) = Z.cwiseProduct(M.col(k));
		}

		if (ODD)
		{
			cplxVectord temp(poles.end()-1,poles.end());
			auto A = poly(temp);
			M.col(polesNumber - 1) = (A[0] + (A[1] * Z).array()).cwiseInverse();
		}

		for(unsigned i = 0;i<nFIR;++i)
		{
			M.col(polesNumber + i) = Z.array(). pow(i);
		}
		
		// least squares solution - using Eigen //
		//--------------------------------------//

		Eigen::Map<Eigen::VectorXcd> H1(H.data(), H.size());

		Eigen::MatrixXd A = (M.transpose().conjugate() * M).real();
		Eigen::MatrixXd b = (M.transpose().conjugate() * H1).real();
		
		Eigen::MatrixXd par;
		solve(A, b, par);

		// constructing Bm and Am coefficients //
		//-------------------------------------//

		Vectord2D Am( int(ceil(double( polesNumber)/2)), Vectord(3)), Bm(int(ceil(double(polesNumber) / 2)), Vectord(2));
		
		for(auto i = 0;i<evenPolesNumber/2;++i)
		{
			auto tempPoles = cplxVectord(&poles[2 * i], &poles[2 * i + 1]+1);
			auto tempA = poly(tempPoles);
			Am[i] = real(tempA);
			Bm[i] = Vectord(&par(2 * i), &par(2 * i + 1) + 1);
		}

		if (ODD)
		{
			auto tempPoles = cplxVectord(poles.end()-1,poles.end());
			auto tempA = poly(tempPoles);
			auto rTempA = real(tempA);
			rTempA.push_back(0);
			Am[evenPolesNumber/2] = rTempA;

			Vectord tempPar{ par(polesNumber - 1) };
			tempPar.push_back(0);
			Bm[evenPolesNumber / 2] = tempPar;
		}

		Vectord FIRcoeffs(nFIR);
		if(nFIR>0)
		{
			for (auto i = 0; i < nFIR; ++i)
				FIRcoeffs[i] = par(polesNumber + i);
		}

		return tuple<vector<Vectord>, vector<Vectord>,Vectord>(Bm,Am,FIRcoeffs);
	}


	std::tuple < Vectord, Vectord> lsidFr(Vectord& W, cplxVectord& X, cplxVectord& Y, unsigned nB, unsigned nA,
		unsigned iter)
	{
		auto PARL = nA + nB + 1;
		auto L = Y.size();

		cplxVectord Z(L);
		for (auto i = 0; i < W.size(); i++)
		{
			Z[i] = exp(complex<double>(0.0, -W[i]));
		}

		Eigen::MatrixXcd M(L, PARL);
		for (auto k = 0; k <= nB; ++k)
		{
			for (auto i = 0; i < L; ++i)
			{
				M(i,k) = X[i] * pow(Z[i], k);
			}
		}
		for (auto k = 1; k <= nA; ++k)
		{
			for (auto i = 0; i < L; ++i)
			{
				M(i,k + nB) = -Y[i] * pow(Z[i], k);
			}
		}

		//transform Y into Eigen type
		Eigen::Map<Eigen::VectorXcd> Y1(Y.data(), Y.size());
			   
		Eigen::MatrixXd Am = (M.transpose().conjugate() * M).real();
		Eigen::MatrixXd b = (M.transpose().conjugate() * Y1).real();

		//LLT nebo LDLT
		Eigen::MatrixXd par;
		solve(Am, b, par);

		Vectord B(nB + 1);
		Vectord A(nA + 1);

		for (auto i = 0; i < nB; ++i)
		{
			B[i] = par(i);
		}
		for (auto i = 1; i <= nA; ++i)
		{
			A[i] = par(nB + i);
		}
		A[0] = 1.0;


		// Iterations //
		//------------//
		auto MW = M;
		Eigen::VectorXcd YW(Y1.size());	

		Eigen::VectorXcd AWt(Z.size());
		Eigen::Map<Eigen::VectorXcd> Z1(Z.data(), Z.size());

		for (auto k = 1; k <= iter; ++k)
		{
			AWt.setZero();
			auto tempA(A);
			tempA.insert(tempA.begin(), 1.0);
			polyval(tempA, Z1,AWt);		
			AWt = AWt.array().abs().real().inverse();

			for (auto col = 0; col < MW.cols(); ++col)
			{
				MW.col(col) = M.col(col).cwiseProduct(AWt);
			}
			YW = Y1.cwiseProduct(AWt);

			Am = (MW.transpose().conjugate() * MW).real();
			b = (MW.transpose().conjugate() * YW).real();

			solve(Am, b, par);

			for (auto i = 0; i <= nB; ++i)
			{
				B[i] = par(i);
			}
			for (auto i = 1; i <= nA; ++i)
			{
				A[i] = par(nB + i);
			}
		}

		return tuple<Vectord, Vectord>(B, A);
	}

	cplxVectord warpPolesFr(Vectord &W, cplxVectord &X, cplxVectord &Y, unsigned nPoles, double lambda, unsigned iter)
	{
		// Transforming frequency axis //
		//-----------------------------//
		Vectord warpedW(W.size());
		auto lambda2 = lambda * lambda;
		for (auto i = 0; i < W.size(); ++i)
		{
			warpedW[i] = atan2((1-lambda2)*sin(W[i]),(1+ lambda2)*cos(W[i])-2*lambda);
		}

		auto pars = lsidFr(warpedW, X, Y, nPoles, nPoles, iter);
		//auto pars = lsidFr(warpedW, X, Y, nPoles, nPoles, iter);

		auto pwp = roots(get<1>(pars));

		// Unwrapping //
		//------------//
		auto p(pwp);
		for (auto i = 0; i < p.size(); ++i)
			p[i] = (p[i] + lambda) / (1.0 + lambda * p[i]);

		return p;
	}


	cplxVectord dualWarpPolesFr(Vectord &W, cplxVectord &X, const Vectord &Y, double crossOverFreqHz,
		unsigned crossOverLengthSamples, double lambda1, double lambda2, unsigned nPoles1, unsigned nPoles2,
		double sampleRateHz, unsigned iter, bool useNAKspline)
	{
		auto HL = X.size();
		auto C = firstOfVal(W, 2 * M_PI * crossOverFreqHz / sampleRateHz);		//index of cross freq in vector W
		auto WL = 2 * int(round(crossOverLengthSamples / 2.0));

		if ((int(C) - WL / 2) <0)
			throw out_of_range("Low freq, cross freq and cross length are not compatible!");

		Vectord tmpHan(2 * WL + 1); hanning(tmpHan);
		Vectord window(&tmpHan[0], &tmpHan[WL]);

		Vectord WHF(X.size());	//instead of zeros and insert
		std::copy(window.begin(), window.end(), WHF.begin() + C - WL / 2+1);
		std::fill(WHF.begin() + C + WL / 2 +1, WHF.end(),1.0);		//instead of resize


		Vectord WLF(X.size(),1.);
		std::transform(WLF.begin(), WLF.end(), WHF.begin(), WLF.begin(), std::minus<>());	//1-WHF

		// windowing the input and output //
		//auto X1 = minPhaseN(VectorExtensions::operator*(real(X), WLF) + VectorExtensions::operator*(Vectord(HL, real(X[C])), WHF), W, pow(2, 15));	//only if Y is the target!!!
		//auto X2 = minPhaseN(VectorExtensions::operator*(real(X), WHF) + VectorExtensions::operator*(Vectord(HL, real(X[C])), WLF), W, pow(2, 15));
		// we can skip this only if desing a filter! not a compensation!

		Vectord Y1(Y); vectMult(Y1, WLF);
		Vectord Y2(Y); vectMult(Y2, WHF);


		multiplyVector(WHF, Y[C]);
		multiplyVector(WLF, Y[C]);

		std::transform(Y1.begin(), Y1.end(), WHF.begin(), Y1.begin(), std::plus<>());
		std::transform(Y2.begin(), Y2.end(), WLF.begin(), Y2.begin(), std::plus<>());

		cplxVectord Y1m(Y1.size()), Y2m(Y1.size());
		minPhaseN(Y1, W, Y1m,useNAKspline);
		minPhaseN(Y2, W, Y2m, useNAKspline);

		// first band //
		auto pWarp = warpPolesFr(W, X, Y1m, nPoles1, lambda1, iter);		
		// second band //
		auto pWarp2 = warpPolesFr(W, X, Y2m, nPoles2, lambda2, iter);	

		pWarp.insert(pWarp.end(), pWarp2.begin(), pWarp2.end());

		return pWarp;	
	}

	void minPhaseN(Vectord& Magn, Vectord& Wspec, cplxVectord& output, bool useNAKspline, size_t nPoints)
	{
		// resample to linear frequency scale
		Vectord LinW(nPoints + 1); linspace(LinW, 0, M_PI, LinW.size());
		//warp transform
		double lambda = -0.9;
		auto lambda2 = lambda * lambda;
		for (double& i : LinW)
		{
			i = atan2((1 - lambda2) * sin(i), (1 + lambda2) * cos(i) - 2 * lambda);
		}

		Vectord LinMagn(LinW.size());

		if(useNAKspline)
			absSplineNaK(Wspec, Magn, LinW, LinMagn);			//not a knot spline, but the results the same now
		else
			absSpline(Wspec, Magn, LinW, LinMagn);

		const auto N = LinMagn.size() - 1;

		// make LinMagn symmetric and periodic [reverse(LinMagn) LinMagn reverse(LinMagn) LinMagn] to avoid convolution probloms at boundaries
		//little cheaty but half the price [LinMagn reverse(LinMagn)]
		Vectord PerLinMagn(2*nPoints);
		std::copy(LinMagn.begin(), LinMagn.end()-1, PerLinMagn.begin());			// .. LinMagn
		std::copy(LinMagn.begin() + 1, LinMagn.end(), PerLinMagn.begin() + LinMagn.size() - 1);
		reverse(PerLinMagn.begin() + LinMagn.size()-1 , PerLinMagn.end());						

		// calculating the Phase with Hilbert transform and taking part corresponding to [0:pi]
		log(PerLinMagn);

		cplxVectord cplxPerLinMagn(PerLinMagn.size());
		hilbert(PerLinMagn, cplxPerLinMagn);
		
		Vectord LinPhase(N+1);
		negImag(cplxPerLinMagn, LinPhase, N+1);				//negative imag value

		Vectord Phase(Wspec.size());
		if (useNAKspline)
			splineNaK(LinW, LinPhase, Wspec, Phase);		//not a knot spline
		else
			spline(LinW, LinPhase, Wspec, Phase);

		// The minphase version of original signal is Magn.*exp(j*Phase)
		multiplyVector(Phase, complex<double>{0., 1.},output);
		for (auto i = 0; i < output.size(); ++i)
			output[i] = std::exp(output[i])* Magn[i];
	}

	void parFiltFreqResp(Vectord2D& Bm, Vectord2D& Am, Vectord& FIR, Vectord& w, cplxVectord& outResp)
	{
		int nSec = Bm.size();
		int nFIR = FIR.size();

		Eigen::VectorXcd Z(w.size());

		for (int i = 0; i < w.size(); i++)
		{
			Z(i) = exp(complex<double>(0.0, -w[i]));
		}

		Eigen::VectorXcd Z2 = Z.cwiseProduct(Z);

		Eigen::VectorXcd H(w.size());
		H.setZero();
		Eigen::VectorXcd num(w.size()), den(w.size());

		for (auto i = 0; i < nSec; ++i)
		{
			num = Bm[i][0] + (Bm[i][1] * Z).array();
			den = (Am[i][0] + (Am[i][1] * Z + Am[i][2] * Z2).array()).cwiseInverse();
			num = num.cwiseProduct(den);
			H += num;

		}

		for (auto i = 0; i < nFIR; ++i)
		{
			num = FIR[i] * Z.array().pow(i);
			H += num;
		}

		for (auto i = 0; i < w.size(); ++i)
			outResp[i] = H(i);
	}

	void parFiltFreqResp(Vectord2D& Bm, Vectord2D& Am, Vectord& FIR, Vectord& w, Vectord& outMagn)
	{
		int nSec = Bm.size();
		int nFIR = FIR.size();

		Eigen::VectorXcd Z(w.size());

		for (int i = 0; i < w.size(); i++)
		{
			Z(i) = exp(complex<double>(0.0, -w[i]));
		}

		Eigen::VectorXcd Z2 = Z.cwiseProduct(Z);

		Eigen::VectorXcd H(w.size());
		H.setZero();
		Eigen::VectorXcd num(w.size()), den(w.size());

		for(auto i = 0;i< nSec;++i)
		{
			num =  Bm[i][0] + (Bm[i][1] * Z).array();
			den = (Am[i][0] + (Am[i][1] * Z + Am[i][2] * Z2).array()).cwiseInverse();
			num = num.cwiseProduct(den);
			H += num;
			
		}

		for(auto i = 0;i<nFIR; ++i)
		{
			num = FIR[i] * Z.array().pow(i);
			H += num;
		}

		for (auto i = 0; i < w.size(); ++i)
			outMagn[i] = abs(H(i));
	}
}





