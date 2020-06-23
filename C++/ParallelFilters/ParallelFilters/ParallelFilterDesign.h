/*
* Copyright (c) 2020 Matouš Vrbík
*
* URL	https://github.com/TheMates/Diploma-thesis-public
* 06/2020
* matousvrbik[at]gmail.com
*
* Distributed under CC BY-NC-SA 4.0 licence see LICENCE for details.
*/

#pragma once
#ifndef _PARALLEL_FILTER_DESIGN_H_
#define _PARALLEL_FILTER_DESIGN_H_ 


#include <complex>
#include <string>
#include <vector>
#include <iterator>
#include "Extensions.h"


namespace ParallelFilters
{

	/// <summary> 
	/// Pole set generation for parallel filter design. \n
	/// [P] = freqpoles(fr, Fs) creates poles from the given frequency vector 
	///	   FR with sampling frequency Fs such that the frequency responses  \n
	///	   of the second - order basis functions cross approximatelly at their -3 dB point. \n
	///	 (See the equations in my 128th AES Convention paper.) \n
	///	 \n
	///	[P] = freqpoles(fr, Fs, Q) creates poles from the given frequency vector 
	///	   FR with sampling frequency Fs but it uses the quality factors \n
	///	   given in Q vector to set the radius of the poles. (-3dB bandwidth is DeltaF = FR / Q.) \n
	///	 \n
	///	   http://www.mit.bme.hu/~bank/parfilt \n
	/// </summary>
	/// <param name="fr"> Set of frequencies in Hz </param>
	/// <param name="output"> Output vector. </param>
	/// <param name="Fs"> Sample rate </param>
	/// <param name="Q"> Quality factor  </param>
	void freqpoles(const Vectord &fr, cplxVectord& output, const double Fs = 44100,const double Q = -1);

	/// <summary> Solves equation system Ax=b for x. Utilizing Eigen library and LDLT solver. </summary>
	/// <param name="A"> Matrix of left side. </param>
	/// <param name="b"> Matrix of right side. </param>
	/// <param name="x"> Solution parameters. </param>
	/// <returns> Solution x. </returns>
	void solve(Eigen::MatrixXd &A, Eigen::MatrixXd&b, Eigen::MatrixXd &x);


	/// <summary> Direct design of second-order parallel filters for a given pole set in the frequency domain. \n
	///This function is based on the paper
	/// \n
	///	Balázs Bank, "Logarithmic Frequency Scale Parallel Filter Design with Complex
	///	and Magnitude - Only Specifications, " IEEE Signal Processing Letters,
	///	Vol. 18, No. 2, pp. 138 - 141, Feb. 2011.
	/// \n
	///	http://www.mit.bme.hu/~bank/parfilt
	/// </summary>
	/// <param name="W"> Angular frequencies in radians from 0 to pi. </param>
	/// <param name="H"> Transfer function. </param>
	/// <param name="poles"> Poles of filters. Expecting result from freqpoles() function. </param>
	/// <param name="nFIR"> The number of taps in parralel FIR filter. Default is 1 - in this case the FIRcoeff is simply gain. </param>
	/// <returns> Tuple<Vectord <b>Bm</b>, Vectord <b>Am</b>, Vectord <b>FIRcoeff</b>> where <b>Bm</b> are numerator coefficients of second-order filters and <b>Am</b> are their denominators. <b>FIRcoeff</b> are taps of FIR filter. </returns>
	std::tuple<std::vector<Vectord>, std::vector<Vectord>,Vectord> parFiltDesignFr
	(Vectord W, cplxVectord H, cplxVectord poles, unsigned int nFIR = 1);	//in matlab parameter Wt, but not doing it now


	/// <summary> 
	///  LSIDFR[B, A] = lsidfr(W, X, Y, NB, NA, ITER, WEIGHT); estimates a discrete - time
	///	system which produces an output Y
	///	to the input X, given at angular frequencies W, ranging from 0 to pi(one sided specification).
	///	The numerator and denominator are given in B and A, and their order
	///	is set by NB and NA, respectively.
	/// \n\n
	///	The estimation is done in the frequency domain.If ITER = 0 (or , not given),
	///	we minimize the equation error :
	/// | Y * A - B * X | ^ 2
	///	If ITER > 0, we make an iterative weighting, corresponding to the
	///	frequency - domain Steiglitz - McBride algorithm, minimizing the true
	///	error
	///	| Y - B / A * X | ^2.
	/// \n\n
	///	If X is a scalar, then we simply design a filter for the
	///	specification H = Y / X, where the code extends the dimension of X to
	///	match that of Y. (Typically, use X = 1 and Y = H as the filter specification.)
	/// \n\n
	///	If WEIGHT is given, we minimize
	///	WEIGHt* | Y * A - B * X | ^ 2
	///	or
	///	WEIGHT * | Y - B / A * X | ^2.
	/// \n\n
	///	For the theory, see
	///	L.B.Jackson, “Frequency - domain Steiglitz - McBride method for least -
	///	squares filter design, ARMA modeling, and periodogram smoothing, ”
	///	IEEE Signal Process.Lett., vol. 15, pp. 49–52, 2008.
	/// </summary>
	/// <param name="W"> Angular frequencies in radians from 0 to pi. </param>
	/// <param name="X"> Input transfer function.  </param>
	/// <param name="Y"> Output transfer function. </param>
	/// <param name="nB"> Order of numerator polynom. </param>
	/// <param name="nA"> Order of denominator polynom. </param>
	/// <param name="iter"> Number of iterations. </param>
	/// <returns> tuple<Vectord <b>B</b>, Vectord <b>A</b>> , where <b>B</b> is the numerator and <b>A</b> is the denominator of discrete time system, which produces the output Y. </returns>
	std::tuple < Vectord, Vectord> lsidFr(Vectord& W, cplxVectord& X, cplxVectord& Y, unsigned int nB, unsigned int nA, unsigned int iter = 5);	//no weights now


	/// <summary> 
	///  Pole set generation for parallel filters based on a warped IIR design.
	///	[p] = warppolesfr(W, X, Y, PNUM, LAMBDA, ITER, Wt) designs a warped IIR filter of the
	///	order PNUM based on the input and output spectrum X and Y specified at
	///	W angular frequencies(ranging from 0 to pi) with the warping parameter
	///	LAMBDA.Then it finds the poles and maps them back to linear frequency
	///	resolution.The poles are given in vector P and can be directly used by
	///	the PARFILTDESFR or PARFILTIDFR commands.
	/// \n\n
	///	The IIR filter design is based on the frequency - domain
	///	Steiglitz - McBride algorithm(see the help of LSIDFR) : ITER is the
	///	number of iterations(default is 0, meaning equation error
	///		minimization), and Wt gives weights to the frequency points.There last
	///	two parameters are optional.
	/// \n\n
	///	For normal filter design, set X = 1 and Y as the target frequency response of the
	///	filter.For equalizer design, X is the system response, while Y is the
	///	target response. (This is basically a system identification between input X
	///		and output Y).
	/// \n\n
	///	For related publications and matlab code see
	///	http ://www.mit.bme.hu/~bank/parfilt
	/// </summary>
	/// <param name="W"> Vector of angular frequencies. </param>
	/// <param name="X"> Transfer function of input. </param>
	/// <param name="Y"> Transfer function of output. </param>
	/// <param name="nPoles"> Number of poles to be generated. </param>
	/// <param name="lambda"> Lambda warping parameter. </param>
	/// <param name="iter"> Number of iterations of optimalization. </param>
	/// <returns> Set of poles desined in warped frequency resolution. </returns>
	cplxVectord warpPolesFr(Vectord &W, cplxVectord &X, cplxVectord &Y, unsigned int nPoles, double lambda, unsigned int iter = 5);

	/// <summary> 
	///  Pole set generation for parallel filter design based on
	///	dual - band warped filter design.
	/// \n\n
	///	[p] = dualwarppolesfr(W, X, Y, CrossFr, CrossFade, Lambda1, Lambda2, P1, P2);
	/// designs two warped IIR filters  based on the input and output spectrum
	///	X and Y specified at W angular frequencies(ranging from 0 to pi).
	/// \n\n
	///	The design is based on splitting the frequency range to two parts and
	///	designing separate warped IIR filters with Lambda1 and Lambda2 warping
	///	parameters, and P1 and P2 orders.
	/// \n\n
	///	The crossover frequency between the two bands is given by the CrossFr parameter,
	///	which is the crossover frequency.
	///	The out - of band regions of X and Y are crossfaded to constant gain, the width
	///	of the crossfade region is defined by CrossFade in samples of the the W
	///	vector. Thus, the crossfade starts at angular frequency
	///	W(CrossFr - CrossFade / 2) and stops at W(CrossFr + CrossFade / 2).
	/// \n\n
	///	After separating the bands, they are made minimum - phase(this is
	///		necassary because we crossfaded to constant gain - and there are no
	///		drawbacks, since we are interested in the poles only), and warped IIR
	///	filters are designed using WARPPOLESFR. Finally, the pole sets are united.
	/// \n\n
	///	For normal filter design, set X = 1 and Y as the target frequency response of the
	///	filter.For equalizer design, X is the system response, while Y is the
	///	target response. (This is basically a system identification between input X and output Y).
	/// \n\n
	///	This code is based on the paper
	///	Balázs Bank and Germán Ramos, "Improved Pole Positioning for Parallel Filters
	///	Based on Spectral Smoothing and Multiband Warping, " IEEE Signal Processing Letters,
	///	vol. 18, no. 5, pp. 299 - 302, Mar. 2011.
	/// \n\n
	///	http://www.mit.bme.hu/~bank/parfilt
	/// </summary>
	/// <param name="W"> Vector of angular frequencies. </param>
	/// <param name="X"> Transfer function of input. </param>
	/// <param name="Y"> Transfer function of output. Not minphase. </param>
	/// <param name="crossOverFreqHz"> Crossover frequency in Hz. </param>
	/// <param name="crossOverLengthSamples"> Length of window for crossover in samples. </param>
	/// <param name="nPoles1"> Number of poles in low section. </param>
	/// <param name="nPoles2"> Number of poles in high section. </param>
	/// <param name="lambda1"> Lambda warping parameter in low section. </param>
	/// <param name="lambda2"> Lambda warping parameter in high section. </param>
	/// <param name="sampleRateHz"> Sampling frequency. </param>
	/// <param name="iter"> Number of iterations of optimalization. </param>
	/// <param name="useNAKspline"> If minphasen uses not-a-knot extrapolation or linear. </param>
	/// <returns> Set of poles designed in warped frequency resolution. </returns>
	cplxVectord dualWarpPolesFr(Vectord &W, cplxVectord &X, const Vectord &Y, double crossOverFreqHz, unsigned int crossOverLengthSamples,
		double lambda1, double lambda2, unsigned int nPoles1, unsigned int nPoles2, double sampleRateHz = 44100.0, unsigned int iter = 0, bool useNAKspline = false);


	/// <summary> 
	/// Calculates the transfer function H = Magn * exp(j*Phase)
	///	of a minimum phase system from its magnitude by Hilbert
	///	transform.
	/// \n\n
	///	Phase = -Hilbert{ ln(Magn) }
	/// \n\n
	///	H = minphasen(Magn) calculates the transfer function "H" from the
	///	specification "Magn", where the frequency points are
	///	assumed to be linearly spaced between 0 and pi,
	///	i.e., W = pi * [0:length(Magn)]. / length(Magn)
	/// \n\n
	///	H = minphasen(Magn, Wspec, LinN) calculates the transfer function "H" at
	///	the frequency points "Wspec" from the magnitude values "Magn"
	///	at the same frequency points. "Wspec" should begin with 0
	///	and end with pi.The magnitude on the linear frequency
	///	scale(which is needed for the Hilbert transform) is calculated
	///	by cubic spline interpolation, and LinN linear frequency points
	///	are used.Default is LinN = 2 ^ 14.
	/// </summary>
	/// <param name="Magn">Magnitude of transfer function. In abs value, not in dB.</param>
	/// <param name="Wspec"> Vector of frequencies to be calculated. </param>
	/// <param name="output"> Preallocated complex output vector. </param>
	/// <param name="useNAKspline"> Use not-a-knot spline extrapolation. Otherwise linear extrapolation. </param>
	/// <param name="nPoints"> Number of points to calculate the linear frequency points for Hilbert transfrom. </param>
	/// <returns> Complex minimum phase transfer function. </returns>
	void minPhaseN(Vectord &Magn, Vectord &Wspec, cplxVectord & output, bool useNAKspline, size_t nPoints = 1024);


	/// <summary> 
	///PARFILTDFRESP - Frequency response of second - order parallel filters.
	/// \n
	///	H = parfiltfresp(Bm, Am, FIR, W, out); computes the frequency response of the
	///	parallel filter with[Bm, Am] second - order section coefficients
	///	and the coefficients of the FIR part(FIRcoeff) at given W frequencies
	///	(in radians from 0 to pi).
	/// \n\n
	///	http://www.mit.bme.hu/~bank/parfilt
	/// </summary> 
	/// <param name="Bm">Numerator coefficients of the filter.</param>
	/// <param name="Am">Denominator coefficients of the filter.</param>
	/// <param name="FIR">FIR coefficients of the filter.</param>
	/// <param name="w">Angular frequencies.</param>
	/// <param name="outResp">Output frequency response of the filter response.</param>
	void parFiltFreqResp(Vectord2D& Bm, Vectord2D& Am, Vectord& FIR, Vectord& w, cplxVectord& outResp);
	/// <summary> 
	///PARFILTDFRESP - Frequency response of second - order parallel filters, only magnitude.
	/// \n
	///	H = parfiltfresp(Bm, Am, FIR, W, out); computes the frequency response of the
	///	parallel filter with[Bm, Am] second - order section coefficients
	///	and the coefficients of the FIR part(FIRcoeff) at given W frequencies
	///	(in radians from 0 to pi).
	/// \n\n
	///	http://www.mit.bme.hu/~bank/parfilt
	/// </summary> 
	/// <param name="Bm">Numerator coefficients of the filter.</param>
	/// <param name="Am">Denominator coefficients of the filter.</param>
	/// <param name="FIR">FIR coefficients of the filter.</param>
	/// <param name="w">Angular frequencies.</param>
	/// <param name="outMagn">Output magnitude of the filter response.</param>
	void parFiltFreqResp(Vectord2D& Bm, Vectord2D& Am, Vectord& FIR, Vectord& w, Vectord& outMagn);


}

#endif 
