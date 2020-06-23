#ifndef _PARALLELFILTER_EXTENSIONS_H_
#define _PARALLELFILTER_EXTENSIONS_H_

#include "simple_fft/fft.h"		//Dmitry Ivanov https://github.com/d1vanov/Simple-FFT
#include "Eigen/Eigen"			//https://gitlab.com/libeigen/eigen
#include "rapidcsv.h"			//Kristofer Berggren https://github.com/d99kris/rapidcsv


typedef std::vector<std::complex<double>> cplxVectord;
typedef std::vector<std::vector<std::complex<double>>> cplxVectord2D;
typedef std::vector<std::vector<double>> Vectord2D;
typedef std::vector<double> Vectord;


namespace ParallelFilters
{
	namespace Extensions
	{
		/// <summary> Parses string to complex number. Accepts 'a+bi' or 'a, bi'. 'i' or 'j' can be used.</summary>
		/// <param name="inputString"> Complex number in string format. </param>
		/// <returns> Complex number of template type. </returns>
		template<typename T>
		std::complex<T> parseStringToComplex(const std::string& inputString)
		{
			return std::complex<T>(std::stod(inputString.substr(0, inputString.find_first_of("+,"))),
				std::stod(inputString.substr(inputString.find_first_of("+,") + 1, inputString.find_first_of("ij") - 1)));
		}

		/// <summary> Loads vector stored as column in csv file to output vector. </summary>
		/// <param name="filePath"> Path of file. </param>
		/// <param name="output"> Output vector. </param>
		void loadRealFromCsvVec(const std::string filePath, Vectord & output);

		/// <summary> Loads vector stored as columns in csv file to output vector. </summary>
		/// <param name="filePath"> Path of file. </param>
		/// <param name="output"> Output vector. </param>
		void loadRealFromCsvVec2D(const std::string filePath, Vectord2D& output);


		/// <summary> Writes 2D vector to csv file. </summary>
		/// <param name="inpVector"> Input vector. </param>
		/// <param name="path"> Path of file to be saved including ".csv". </param>
		void writeToCsv(cplxVectord2D &inpVector, std::string path);
		/// <summary> Writes 2D vector to csv file. </summary>
		/// <param name="inpVector"> Input vector. </param>
		/// <param name="path"> Path of file to be saved including ".csv". </param>
		void writeToCsv(Vectord2D &inpVector, std::string path);
		/// <summary> Writes vector to csv file. </summary>
		/// <param name="inpVector"> Input vector. </param>
		/// <param name="path"> Path of file to be saved including ".csv". </param>
		void writeToCsv(Vectord &inpVector, std::string path);
		/// <summary> Writes complex vector to csv file. First column real, second column imag. </summary>
		/// <param name="inpVector"> Input vector. </param>
		/// <param name="path"> Path of file to be saved including ".csv". </param>
		void writeToCsv(cplxVectord& inpVector, std::string path);



		/// <summary> Calculates next 2^n greater than given length. </summary>
		/// <param name="length"> Length of signal </param>
		/// <returns> Next greater 2^n number. </returns>
		static size_t getNfft(unsigned length)
		{
			const unsigned int originalLength = length;
			unsigned int bit = 1;

			while (length >>= 1)
				bit <<= 1;
			//returns highest set bit
			return originalLength > bit ? 2 * bit : bit;
		}

#pragma region MATLAB

		/// <summary> Function returns vector of lineary spaced numbers between given borders. </summary>
		/// <param name="output"> Preallocated output </param>
		/// <param name="min"> First number </param>
		/// <param name="max"> Last Number </param>
		/// <param name="n"> Size of vector </param>
		void linspace(Vectord& output, double min, double max, unsigned int n = 101);

		/// <summary> Function returns vector of logarithmicaly spaced numbers between given borders. </summary>
		/// <param name="min"> First number </param>
		/// <param name="max"> Last Number </param>
		/// <param name="n"> Size of vector </param>
		/// <param name="base"> Base of logarithm </param>
		/// <returns> Vector of logarithmicaly spaced numbers between given borders. </returns>
		Vectord logspace(double min, double max, int n = 101, double base = 10.0f);
		/// <summary> Function returns vector of double logarithmicaly spaced numbers between given borders. Function is base^(base^t)) where t is lineary spaced vector between double logarithm of given borders. </summary>
		/// <param name="min"> First number </param>
		/// <param name="max"> Last Number </param>
		/// <param name="n"> Size of vector </param>
		/// <param name="base"> Base of logarithm </param>
		/// <returns> Vector of double logarithmicaly spaced numbers between given borders. </returns>
		Vectord loglogspace(double min, double max, int n = 101, double base = 10.0f);


		/// <summary> Creates hanning window of given length of the input vector. Changing the input vector. </summary>
		/// <param name="input">Preallocated input.</param>
		void hanning(Vectord& input);

		/// <summary> Calculates hilbert transfrom of input vector, which should be 2^n lenth. If not the length, padded with zeros. 
		/// \n\n
		///	FFT of input signal is calculated. Values between 1 : n/2 are multiplied by 2. Values between n/2+1 : n are set to 0. IFFT is calculated. 
		/// </summary>
		void hilbert(Vectord& inputVector, cplxVectord &outputVector);


		/// <summary> Tries to find a complex conjugate match for given complex number from complex numbers container. </summary>
		/// <param name="first"> Iterator to the first element in container to search. </param>
		/// <param name="last"> Iterator to the last element in container to search. </param>
		/// <param name="match"> Complex number, which complex conjugate is to be found. </param>
		/// <returns> Iterator to element, which is complex conjugate to match. If no match is found, throws out_of_range error. </returns>
		cplxVectord::iterator findCplxPair(cplxVectord::iterator first, cplxVectord::iterator last,
			const std::complex<double> match);

		/// <summary> Sorts complex numbers into complex conjugate pairs, grouping together complex conjugate pairs. Conjugate pairs are ordered 
		/// by increasing real part. Within a pair, the element with negative imaginary part comes first.</summary>
		/// <param name="poles"> Set of poles. </param>
		/// <returns> Reordered set of poles. </returns>
		void cplxpair(cplxVectord& poles);

		/// <summary> Returns vector of polynomial coefficients in order from highest exponent to lowest. </summary>
		/// <param name="roots"> Vector of roots. </param>
		/// <returns> Vector of polynomial coefficients.  </returns>
		cplxVectord poly(cplxVectord& roots);
		/// <summary> Returns vector of roots of polynomial. </summary>
		/// <param name="coeffs"> Vector of polynomial coefficients in order from highest exponent to lowest. Utilizing Eigen library - Eigen values of matrix. </param>
		/// <returns> Vector of polynomial roots.  </returns>
		cplxVectord roots(cplxVectord& coeffs);
		/// <summary> Returns vector of roots of polynomial. </summary>
		/// <param name="coeffs"> Vector of polynomial coefficients in order from highest exponent to lowest. Utilizing Eigen library - Eigen values of matrix. </param>
		/// <returns> Vector of polynomial roots.  </returns>
		cplxVectord roots(Vectord& coeffs);

		/// <summary> Polynomial evaluation at given values. </summary>
		/// <param name="coeffs"> Vector of polynomial coefficients in order from highest exponent to lowest. </param>
		/// <param name="values"> Vector of values at which the polynomial will be evaluated. </param>
		/// <returns> Values of polynomial at given points. </returns>
		Vectord polyval(Vectord& coeffs, Vectord& values);
		/// <summary> Polynomial evaluation at given values. </summary>
		/// <param name="coeffs"> Vector of polynomial coefficients in order from highest exponent to lowest. </param>
		/// <param name="values"> Vector of values at which the polynomial will be evaluated. </param>
		/// <returns> Values of polynomial at given points. </returns>
		cplxVectord polyval(cplxVectord& coeffs, cplxVectord& values);
		/// <summary> Polynomial evaluation at given values. </summary>
		/// <param name="coeffs"> Vector of polynomial coefficients in order from highest exponent to lowest. </param>
		/// <param name="values"> Vector of values at which the polynomial will be evaluated. </param>
		/// <returns> Values of polynomial at given points. </returns>
		cplxVectord polyval(Vectord& coeffs, cplxVectord& values);
		/// <summary> Polynomial evaluation at given values. </summary>
		/// <param name="coeffs"> Vector of polynomial coefficients in order from highest exponent to lowest. </param>
		/// <param name="values"> Vector of values at which the polynomial will be evaluated. </param>
		/// <returns> Values of polynomial at given points. </returns>
		cplxVectord polyval(cplxVectord& coeffs, Vectord& values);
		/// <summary> Polynomial evaluation at given values. </summary>
		/// <param name="coeffs"> Vector of polynomial coefficients in order from highest exponent to lowest. </param>
		/// <param name="input"> Vector of values at which the polynomial will be evaluated. </param>
		/// <param name="output"> Output vector that will be filled with values of polynomial at given points. </param>
		void polyval(Vectord& coeffs, Eigen::Map<Eigen::VectorXcd>& input, Eigen::VectorXcd& output);



		/// <summary> Calculates cubic spline interpolation from vector X and Y in given points. 
		/// \n\n
		/// Utilizing spline library https://github.com/ttk592/spline/ </summary>
		void spline(Vectord &X, Vectord &Y, Vectord &points, Vectord &output);
		/// <summary> Calculates not-a-knot cubic spline interpolation from vector X and Y in given points. 
		/// \n\n
		/// Utilizing Eigen library. </summary>
		void splineNaK(Vectord& X, Vectord& Y, Vectord& points, Vectord& output);
		/// <summary> Calculates cubic spline interpolation from vector X and Y in given points and also does abs() of the value - local speedup.
		/// \n\n
		/// Utilizing spline library https://github.com/ttk592/spline/ </summary>
		void absSpline(Vectord& X, Vectord& Y, Vectord& points, Vectord& out);
		/// <summary> Calculates cubic spline interpolation from vector X and Y in given points and also does abs() of the value - local speedup.
		/// \n\n
		/// Utilizing Eigen library. </summary>
		void absSplineNaK(Vectord& X, Vectord& Y, Vectord& points, Vectord& out);


		/// <summary> Makes the imaginary part negative. </summary>
		/// <param name="matrix"> Input matrix. </param>
		void conj(cplxVectord2D& matrix);

#pragma endregion


	}
}

#endif