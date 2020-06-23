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
#ifndef _VECTOR_EXTENSIONS_H_
#define _VECTOR_EXTENSIONS_H_

#include <vector>
#include <complex>

typedef std::vector<std::complex<double>> cplxVectord;

namespace VectorExtensions
{

	inline std::ostream & operator << (std::ostream &out, const std::vector<double> &vec)
	{
		for (auto value : vec)
		{
			out << value << std::endl;
		}
		return out;
	}

	/// <summary> Multiplies given vector with coefficient changing the original. </summary>
	template<typename T>
	void multiplyVector(std::vector<T>& inputVector, T coefficient)
	{
		transform(inputVector.begin(), inputVector.end(), inputVector.begin(),
			bind(std::multiplies<T>(), std::placeholders::_1, coefficient));
	}
	/// <summary> Multiplies given vector with complex coeficient. </summary>
	template<typename T, typename U>
	void multiplyVector(std::vector<T>& inputVector, U coefficient, cplxVectord &outputVector)
	{
		for(auto i = 0;i< outputVector.size();++i)
		{
			outputVector[i] = inputVector[i] * coefficient;
		}
	}



	template<typename T>
	void vectMult (const std::vector<T> lhf, const std::vector<T> rhf, std::vector<T> &output)
	{
		transform(lhf.begin(), lhf.end(), rhf.begin(), output.begin(), std::multiplies<>());
	}

	/// <summary> Piecewise multiplication of vectors, but the first is changed.</summary>
	template<typename T>
	void vectMult(std::vector<T> &lhf, const std::vector<T> &rhf)
	{
		transform(lhf.begin(), lhf.end(), rhf.begin(), lhf.begin(), std::multiplies<>());
	}


	/// <summary> Calculates logarithm of given vector changing the original. </summary>
	template <typename T>
	void log(std::vector<T>& inputVector)
	{
		for (auto i = 0; i <inputVector.size(); ++i)
			inputVector[i] = std::log(inputVector[i]);
		
	}


	/// <summary> Calculates exp() of given vector, changing the original. </summary>
	template <typename T>
	void exp(std::vector<T>& inputVector)
	{
		for (auto i = 0; i <= inputVector.size(); ++i)
			inputVector[i] = std::exp(inputVector[i]);
	}


	/// <summary> Returns real value of given vector </summary>
	inline std::vector<double> real(cplxVectord& inputVector)
	{
		std::vector<double> result(inputVector.size());
		for (auto i = 0; i < inputVector.size(); ++i)
			result[i] = real(inputVector[i]);
		return result;
	}

	/// <summary> Returns imag value of given vector </summary>
	inline void negImag(cplxVectord& inputVector, Vectord &outputVector, int howMany)
	{
		for (auto i = 0; i < howMany; ++i)
			outputVector[i] = -imag(inputVector[i]);
	}

	inline void posImag(cplxVectord& inputVector, Vectord& outputVector, int howMany)
	{
		for (auto i = 0; i < howMany; ++i)
			outputVector[i] = imag(inputVector[i]);
	}

	/// <summary> Returs index of first value greater than parameter. </summary>
	template <typename T>
	inline unsigned int firstOfVal(const std::vector<T>& inpVector, T value)
	{
		return (lower_bound(inpVector.begin(), inpVector.end(), value) - inpVector.begin());
	}


}


#endif
