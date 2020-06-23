// Parallel filter bank. Based on Balázs' Bank approach for modeling desired frequency response.
// Second order sections in biquad form without B2 coefficient.
// 
// Matouš Vrbík
// matousvrbik[at]gmail.com
// 04/2020


#pragma once
#include "denormals.h"
#include <deque>
#include <memory>

class FilterBank
{
public:
	FilterBank(float *A_coeffs1, float *A_coeffs2, float *B_coeffs1, float *B_coeffs2, int numFilters, float *FIR_coeffs, int FIR_size, float outGain);
	FilterBank(const FilterBank& filter_bank);
	~FilterBank();

	/// Procesess input sample throught parallel filters.
	float processSample(float input);

	/// Resets state variables and input buffer to zeros.
	void resetMemory();

	/// Sets all coeffs pointers to nullptr and clears state variables and buffers.
	void Reset();

	/// Sets all the coeffs pointers and filter dimensions. Also clears the current state variables
	void SetCoeffsPtrs(float *A_coeffs1, float *A_coeffs2, float *B_coeffs1, float *B_coeffs2, int numFilters, float *FIR_coeffs, int FIR_size, float outGain);


private:
	float *A1, *A2, *B0, *B1, *FIR, outGain;			//pointers to data, either static from header, or loaded from file. Mostly pointing to SimulationCoefficients/Static/ object.
	int numFilters, FIR_size;
	float * d1, *d2, *out_tmp;			// canonical filter memory
	std::deque<float> input_mem;

	float CanonicalForm(float input);

};

inline FilterBank::FilterBank(float* A_coeffs1, float* A_coeffs2, float* B_coeffs1, float* B_coeffs2,
	int numFilters, float* FIR_coeffs, int FIR_size, float outGain)
{
	A1 = A_coeffs1;
	B0 = B_coeffs1;
	A2 = A_coeffs2;
	B1 = B_coeffs2;
	this->numFilters = numFilters;
	this->FIR_size = FIR_size;

	//init delay registers
	
	d1 = new float[numFilters]();
	d2 = new float[numFilters]();
	out_tmp = new float[numFilters]();

	FIR = FIR_coeffs;
	this->outGain = outGain;
	for (int i = 0; i < FIR_size - 1; ++i)
		input_mem.push_front(0.f);
}

inline FilterBank::FilterBank(const FilterBank& filter_bank)
{
	A1 = filter_bank.A1;
	B0 = filter_bank.B0;
	A2 = filter_bank.A2;
	B1 = filter_bank.B1;
	numFilters = filter_bank.numFilters;
	FIR_size = filter_bank.FIR_size;

	//init delay registers
	d1 = new float[numFilters]();
	d2 = new float[numFilters]();
	out_tmp = new float[numFilters]();


	FIR = filter_bank.FIR;
	outGain = filter_bank.outGain;
	for (int i = 0; i < FIR_size - 1; ++i)
		input_mem.push_front(0.f);

}

inline FilterBank::~FilterBank()
{
	delete [] d1;
	delete [] d2;
	delete [] out_tmp;
}

inline float FilterBank::processSample(float input)
{
	//return DirectForm(input);
	return CanonicalForm(input);
}


inline void FilterBank::resetMemory()
{
	memset(d1, 0, sizeof(float)*numFilters);
	memset(d2, 0, sizeof(float)*numFilters);

	input_mem.clear();
	input_mem.resize(FIR_size);		//snad to tam dá nuly
}

inline void FilterBank::Reset()
{
	B0 = nullptr;
	B1 = nullptr;
	A1 = nullptr;
	A2 = nullptr;
	FIR = nullptr;
	outGain = 0.f;
	numFilters = FIR_size = 0;
	resetMemory();
}

inline void FilterBank::SetCoeffsPtrs(float* A_coeffs1, float* A_coeffs2, float* B_coeffs1, float* B_coeffs2,
	int numFilters, float* FIR_coeffs, int FIR_size, float outGain)
{
	A1 = A_coeffs1;
	B0 = B_coeffs1;
	A2 = A_coeffs2;
	B1 = B_coeffs2;
	this->numFilters = numFilters;
	this->FIR_size = FIR_size;

	//init delay registers
	d1 = new float[numFilters]();
	d2 = new float[numFilters]();
	out_tmp = new float[numFilters]();


	FIR = FIR_coeffs;
	this->outGain = outGain;
	resetMemory();
}


inline float FilterBank::CanonicalForm(float input)
{
	//
	// x(n) ---> B0 ---> Σ -----------------> y(n)
	//       |           ^			  |
	//	     |         [n-1] - d1	  |
	//		 |			 |			  |
	//		 |			 |			  |
	//		 --> B1 ---> Σ <--- A1 <---
	//					 ^			  |
	//					 |			  |
	//				   [n-1] - d2	  |
	//					 |			  |
	//					 Σ <--- A2 <---
	//
	// A1, A2 must be added negatively, if the coeffs are not already negated

	undenormalise(input);

	for (int i = 0; i < numFilters; i++)
	{
		out_tmp[i] = B0[i] * input + d1[i]; // out = b0 * in + d1;
	}
	for (int i = 0; i < numFilters; i++)
	{
		d1[i] = B1[i] * input + d2[i];		// d1 = b1 * in + d2;
	}
	for (int i = 0; i < numFilters; i++)
	{
		d1[i] -= (A1[i] * out_tmp[i]);      // d1 += a1 * out;	
	}
	for (int i = 0; i < numFilters; i++)
	{
		d2[i] = (-A2[i] * out_tmp[i]);		// d2 = a2 * out;
	}
	for (int i = 1; i < numFilters; i++)	// sum all filters outputs
	{
		out_tmp[0] += out_tmp[i];
	}


	if (FIR != nullptr)
	{
		float firOut = input * FIR[0];

		for (int i = 1; i < FIR_size; i++)		//for FIR of length <1, NOT TESTED
		{
			firOut += FIR[i] * input_mem[FIR_size - i - 1];		//in reverse, bcs its deque
		}
		if (FIR_size > 1)
		{
			input_mem.push_back(input);		//if FIR_size == 1 it pops the same as the pushed, so it makes no sense, and sometimes it causes undefined behavior in xutility orphan_all()
			input_mem.pop_front();
		}
		out_tmp[0] += firOut;
	}
	return out_tmp[0] * outGain;

}
