/*
* Copyright (c) 2020 Matouš Vrbík
*
* URL	https://github.com/TheMates/Diploma-thesis-public
* 06/2020
* matousvrbik[at]gmail.com
*
* Distributed under MIT licence see LICENCE for details.
*/

#define NOMINMAX	//Windows.h defines min and max as well as std - so it messes everything up
#include "Extensions.h"
#include "VectorExtensions.h"
#include <fstream>
#include <iomanip>
#include "Windows.h"
#include "spline.h"		//https://github.com/ttk592/spline/

#include "SplineNaK.h"


using namespace std;
using namespace ParallelFilters::Extensions;



void ParallelFilters::Extensions::loadRealFromCsvVec(const std::string filePath, Vectord & output)
{
	try
	{
		rapidcsv::Document doc(filePath, rapidcsv::LabelParams(-1, -1));			//so first column and row are counted as data, not headers
		output = doc.GetColumn<double>(0);
	}
	catch (exception &e)
	{
		throw e.what();
	}
}

void ParallelFilters::Extensions::loadRealFromCsvVec2D(const std::string filePath, Vectord2D& output)
{
	char sep;
	GetLocaleInfoA(LOCALE_USER_DEFAULT, LOCALE_SLIST, &sep, sizeof(sep));
	try
	{
		rapidcsv::Document doc(filePath, rapidcsv::LabelParams(-1, -1), rapidcsv::SeparatorParams(sep));			//so first column and row are counted as data, not headers
		int cols = doc.GetColumnCount();
		output = Vectord2D(cols);
		for (auto col = 0; col < cols; ++col)
			output[col] = doc.GetColumn<double>(col);
	}
	catch (exception& e)
	{
		throw;
	}
}

void ParallelFilters::Extensions::writeToCsv(cplxVectord2D& inpVector, std::string path)
{
	ofstream myfile;
	//choose region dependent separator
	char sep;
	GetLocaleInfoA(LOCALE_USER_DEFAULT, LOCALE_SLIST, &sep, sizeof(sep));
	try
	{
		myfile.open(path);
		myfile << std::fixed << std::setprecision(15);
		for (auto row = 0;row< inpVector.size();++row)
		{
			for (auto col = 0;col< inpVector[0].size();++col)
			{
				myfile << noshowpos << inpVector[row][col].real() << showpos << inpVector[row][col].imag() << noshowpos << "i"<<sep;
			}
			myfile << "\n";
		}
	}
	catch (exception &e)
	{
		cout << e.what();
	}
	
	myfile.close();
}
void ParallelFilters::Extensions::writeToCsv(Vectord2D& inpVector, const std::string path)
{
	ofstream myfile;
	//choose region dependent separator
	char sep;
	GetLocaleInfoA(LOCALE_USER_DEFAULT, LOCALE_SLIST, &sep, sizeof(sep));

	try
	{
		myfile.open(path);
		myfile << std::fixed << std::setprecision(15);
		for (auto row = 0; row < inpVector.size(); ++row)
		{
			for (auto col = 0; col < inpVector[0].size(); ++col)
			{
				myfile << inpVector[row][col];
				if(col != inpVector[0].size()-1)
					myfile << sep;
			}
			myfile << "\n";
		}
	}
	catch (exception &e)
	{
		cout << e.what();
	}

	myfile.close();
}
void ParallelFilters::Extensions::writeToCsv(Vectord& inpVector, const std::string path)
{
	ofstream myfile;

	try
	{
		myfile.open(path);
		myfile << std::fixed << std::setprecision(15);
		for (auto row = 0; row < inpVector.size()-1; ++row)
		{
			myfile << inpVector[row] << "\n";
		}
		myfile << inpVector[inpVector.size() - 1];	//last one without separator
	}
	catch (exception &e)
	{
		cout << e.what();
	}

	myfile.close();
}
void ParallelFilters::Extensions::writeToCsv(cplxVectord& inpVector, std::string path)
{
	ofstream myfile;

	try
	{
		myfile.open(path);
		myfile << std::fixed << std::setprecision(15);
		for (auto row = 0; row < inpVector.size(); ++row)
		{
			myfile << inpVector[row].real() << "," << inpVector[row].imag();
			myfile << "\n";
		}
	}
	catch (exception& e)
	{
		cout << e.what();
	}

	myfile.close();

}


#pragma region MATLAB

void ParallelFilters::Extensions::linspace(Vectord& output, double min, double max, unsigned n)
{
	if (max < min) swap(min, max);
	const double delta = (max - min) / static_cast<float>(n - 1);
	//for (auto i = 0; i < n; i++)
	//	output[i] = min + i * delta;
	std::vector<double>::iterator o;
	double val;
	for (o = output.begin(), val = min; o != output.end(); ++o, val += delta)
	{
		*o = val;
	}
}

Vectord ParallelFilters::Extensions::logspace(double min, double max, int n, double base)
{
	if (base <= 0.0) throw out_of_range("Parameter base must be a positive number!");
	if (n < 0) throw out_of_range("Parameter n must be a positive number!");
	if (max < min) swap(min, max);
	if (min <= 0.0) throw out_of_range("Minimum must be greater than 0!");

	const double logMin = log(min) / log(base);
	const double logMax = log(max) / log(base);
	const double delta = (logMax - logMin) / float(n - 1);

	Vectord coeff(n); linspace(coeff, logMin, logMax, n);

	Vectord result(n);
	for (auto i = 0; i < n; ++i)
	{
		result.operator[](i) = powf(base, coeff[i]);
	}
	return result;
}

Vectord ParallelFilters::Extensions::loglogspace(double min, double max, int n, double base)
{
	if (base <= 0.0) throw out_of_range("Parameter base must be a positive number!");
	if (n < 0) throw out_of_range("Parameter n must be a positive number!");
	if (max < min) swap(min, max);
	if (min <= 1.0) throw out_of_range("Minimum must be greater than 1!");

	const double logMin = log(log(min) / log(base)) / log(base);
	const double logMax = log(log(max) / log(base)) / log(base);
	const double delta = (logMax - logMin) / float(n - 1);

	Vectord coeff(n); linspace(coeff,logMin, logMax, n);

	vector<double>result(n);
	for (auto i = 0; i < n; ++i)
	{
		result.operator[](i) = pow(base, pow(base, coeff[i]));
	}
	return result;
}


void ParallelFilters::Extensions::hanning(Vectord& input)
{
	auto length = input.size();
	int half = length / 2;
	if (length % 2 == 1)
		half = (length / 2) + 1;
	for (auto i = 0; i < half; i++)
	{
		input[i] = 0.5 * (1.0 - cos(2 * M_PI * (i + 1) / float(length + 1)));
		input[length - i - 1] = 0.5 * (1.0 - cos(2 * M_PI * (i + 1) / float(length + 1)));
	}
}

void ParallelFilters::Extensions::hilbert(Vectord& inputVector, cplxVectord& outputVector)
{
	//here we have 2^N signal, so no zero padding
	const char* error = 0;
	simple_fft::FFT(inputVector, outputVector, inputVector.size(), error);
	
	//Vectord inputImag(inputVector.size());
	auto n = inputVector.size();

	// ---------- Simple FFT library
	
	cplxVectord In(inputVector.size());
	simple_fft::FFT(inputVector, In, inputVector.size(), error);

	for (auto i = 1; i < n / 2; ++i)
		outputVector[i] *= 2.0;
	for (auto i = n / 2 + 1; i < n; ++i)
		outputVector[i] = 0.0;
	
	simple_fft::IFFT(outputVector, outputVector, outputVector.size(), error);

}


cplxVectord::iterator ParallelFilters::Extensions::findCplxPair(cplxVectord::iterator first, cplxVectord::iterator last,
	const std::complex<double> match)
{
	bool matchFound = false;
	double localEps = std::numeric_limits<double>::epsilon()* 100.0 / abs(match);
	//predicate returns true, if difference betew reals and imags are less than eps. Also flips matchFound flag
	auto it = find_if(first, last, [&](const complex<double> &actual)
	{
		if ((std::abs(actual.real() - match.real()) < localEps) && (std::abs(actual.imag() + match.imag()) < localEps))
		{
			matchFound = true; return true;
		}
		else
			return false;
	});
	return matchFound ? it : throw out_of_range("Match was not found!");
}

void ParallelFilters::Extensions::cplxpair(cplxVectord& poles)
{
	vector<int> notMatched;

	sort(poles.begin(), poles.end(), [](complex<double>a, complex<double> b) { return a.real() < b.real(); });

	for (auto i = 0; i < poles.size(); i++)						//sort conjugates together
	{
		try
		{
			auto mypair = findCplxPair(poles.begin() + i + 1, poles.end(), poles[i]);	//throws if pair not found
			//if found, switch with the one after

			if (poles.begin() + i + 1 != mypair)					//if its right next to its pair, do nothing
			{
				iter_swap(mypair, poles.begin() + i + 1);
			}

			if (poles[i].imag() > poles[i + 1].imag())				//swap to negative imag first
				swap(poles[i], poles[i + 1]);
			++i;
		}
		catch (out_of_range& ex)
		{
			notMatched.push_back(i);
		}
	}
	if (!notMatched.empty())											//place the mismatched at the end
	{
		cplxVectord temp;
		for (auto i = 0; i < notMatched.size(); ++i)
		{
			temp.push_back(poles[notMatched[i]]);
		}
		for (auto i = 0; i < notMatched.size(); ++i)
		{
			poles.erase(poles.begin() + (notMatched[i] - i));
		}
		poles.insert(poles.end(), temp.begin(), temp.end());
	}

}

cplxVectord ParallelFilters::Extensions::poly(cplxVectord& roots)
{
	//https://stackoverflow.com/a/42032395
	//sometimes it can happen, that result is negative zero, we will se, if it is a problem later...
	int N = roots.size();
	vector<complex<double>> coefs(N + 1);
	coefs[0] = -roots[0];
	coefs[1] = 1.0;

	for (int k = 2; k <= N; k++)
	{
		coefs[k] = 1.0;
		for (int i = k - 2; i >= 0; i--)
		{
			coefs[i + 1] = coefs[i] - roots[k - 1] * coefs[i + 1];
		}
		coefs[0] *= -roots[k - 1];

		if (remainder(k, 2) == 1)
			coefs[k] = -coefs[k];
	}
	reverse(coefs.begin(), coefs.end());
	return coefs;

}

cplxVectord ParallelFilters::Extensions::roots(cplxVectord& coeffs)
{
	// matrix A - n = 4 (coeffs.size = 5)
	// c1 c2 c3 c4 
	// 1  0  0  0 
	// 0  1  0  0 
	// 0  0  1  0 

	//and each c is divided by c0

	if (coeffs.size() < 2)
		return cplxVectord();
	const auto n = coeffs.size() - 1;
	Eigen::MatrixXcd A(n, n);
	for (auto row = 0; row < A.rows(); ++row)
	{
		for (auto col = 0; col < A.cols(); ++col)
			A(row, col) = complex<double>{ 0.0 };
	}
	for (auto i = 1; i < n; i++)
	{
		A(i, i - 1) = 1.0;
	}
	for (auto i = 0; i < n; i++)
		A(0, i) = -coeffs[i + 1] / coeffs[0];

	auto rts = A.eigenvalues();
	cplxVectord roots(n);
	for (auto i = 0; i < n; i++)
		roots[i] = rts(i);
	return roots;
}

cplxVectord ParallelFilters::Extensions::roots(vector<double>& coeffs)
{
	cplxVectord cfs(coeffs.size());
	for (auto i = 0; i < coeffs.size(); ++i)
		cfs[i] = complex<double>{ coeffs[i] };
	return roots(cfs);
}

Vectord ParallelFilters::Extensions::polyval(Vectord& coeffs, Vectord& values)
{
	Vectord result(values.size());
	for (auto i = 0; i < result.size(); ++i)
	{
		for (auto j = 0; j < coeffs.size(); ++j)
		{
			result[i] += coeffs[j] * pow(values[i], coeffs.size() - j - 1);
		}
	}
	return result;
}

cplxVectord ParallelFilters::Extensions::polyval(cplxVectord& coeffs, cplxVectord& values)
{
	cplxVectord result(values.size());
	for (auto i = 0; i < result.size(); ++i)
	{
		for (auto j = 0; j < coeffs.size(); ++j)
		{
			result[i] += coeffs[j] * pow(values[i], coeffs.size() - j - 1);
		}
	}
	return result;
}

cplxVectord ParallelFilters::Extensions::polyval(Vectord& coeffs, cplxVectord& values)
{
	cplxVectord temp(coeffs.begin(), coeffs.end());
	return polyval(temp, values);
}

cplxVectord ParallelFilters::Extensions::polyval(cplxVectord& coeffs, Vectord& values)
{
	cplxVectord temp(values.begin(), values.end());
	return polyval(coeffs, temp);
}

void ParallelFilters::Extensions::polyval(Vectord& coeffs, Eigen::Map<Eigen::VectorXcd>& input, Eigen::VectorXcd& output)
{
	Eigen::VectorXcd temp(input.size());

	for (auto j = 0; j < coeffs.size(); ++j)
	{
		temp = input.array(). pow(coeffs.size() - j - 1) *coeffs[j];
		output += temp;
	}
}



void ParallelFilters::Extensions::spline(Vectord& X, Vectord& Y, Vectord& points, Vectord& output)
{
	if (X.size() != Y.size()) throw length_error("X and Y length mismatch!");
	tk::spline s;
	s.set_points(X, Y);

	for (auto i = 0; i < output.size(); ++i)
		output[i] = s(points[i]);
}
void ParallelFilters::Extensions::splineNaK(Vectord& X, Vectord& Y, Vectord& points, Vectord& output)
{
	SplineNaK::Spline s;
	s.setPoints(X, Y);
	for (auto i = 0; i < output.size(); ++i)
		output[i] = s(points[i]);
}
void ParallelFilters::Extensions::absSpline(Vectord& X, Vectord& Y, Vectord& points, Vectord &out)
{
	if (X.size() != Y.size()) throw length_error("X and Y length mismatch!");
	tk::spline s;
	s.set_points(X, Y);
	
	for (auto i = 0; i < out.size(); ++i)
		out[i] = abs(s(points[i]));
}
void ParallelFilters::Extensions::absSplineNaK(Vectord& X, Vectord& Y, Vectord& points, Vectord& out)
{
	SplineNaK::Spline s;
	s.setPoints(X, Y);
	for (auto i = 0; i < out.size(); ++i)
		out[i] = abs(s(points[i]));
}

void ParallelFilters::Extensions::conj(cplxVectord2D& matrix)
{
	const complex<double>neg(1., -1.);
	for(auto row =0 ;row <matrix.size();++row)
	{
		for (auto col = 0; col < matrix[0].size(); ++col)
			matrix[row][col] *= neg;
	}
}

#pragma endregion



