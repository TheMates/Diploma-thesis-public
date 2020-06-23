/*
* Copyright (c) 2020 Matouš Vrbík
*
* URL	https://github.com/TheMates/Diploma-thesis-public
* 06/2020
* matousvrbik[at]gmail.com
*
* Distributed under MIT licence see LICENCE for details.
*/

using System;
using System.ComponentModel;
using System.IO;
using System.Runtime.InteropServices;


namespace ParFilt
{
    /// <summary>
    /// Path to ParFiltDllWrapper.dll file.
    /// </summary>
    public static class DLLPATH
    {
        public const string PATH = "ParFiltDllWrapper.dll";
    }

    /// <summary>
    /// This is the API for the Parallel Filter Design written in C++. It uses functions from compiles .dll file.
    /// </summary>
    public static class ParFiltDesign
    {
        /// <summary>
        /// Loads reponse from csv file. First column is angular frequency vector, second column is the response magnitude (not in dB).
        /// w and H must be preallocated.
        /// Returns false if loaded data were longer than preallocated, than a new allocation must be done and called again.
        /// </summary>
        /// <param name="w">Angular frequency vector.</param>
        /// <param name="H">Reponse magnitude vector.</param>
        /// <param name="size"></param>
        /// <param name="path"></param>
        /// <returns>True if loading was successfull.</returns>
        [DllImport( ParFilt.DLLPATH.PATH, CallingConvention = CallingConvention.Cdecl)]
        public static extern bool loadResponse(double[] w, double[] H, out int size, string path);


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
        [DllImport(ParFilt.DLLPATH.PATH, CallingConvention = CallingConvention.Cdecl)]
        public static extern void computeResponse(double[] response, double[] w, double[] target, out int size, int nPoles1, int nPoles2, double crossFreq, int crossLength, double lambda1, double lambda2, double sampleRate, int NFIR, bool useNAK);

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
        [DllImport(ParFilt.DLLPATH.PATH, CallingConvention = CallingConvention.Cdecl)]
        public static extern void exportToCsv(double[] w, double[] target, out int size,
            int nPoles1, int nPoles2, double crossFreq, int crossLength, double lambda1, double lambda2,
            double sampleRate, int NFIR, bool useNAK, string fname);

        /// <summary>
        /// Computes the parallel filter coefficients Bm, Am and FIR.
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
        /// <returns>Tuple (Bm,Am,FIR) of coefficients. Bm Am are column vectorw with 2 resp 3 columns, each row is one parallel section.</returns>
        public static (double[,] Bm, double[,] Am, double[] FIR) ComputeCoeffs(double[] w, double[] target, int size, int nPoles1, int nPoles2, double crossFreq, int crossLength, double lambda1, double lambda2, double sampleRate, int NFIR, bool useNAK)
        {
            int nCoeffs = Convert.ToInt32(Math.Ceiling(Convert.ToDecimal(Convert.ToDouble(nPoles1 + nPoles2) / 2.0)));  //wtf mate

            double[] AmCat = new double[2 * nCoeffs];
            double[] BmCat = new double[2 * nCoeffs];
            double[] FIR = new double[NFIR];

            computeCoeffs(AmCat, BmCat, FIR, w, target, size, nPoles1, nPoles2, crossFreq, crossLength, lambda1,
                lambda2, sampleRate, NFIR, useNAK);

            double[,] Am = new double[nCoeffs, 3];
            double[,] Bm = new double[nCoeffs, 2];

            for (var i = 0; i < nCoeffs; ++i)
            {
                Bm[i, 0] = BmCat[i];
                Bm[i, 1] = BmCat[i + nCoeffs];
                Am[i, 0] = 1;
                Am[i, 1] = AmCat[i];
                Am[i, 2] = AmCat[i + nCoeffs];
            }
            return (Bm, Bm, FIR);
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
        [DllImport(ParFilt.DLLPATH.PATH, CallingConvention = CallingConvention.Cdecl)]
        private static extern void computeCoeffs(double[] Am, double[] Bm, double[] FIR, double[] w, double[] target, int size, int nPoles1, int nPoles2, double crossFreq, int crossLength, double lambda1, double lambda2, double sampleRate, int NFIR, bool useNAK);

    }
}

