
#include <Windows.h>
#include "mex.h"

// Returns local csv separator.
void mexFunction(int nlhs, mxArray* plhs[],
        int nrhs, const mxArray* prhs[]) {

    char sep[2];
	GetLocaleInfoA(LOCALE_USER_DEFAULT, LOCALE_SLIST, sep, sizeof(sep));
    const char * out = &sep[0];
    plhs[0] = mxCreateString(out);    
}