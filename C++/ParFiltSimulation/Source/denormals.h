// Macro for killing denormalled numbers
//
// Written by David Oboril, December 2000
//
// Based on IS_DENORMAL macro by Jon Watte
// This code is public domain

#ifndef _denormals_
#define _denormals_

#define undenormalise(sample) if(((*(unsigned int*)&sample)&0x7f800000)==0) sample=0.0f

#endif //_denormals_

//ends
