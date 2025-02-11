#ifndef XYMATH_H
#define XYMATH_H

#pragma warning (disable: 4305)
#include <cmath>

#define INFINITYCHECK
#define NANCHECK
#define SIGNCHECK


#ifdef WIN32
	typedef __int64 Integer64;
#else 
	#include <stdint.h>
	typedef int64_t Integer64;
#endif
template <class Real>
class CXYMath
{
public:
	// constructor
	CXYMath(void);

	// destructor
	~CXYMath(void);
	
	// traditional algorithm
	static Real Sqrt (Real fValue);
	static Real FAbs (Real fValue);
	
	// fast algorithm
	static Real FastInvSqrt (Real fValue);
	
	// constant
	static const Real k_B; // Boltzmann constant
	static const Real PI; // pi
	static const Real ZERO_TOLERANCE;
	static const Real MAX_REAL;
	
	static bool AlmostEqual2sComplement(float A, float B, int maxUlps);

};
#include "XYMath.inl"
template<> float CXYMath<float>::FastInvSqrt (float fValue);
template<> double CXYMath<double>::FastInvSqrt (double dValue);
#endif
