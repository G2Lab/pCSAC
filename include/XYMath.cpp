#include <limits.h>
#include "XYMath.h"
#include "XYPoint3D.h"


template<> const float CXYMath<float>::ZERO_TOLERANCE = 1e-06f;
template<> const double CXYMath<double>::ZERO_TOLERANCE = 1e-08;

template<> const float CXYMath<float>::k_B = 1.3806505e-23;
template<> const float CXYMath<float>::PI = 3.1415926;

#ifdef WIN32
	template<> const float CXYMath<float>::MAX_REAL = LONG_MAX;
	template<> const double CXYMath<double>::MAX_REAL = _I64_MAX;
#else
	template<> const float CXYMath<float>::MAX_REAL = LONG_MAX;
	template<> const double CXYMath<double>::MAX_REAL = LLONG_MAX;
#endif
//----------------------------------------------------------------------------
template <>
float CXYMath<float>::FastInvSqrt (float fValue)
{
    float fHalf = 0.5f*fValue;
    int i  = *(int*)&fValue;
    i = 0x5f3759df - (i >> 1);
    fValue = *(float*)&i;
    fValue = fValue*(1.5f - fHalf*fValue*fValue);
    return fValue;
}
//----------------------------------------------------------------------------
template <>
double CXYMath<double>::FastInvSqrt (double dValue)
{
    double dHalf = 0.5*dValue;
    Integer64 i  = *(Integer64*)&dValue;
/*
#if defined(WM4_USING_VC70) || defined(WM4_USING_VC6)
    i = 0x5fe6ec85e7de30da - (i >> 1);
#else
    i = 0x5fe6ec85e7de30daLL - (i >> 1);
#endif
*/
    i = 0x5fe6ec85e7de30daLL - (i >> 1);
    dValue = *(double*)&i;
    dValue = dValue*(1.5 - dHalf*dValue*dValue);
    return dValue;
}
//----------------------------------------------------------------------------

// http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm

// Usable AlmostEqual function
//template<>
//bool CXYMath<float>::AlmostEqual2sComplement(float A, float B, int maxUlps)
//{
//	// Make sure maxUlps is non-negative and small enough that the
//	// default NAN won't compare as equal to anything.
//	assert(maxUlps > 0 && maxUlps < 4 * 1024 * 1024);
//	int aInt = *(int*)&A;
//	// Make aInt lexicographically ordered as a twos-complement int
//	if (aInt < 0)
//		aInt = 0x80000000 - aInt;
//	// Make bInt lexicographically ordered as a twos-complement int
//	int bInt = *(int*)&B;
//	if (bInt < 0)
//		bInt = 0x80000000 - bInt;
//	int intDiff = abs(aInt - bInt);
//	if (intDiff <= maxUlps)
//		return true;
//	return false;
//}

// Support functions and conditional compilation directives for the
// master AlmostEqual function.

//inline bool CXYMath<float>::IsInfinite(float A)
//{
//	const kInfAsInt = 0x7F800000;
//	
//	// An infinity has an exponent of 255 (shift left 23 positions) and
//	// a zero mantissa. There are two infinities - positive and negative.
//	if ((*(int*)&A & 0x7FFFFFFF) == kInfAsInt)
//		return true;
//	return false;
//}
//
//inline bool CXYMath<float>::IsNan(float A)
//{
//	// A NAN has an exponent of 255 (shifted left 23 positions) and
//	// a non-zero mantissa.
//	int exp = *(int*)&A & 0x7F800000;
//	int mantissa = *(int*)&A & 0x007FFFFF;
//	if (exp == 0x7F800000 && mantissa != 0)
//		return true;
//	return false;
//}
//
//inline int CXYMath<float>::Sign(float A)
//{
//	// The sign bit of a number is the high bit.
//	return (*(int*)&A) & 0x80000000;
//}
//
//// This is the 'final' version of the AlmostEqualUlps function.
//// The optional checks are included for completeness, but in many
//// cases they are not necessary, or even not desirable.
//template<>
//bool CXYMath<float>::AlmostEqualUlpsFinal(float A, float B, int maxUlps)
//{
//	// There are several optional checks that you can do, depending
//	// on what behavior you want from your floating point comparisons.
//	// These checks should not be necessary and they are included
//	// mainly for completeness.
//	
//#ifdef  INFINITYCHECK
//	// If A or B are infinity (positive or negative) then
//	// only return true if they are exactly equal to each other -
//	// that is, if they are both infinities of the same sign.
//	// This check is only needed if you will be generating
//	// infinities and you don't want them 'close' to numbers
//	// near FLT_MAX.
//	if (IsInfinite(A) || IsInfinite(B))
//		return A == B;
//#endif
//	
//#ifdef  NANCHECK
//	// If A or B are a NAN, return false. NANs are equal to nothing,
//	// not even themselves.
//	// This check is only needed if you will be generating NANs
//	// and you use a maxUlps greater than 4 million or you want to
//	// ensure that a NAN does not equal itself.
//	if (IsNan(A) || IsNan(B))
//		return false;
//#endif
//	
//#ifdef  SIGNCHECK
//	// After adjusting floats so their representations are lexicographically
//	// ordered as twos-complement integers a very small positive number
//	// will compare as 'close' to a very small negative number. If this is
//	// not desireable, and if you are on a platform that supports
//	// subnormals (which is the only place the problem can show up) then
//	// you need this check.
//	// The check for A == B is because zero and negative zero have different
//	// signs but are equal to each other.
//	if (Sign(A) != Sign(B))
//		return A == B;
//#endif
//	
//	int aInt = *(int*)&A;
//	// Make aInt lexicographically ordered as a twos-complement int
//	if (aInt < 0)
//		aInt = 0x80000000 - aInt;
//	// Make bInt lexicographically ordered as a twos-complement int
//	int bInt = *(int*)&B;
//	if (bInt < 0)
//		bInt = 0x80000000 - bInt;
//	
//	// Now we can compare aInt and bInt to find out how far apart A and B
//	// are.
//	int intDiff = abs(aInt - bInt);
//	if (intDiff <= maxUlps)
//		return true;
//	return false;
//}
