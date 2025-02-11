#include <assert.h>
#include <stdlib.h>

//----------------------------------------------------------------------------
template <class Real>
CXYMath<Real>::CXYMath(void)
{
}

//----------------------------------------------------------------------------
template <class Real>
CXYMath<Real>::~CXYMath(void)
{
}
//----------------------------------------------------------------------------
template <class Real>
Real CXYMath<Real>::Sqrt (Real fValue)
{
    return (Real)sqrt((double)fValue);
}
//----------------------------------------------------------------------------
template <class Real>
Real CXYMath<Real>::FAbs (Real fValue)
{
    return (Real)fabs((double)fValue);
}
//----------------------------------------------------------------------------
template<class Real>
bool CXYMath<Real>::AlmostEqual2sComplement(float A, float B, int maxUlps)
{
	// Make sure maxUlps is non-negative and small enough that the
	// default NAN won't compare as equal to anything.
	assert(maxUlps > 0 && maxUlps < 4 * 1024 * 1024);
	int aInt = *(int*)&A;
	// Make aInt lexicographically ordered as a twos-complement int
	if (aInt < 0)
		aInt = 0x80000000 - aInt;
	// Make bInt lexicographically ordered as a twos-complement int
	int bInt = *(int*)&B;
	if (bInt < 0)
		bInt = 0x80000000 - bInt;
	int intDiff = abs(aInt - bInt);
	if (intDiff <= maxUlps)
		return true;
	return false;
}

//----------------------------------------------------------------------------
