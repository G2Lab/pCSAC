#ifndef XYVECTOR_H
#define XYVECTOR_H

#include <string>
#include <cassert>
#include "XYMath.h"

template <class Real>
class CXYVector
{
public:
	// construction
	CXYVector(int iSize = 0);
	CXYVector (int iSize, const Real* afTuple);
	CXYVector (const CXYVector& rkV);
	~CXYVector(void);

    // coordinate access
    void SetSize (int iSize);
    int GetSize () const;
    operator const Real* () const;
    operator Real* ();
    Real operator[] (int i) const;
    Real& operator[] (int i);

	// assignment
    CXYVector& operator= (const CXYVector& rkV);

    // comparison
    bool operator== (const CXYVector& rkV) const;
    bool operator!= (const CXYVector& rkV) const;
    bool operator<  (const CXYVector& rkV) const;
    bool operator<= (const CXYVector& rkV) const;
    bool operator>  (const CXYVector& rkV) const;
    bool operator>= (const CXYVector& rkV) const;

    // arithmetic operations
    CXYVector operator+ (const CXYVector& rkV) const;
    CXYVector operator- (const CXYVector& rkV) const;
    CXYVector operator* (Real fScalar) const;
    CXYVector operator/ (Real fScalar) const;
    CXYVector operator- () const;

	// arithmetic updates
    CXYVector& operator+= (const CXYVector& rkV);
    CXYVector& operator-= (const CXYVector& rkV);
    CXYVector& operator*= (Real fScalar);
    CXYVector& operator/= (Real fScalar);

	// vector operations
    Real Length () const;
    Real SquaredLength () const;
    Real Dot (const CXYVector& rkV) const;
    Real Normalize ();
	Real Sum() const;
	Real Avg() const;

private:
    // support for comparisons
    int CompareArrays (const CXYVector& rkV) const;

	int m_iSize;
    Real* m_afTuple;
};

template <class Real>
CXYVector<Real> operator* (Real fScalar, const CXYVector<Real>& rkV);

#include "XYVector.inl"

typedef CXYVector<float> CXYVectorf;
typedef CXYVector<double> CXYVectord;

#endif
