#ifndef XYPOINT3D_H
#define XYPOINT3D_H

#include <cassert>
#include "XYMath.h"

template <class Real>
class CXYPoint3D
{
public:
	// construction
	CXYPoint3D(void);
	CXYPoint3D(Real fX, Real fY, Real fZ);

	// destruction
	~CXYPoint3D(void);

	// coordinate access
	operator const Real* () const;
	operator Real* ();
	Real operator[] (int i) const;
	Real& operator[] (int i);
	Real X () const;
	Real& X ();
	Real Y () const;
	Real& Y ();
	Real Z () const;
	Real& Z ();

    // assignment
    CXYPoint3D& operator= (const CXYPoint3D& rP);

	// comparison
    bool operator== (const CXYPoint3D& rP) const;
    bool operator!= (const CXYPoint3D& rP) const;
    bool operator<  (const CXYPoint3D& rP) const;
    bool operator<= (const CXYPoint3D& rP) const;
    bool operator>  (const CXYPoint3D& rP) const;
    bool operator>= (const CXYPoint3D& rP) const;

	// arithmetic operations
    CXYPoint3D operator+ (const CXYPoint3D& rP) const;
    CXYPoint3D operator- (const CXYPoint3D& rP) const;
    CXYPoint3D operator* (Real fScalar) const;
    CXYPoint3D operator/ (Real fScalar) const;
    CXYPoint3D operator- () const;

	// arithmetic updates
	CXYPoint3D& operator+= (const CXYPoint3D& rP);
	CXYPoint3D& operator-= (const CXYPoint3D& rP);
	CXYPoint3D& operator*= (Real fScalar);
	CXYPoint3D& operator/= (Real fScalar);

    // point operations
    Real Length () const;
	Real SquaredLength () const;
	//Real Distance(const CXYPoint3D<Real>& rP) const;

private:
	int CompareArrays (const CXYPoint3D& rP) const;
	Real point_[3];

};

#include "XYPoint3D.inl"

typedef CXYPoint3D<float> CXYPoint3Df;
typedef CXYPoint3D<double> CXYPoint3Dd;
#endif
