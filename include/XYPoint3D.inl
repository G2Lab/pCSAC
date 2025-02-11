//----------------------------------------------------------------------------
template <class Real>
CXYPoint3D<Real>::CXYPoint3D(void)
{
}

//----------------------------------------------------------------------------
template <class Real>
CXYPoint3D<Real>::~CXYPoint3D(void)
{
}

//----------------------------------------------------------------------------
template <class Real>
CXYPoint3D<Real>::CXYPoint3D(Real fX, Real fY, Real fZ)
{
	point_[0] = fX;
	point_[1] = fY;
	point_[2] = fZ;
}
//----------------------------------------------------------------------------
template <class Real>
CXYPoint3D<Real>::operator const Real* () const
{
	return point_;
}
//----------------------------------------------------------------------------
template <class Real>
CXYPoint3D<Real>::operator Real* ()
{
	return point_;
}
//----------------------------------------------------------------------------
template <class Real>
Real CXYPoint3D<Real>::operator[] (int i) const
{
    assert(0 <= i && i <= 2);
    return point_[i];
}
//----------------------------------------------------------------------------
template <class Real>
Real& CXYPoint3D<Real>::operator[] (int i)
{
    assert(0 <= i && i <= 2);
    return point_[i];
}
//----------------------------------------------------------------------------
template <class Real>
Real CXYPoint3D<Real>::X () const
{
	return point_[0];
}
//----------------------------------------------------------------------------
template <class Real>
Real& CXYPoint3D<Real>::X ()
{
	return point_[0];
}
//----------------------------------------------------------------------------
template <class Real>
Real CXYPoint3D<Real>::Y () const
{
	return point_[1];
}
//----------------------------------------------------------------------------
template <class Real>
Real& CXYPoint3D<Real>::Y ()
{
	return point_[1];
}
//----------------------------------------------------------------------------
template <class Real>
Real CXYPoint3D<Real>::Z () const
{
	return point_[2];
}
//----------------------------------------------------------------------------
template <class Real>
Real& CXYPoint3D<Real>::Z ()
{
	return point_[2];
}
//----------------------------------------------------------------------------
template <class Real>
CXYPoint3D<Real>& CXYPoint3D<Real>::operator= (const CXYPoint3D& rP)
{
	point_[0] = rP.point_[0];
	point_[1] = rP.point_[1];
	point_[2] = rP.point_[2];
	return *this;

}
//----------------------------------------------------------------------------
template <class Real>
int CXYPoint3D<Real>::CompareArrays (const CXYPoint3D& rP) const
{
	return memcmp(point_, rP.point_, 3*sizeof(Real));
}
//----------------------------------------------------------------------------
template <class Real>
bool CXYPoint3D<Real>::operator== (const CXYPoint3D& rP) const
{
	return CompareArrays(rP) == 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool CXYPoint3D<Real>::operator!= (const CXYPoint3D& rP) const
{
    return CompareArrays(rP) != 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool CXYPoint3D<Real>::operator< (const CXYPoint3D& rP) const
{
    return CompareArrays(rP) < 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool CXYPoint3D<Real>::operator<= (const CXYPoint3D& rP) const
{
    return CompareArrays(rP) <= 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool CXYPoint3D<Real>::operator> (const CXYPoint3D& rP) const
{
    return CompareArrays(rP) > 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool CXYPoint3D<Real>::operator>= (const CXYPoint3D& rP) const
{
    return CompareArrays(rP) >= 0;
}

//----------------------------------------------------------------------------
template <class Real>
CXYPoint3D<Real> CXYPoint3D<Real>::operator+ (const CXYPoint3D& rP) const
{
	return CXYPoint3D(
		point_[0] + rP.point_[0],
		point_[1] + rP.point_[1],
		point_[2] + rP.point_[2]);
}
//----------------------------------------------------------------------------
template <class Real>
CXYPoint3D<Real> CXYPoint3D<Real>::operator- (const CXYPoint3D& rP) const
{
	return CXYPoint3D(
		point_[0] - rP.point_[0],
		point_[1] - rP.point_[1],
		point_[2] - rP.point_[2]);
}
//----------------------------------------------------------------------------
template <class Real>
CXYPoint3D<Real> CXYPoint3D<Real>::operator* (Real fScalar) const
{
	return CXYPoint3D(
		fScalar*point_[0],
		fScalar*point_[1],
		fScalar*point_[2]);
}
//----------------------------------------------------------------------------
template <class Real>
CXYPoint3D<Real> CXYPoint3D<Real>::operator/ (Real fScalar) const
{
	CXYPoint3D kQuot;
	if (fScalar != (Real)0.0){
		Real fInvScalar = ((Real)1.0)/fScalar;
		kQuot.point_[0] = fInvScalar*point_[0];
		kQuot.point_[1] = fInvScalar*point_[1];
		kQuot.point_[2] = fInvScalar*point_[2];
	}else{
		kQuot.point_[0] = CXYMath<Real>::MAX_REAL;
		kQuot.point_[1] = CXYMath<Real>::MAX_REAL;
		kQuot.point_[2] = CXYMath<Real>::MAX_REAL;
	}
	return kQuot;
}
//----------------------------------------------------------------------------
template <class Real>
CXYPoint3D<Real> CXYPoint3D<Real>::operator- () const
{
	return CXYPoint3D(
		-point_[0],
		-point_[1],
		-point_[2]);

}
//----------------------------------------------------------------------------
template <class Real>
CXYPoint3D<Real>& CXYPoint3D<Real>::operator+= (const CXYPoint3D& rP)
{
	point_[0] += rP.point_[0];
	point_[1] += rP.point_[1];
	point_[2] += rP.point_[2];
	return (*this);
}
//----------------------------------------------------------------------------
template <class Real>
CXYPoint3D<Real>& CXYPoint3D<Real>::operator-= (const CXYPoint3D& rP)
{
	point_[0] -= rP.point_[0];
	point_[1] -= rP.point_[1];
	point_[2] -= rP.point_[2];
	return (*this);
}
//----------------------------------------------------------------------------
template <class Real>
CXYPoint3D<Real>& CXYPoint3D<Real>::operator*= (Real fScalar)
{
	point_[0] *= fScalar;
	point_[1] *= fScalar;
	point_[2] *= fScalar;
	return (*this);
}
//----------------------------------------------------------------------------
template <class Real>
CXYPoint3D<Real>& CXYPoint3D<Real>::operator/= (Real fScalar)
{
	if (fScalar != (Real)0.0){
		Real fInvScalar = ((Real)1.0)/fScalar;
		point_[0] *= fInvScalar;
		point_[1] *= fInvScalar;
		point_[2] *= fInvScalar;
	}else{
		point_[0] = CXYMath<Real>::MAX_REAL;
		point_[1] = CXYMath<Real>::MAX_REAL;
		point_[2] = CXYMath<Real>::MAX_REAL;
	}
	return (*this);
}

//----------------------------------------------------------------------------
template <class Real>
Real CXYPoint3D<Real>::Length () const
{
	return CXYMath<Real>::Sqrt(
		point_[0]*point_[0] +
		point_[1]*point_[1] +
		point_[2]*point_[2]);

}
//----------------------------------------------------------------------------
template <class Real>
Real CXYPoint3D<Real>::SquaredLength () const
{
	return 
		(point_[0]*point_[0] +
		point_[1]*point_[1] +
		point_[2]*point_[2]);
}
//----------------------------------------------------------------------------
//template <class Real>
//Real CXYPoint3D<Real>::Distance(const CXYPoint3D<Real>& rP) const
//{
//	return CXYMath<Real>::Sqrt(
//		(point_[0] - rP.point_[0])*(point_[0] - rP.point_[0]) +
//		(point_[1] - rP.point_[1])*(point_[1] - rP.point_[1]) +
//		(point_[2] - rP.point_[2])*(point_[2] - rP.point_[2]));
//}
