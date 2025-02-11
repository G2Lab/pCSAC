#include <string.h>
//----------------------------------------------------------------------------
template <class Real>
CXYVector<Real>::CXYVector (int iSize)
{
	if (iSize>0)
	{
		m_iSize = iSize;
		m_afTuple = new Real[m_iSize];
		memset(m_afTuple,0,m_iSize*sizeof(Real));
	}
	else
	{
		m_iSize = 0;
		m_afTuple = 0;
	}
}
//----------------------------------------------------------------------------
template <class Real>
CXYVector<Real>::CXYVector (int iSize, const Real* afTuple)
{
	if (iSize >0)
	{
		m_iSize = iSize;
		m_afTuple = new Real[m_iSize];
		size_t uiSize = m_iSize*sizeof(Real);
		memcpy(m_afTuple,afTuple,uiSize);
	}
	else
	{
		m_iSize = 0;
		m_afTuple = 0;
 	}
}
//----------------------------------------------------------------------------
template <class Real>
CXYVector<Real>::CXYVector(const CXYVector& rkV)
{
	m_iSize = rkV.m_iSize;
	if (m_iSize>0)
	{
		m_afTuple = new Real[m_iSize];
		size_t uiSize = m_iSize*sizeof(Real);
		memcpy(m_afTuple,rkV.m_afTuple,uiSize);
	}
	else
	{
		m_afTuple = 0;
	}
}
//----------------------------------------------------------------------------
template <class Real>
CXYVector<Real>::~CXYVector(void)
{
	delete[] m_afTuple;
}
//----------------------------------------------------------------------------
template <class Real>
void CXYVector<Real>::SetSize(int iSize)
{
	delete[] m_afTuple;
	if (iSize > 0)
	{
		m_iSize = iSize;
		m_afTuple = new Real[m_iSize];
		memset(m_afTuple,0,m_iSize*sizeof(Real));
	}
	else
	{
		m_iSize = 0;
		m_afTuple = 0;
	}
}
//----------------------------------------------------------------------------
template <class Real>
int CXYVector<Real>::GetSize() const
{
	return m_iSize;
}
//----------------------------------------------------------------------------
template <class Real>
CXYVector<Real>::operator const Real* () const
{
	return m_afTuple;
}
//----------------------------------------------------------------------------
template <class Real>
CXYVector<Real>::operator Real* ()
{
	return m_afTuple;
}
//----------------------------------------------------------------------------
template <class Real>
Real CXYVector<Real>::operator[] (int i) const
{
    assert(0 <= i && i < m_iSize);
    return m_afTuple[i];
}
//----------------------------------------------------------------------------
template <class Real>
Real& CXYVector<Real>::operator[] (int i)
{
    assert(0 <= i && i < m_iSize);
    return m_afTuple[i];
}
//----------------------------------------------------------------------------
template <class Real>
CXYVector<Real>& CXYVector<Real>::operator= (const CXYVector& rkV)
{
    if (rkV.m_iSize > 0)
    {
        if (m_iSize != rkV.m_iSize)
        {
            delete[] m_afTuple;
            m_iSize = rkV.m_iSize;
            m_afTuple = new Real[m_iSize];
        }
        size_t uiSize = m_iSize*sizeof(Real);
        memcpy(m_afTuple,rkV.m_afTuple,uiSize);
    }
    else
    {
        delete[] m_afTuple;
        m_iSize = 0;
        m_afTuple = 0;
    }
    return *this;
}

//----------------------------------------------------------------------------
template <class Real>
int CXYVector<Real>::CompareArrays (const CXYVector& rkV) const
{
    return memcmp(m_afTuple,rkV.m_afTuple,m_iSize*sizeof(Real));
}
//----------------------------------------------------------------------------
template <class Real>
bool CXYVector<Real>::operator== (const CXYVector& rkV) const
{
    return CompareArrays(rkV) == 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool CXYVector<Real>::operator!= (const CXYVector& rkV) const
{
    return CompareArrays(rkV) != 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool CXYVector<Real>::operator< (const CXYVector& rkV) const
{
    return CompareArrays(rkV) < 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool CXYVector<Real>::operator<= (const CXYVector& rkV) const
{
    return CompareArrays(rkV) <= 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool CXYVector<Real>::operator> (const CXYVector& rkV) const
{
    return CompareArrays(rkV) > 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool CXYVector<Real>::operator>= (const CXYVector& rkV) const
{
    return CompareArrays(rkV) >= 0;
}
//----------------------------------------------------------------------------
template <class Real>
CXYVector<Real> CXYVector<Real>::operator+ (const CXYVector& rkV) const
{
	CXYVector<Real> kSum(m_iSize);
	for (int i = 0; i < m_iSize; i++)
	{
		kSum.m_afTuple[i] = m_afTuple[i] + rkV.m_afTuple[i];
	}
	return kSum;
}

//----------------------------------------------------------------------------
template <class Real>
CXYVector<Real> CXYVector<Real>::operator- (const CXYVector& rkV) const
{
	CXYVector<Real> kDiff(m_iSize);
	for (int i = 0; i < m_iSize; i++)
	{
		kDiff.m_afTuple[i] = m_afTuple[i] - rkV.m_afTuple[i];
	}
	return kDiff;
}

//----------------------------------------------------------------------------
template <class Real>
CXYVector<Real> CXYVector<Real>::operator* (Real fScalar) const
{
    CXYVector<Real> kProd(m_iSize);
    for (int i = 0; i < m_iSize; i++)
    {
        kProd.m_afTuple[i] = fScalar*m_afTuple[i];
    }
    return kProd;
}
//----------------------------------------------------------------------------
template <class Real>
CXYVector<Real> operator* (Real fScalar, const CXYVector<Real>& rkV)
{
    CXYVector<Real> kProd(rkV.GetSize());
    for (int i = 0; i < rkV.GetSize(); i++)
    {
        kProd[i] = fScalar*rkV[i];
    }
    return kProd;
}
//----------------------------------------------------------------------------
template <class Real>
CXYVector<Real> CXYVector<Real>::operator/ (Real fScalar) const
{
	CXYVector<Real> kQuot(m_iSize);
	int i;

	if (fScalar != (Real)0.0)
	{
		Real fInvScalar = ((Real)1.0)/fScalar;
		for (i = 0; i < m_iSize; i++)
		{
			kQuot.m_afTuple[i] = fInvScalar*m_afTuple[i];
		}
	}
	else
	{
		for (i = 0; i < m_iSize; i++)
		{
			kQuot.m_afTuple[i] = CXYMath<Real>::MAX_REAL;
		}
	}
	return kQuot;
}
//----------------------------------------------------------------------------
template <class Real>
CXYVector<Real> CXYVector<Real>::operator -() const
{
	CXYVector<Real> kNeg(m_iSize);
	for (int i = 0; i < m_iSize; i++)
	{
		kNeg.m_afTuple[i] = -m_afTuple[i];
	}
	return kNeg;
}
//----------------------------------------------------------------------------
template <class Real>
CXYVector<Real>& CXYVector<Real>::operator += (const CXYVector& rkV)
{
	for (int i = 0; i < m_iSize; i++)
	{
		m_afTuple[i] += rkV.m_afTuple[i];
	}
	return *this;
}
//----------------------------------------------------------------------------
template <class Real>
CXYVector<Real>& CXYVector<Real>::operator -= (const CXYVector& rkV)
{
	for (int i = 0; i < m_iSize; i++)
	{
		m_afTuple[i] -= rkV.m_afTuple[i];
	}
	return *this;
}
//----------------------------------------------------------------------------
template <class Real>
CXYVector<Real>& CXYVector<Real>::operator *= (Real fScalar)
{
	for (int i = 0; i < m_iSize; i++)
	{
		m_afTuple[i] *= fScalar;
	}
	return *this;
}
//----------------------------------------------------------------------------
template <class Real>
CXYVector<Real>& CXYVector<Real>::operator /= (Real fScalar)
{
	int i;
	if (fScalar != (Real)0.0)
	{
		Real fInvScalar = ((Real)1.0)/fScalar;
		for (i = 0; i < m_iSize; i++)
		{
			m_afTuple[i] *= fInvScalar;
		}
	}
	else
	{
		for (i = 0; i < m_iSize; i++)
		{
			m_afTuple[i] = CXYMath<Real>::MAX_REAL;
		}
	}
}

//----------------------------------------------------------------------------
template <class Real>
Real CXYVector<Real>::Length () const
{
    Real fSqrLen = (Real)0.0;
    for (int i = 0; i < m_iSize; i++)
    {
        fSqrLen += m_afTuple[i]*m_afTuple[i];
    }
    return CXYMath<Real>::Sqrt(fSqrLen);
}
//----------------------------------------------------------------------------
template <class Real>
Real CXYVector<Real>::SquaredLength () const
{
    Real fSqrLen = (Real)0.0;
    for (int i = 0; i < m_iSize; i++)
    {
        fSqrLen += m_afTuple[i]*m_afTuple[i];
    }
    return fSqrLen;
}
//----------------------------------------------------------------------------
template <class Real>
Real CXYVector<Real>::Dot (const CXYVector& rkV) const
{
    Real fDot = (Real)0.0;
    for (int i = 0; i < m_iSize; i++)
    {
        fDot += m_afTuple[i]*rkV.m_afTuple[i];
    }
    return fDot;
}

//----------------------------------------------------------------------------
template <class Real>
Real CXYVector<Real>::Normalize ()
{
    Real fLength = Length();
    int i;

    if (fLength > CXYMath<Real>::ZERO_TOLERANCE)
    {
        Real fInvLength = ((Real)1.0)/fLength;
        for (i = 0; i < m_iSize; i++)
        {
            m_afTuple[i] *= fInvLength;
        }
    }
    else
    {
        fLength = (Real)0.0;
        for (i = 0; i < m_iSize; i++)
        {
            m_afTuple[i] = (Real)0.0;
        }
    }

    return fLength;
}
//----------------------------------------------------------------------------
template <class Real>
Real CXYVector<Real>::Sum () const
{
	Real fSum = Real(0.0);
	for (int i=0; i<m_iSize; i++)
	{
		fSum += m_afTuple[i];
	}
	return fSum;
}
//----------------------------------------------------------------------------
template <class Real>
Real CXYVector<Real>::Avg () const
{
	return Real (Sum()/m_iSize);
}







