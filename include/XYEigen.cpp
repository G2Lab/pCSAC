#include "XYMath.h"
#include "XYVector.h"
#include "XYMatrix.h"
#include "XYEigen.h"
//----------------------------------------------------------------------------
template <class Real>
CXYEigen<Real>::CXYEigen (int iSize)
    :
    m_kMat(iSize,iSize)
{
    assert(iSize >= 2);
    m_iSize = iSize;
    m_afDiag = new Real[m_iSize];
    m_afSubd = new Real[m_iSize];
    m_bIsRotation = false;
}

//----------------------------------------------------------------------------
template <class Real>
CXYEigen<Real>::CXYEigen (const CXYMatrix<Real>& rkM)
    :
    m_kMat(rkM)
{
    m_iSize = rkM.GetRows();
    assert(m_iSize >= 2 && (rkM.GetColumns() == m_iSize));
    m_afDiag = new Real[m_iSize];
    m_afSubd = new Real[m_iSize];
    m_bIsRotation = false;
}
//----------------------------------------------------------------------------
template <class Real>
CXYEigen<Real>::~CXYEigen ()
{
    delete[] m_afSubd;
    delete[] m_afDiag;
}
//----------------------------------------------------------------------------
template <class Real>
Real& CXYEigen<Real>::operator() (int iRow, int iCol)
{
    return m_kMat[iRow][iCol];
}
//----------------------------------------------------------------------------
template <class Real>
CXYEigen<Real>& CXYEigen<Real>::operator= (const CXYMatrix<Real>& rkM)
{
    m_kMat = rkM;
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
Real CXYEigen<Real>::GetEigenvalue (int i) const
{
    return m_afDiag[i];
}
//----------------------------------------------------------------------------
template <class Real>
const Real* CXYEigen<Real>::GetEigenvalues () const
{
    return m_afDiag;
}
//----------------------------------------------------------------------------
template <class Real>
CXYVector<Real> CXYEigen<Real>::GetEigenvector (int i) const
{
    return m_kMat.GetColumn(i);
}
//----------------------------------------------------------------------------
template <class Real>
const CXYMatrix<Real>& CXYEigen<Real>::GetEigenvectors () const
{
    return m_kMat;
}
//----------------------------------------------------------------------------
template <class Real>
void CXYEigen<Real>::TridiagonalN ()
{
    int i0, i1, i2, i3;

    for (i0 = m_iSize-1, i3 = m_iSize-2; i0 >= 1; i0--, i3--)
    {
        Real fH = (Real)0.0, fScale = (Real)0.0;

        if (i3 > 0)
        {
            for (i2 = 0; i2 <= i3; i2++)
            {
                fScale += CXYMath<Real>::FAbs(m_kMat[i0][i2]);
            }
            if (fScale == (Real)0.0)
            {
                m_afSubd[i0] = m_kMat[i0][i3];
            }
            else
            {
                Real fInvScale = ((Real)1.0)/fScale;
                for (i2 = 0; i2 <= i3; i2++)
                {
                    m_kMat[i0][i2] *= fInvScale;
                    fH += m_kMat[i0][i2]*m_kMat[i0][i2];
                }
                Real fF = m_kMat[i0][i3];
                Real fG = CXYMath<Real>::Sqrt(fH);
                if (fF > (Real)0.0)
                {
                    fG = -fG;
                }
                m_afSubd[i0] = fScale*fG;
                fH -= fF*fG;
                m_kMat[i0][i3] = fF-fG;
                fF = (Real)0.0;
                Real fInvH = ((Real)1.0)/fH;
                for (i1 = 0; i1 <= i3; i1++)
                {
                    m_kMat[i1][i0] = m_kMat[i0][i1]*fInvH;
                    fG = (Real)0.0;
                    for (i2 = 0; i2 <= i1; i2++)
                    {
                        fG += m_kMat[i1][i2]*m_kMat[i0][i2];
                    }
                    for (i2 = i1+1; i2 <= i3; i2++)
                    {
                        fG += m_kMat[i2][i1]*m_kMat[i0][i2];
                    }
                    m_afSubd[i1] = fG*fInvH;
                    fF += m_afSubd[i1]*m_kMat[i0][i1];
                }
                Real fHalfFdivH = ((Real)0.5)*fF*fInvH;
                for (i1 = 0; i1 <= i3; i1++)
                {
                    fF = m_kMat[i0][i1];
                    fG = m_afSubd[i1] - fHalfFdivH*fF;
                    m_afSubd[i1] = fG;
                    for (i2 = 0; i2 <= i1; i2++)
                    {
                        m_kMat[i1][i2] -= fF*m_afSubd[i2] +
                            fG*m_kMat[i0][i2];
                    }
                }
            }
        }
        else
        {
            m_afSubd[i0] = m_kMat[i0][i3];
        }

        m_afDiag[i0] = fH;
    }

    m_afDiag[0] = (Real)0.0;
    m_afSubd[0] = (Real)0.0;
    for (i0 = 0, i3 = -1; i0 <= m_iSize-1; i0++, i3++)
    {
        if (m_afDiag[i0] != (Real)0.0)
        {
            for (i1 = 0; i1 <= i3; i1++)
            {
                Real fSum = (Real)0.0;
                for (i2 = 0; i2 <= i3; i2++)
                {
                    fSum += m_kMat[i0][i2]*m_kMat[i2][i1];
                }
                for (i2 = 0; i2 <= i3; i2++)
                {
                    m_kMat[i2][i1] -= fSum*m_kMat[i2][i0];
                }
            }
        }
        m_afDiag[i0] = m_kMat[i0][i0];
        m_kMat[i0][i0] = (Real)1.0;
        for (i1 = 0; i1 <= i3; i1++)
        {
            m_kMat[i1][i0] = (Real)0.0;
            m_kMat[i0][i1] = (Real)0.0;
        }
    }

    // re-ordering if CXYEigen::QLAlgorithm is used subsequently
    for (i0 = 1, i3 = 0; i0 < m_iSize; i0++, i3++)
    {
        m_afSubd[i3] = m_afSubd[i0];
    }
    m_afSubd[m_iSize-1] = (Real)0.0;

    m_bIsRotation = ((m_iSize % 2) == 0);
}
//----------------------------------------------------------------------------
template <class Real>
bool CXYEigen<Real>::QLAlgorithm ()
{
    const int iMaxIter = 32;

    for (int i0 = 0; i0 < m_iSize; i0++)
    {
        int i1;
        for (i1 = 0; i1 < iMaxIter; i1++)
        {
            int i2;
            for (i2 = i0; i2 <= m_iSize-2; i2++)
            {
                Real fTmp = CXYMath<Real>::FAbs(m_afDiag[i2]) +
                    CXYMath<Real>::FAbs(m_afDiag[i2+1]);
                if ( CXYMath<Real>::FAbs(m_afSubd[i2]) + fTmp == fTmp )
                    break;
            }
            if (i2 == i0)
            {
                break;
            }

            Real fG = (m_afDiag[i0+1] - m_afDiag[i0])/(((Real)2.0) *
                m_afSubd[i0]);
            Real fR = CXYMath<Real>::Sqrt(fG*fG+(Real)1.0);
            if (fG < (Real)0.0)
            {
                fG = m_afDiag[i2]-m_afDiag[i0]+m_afSubd[i0]/(fG-fR);
            }
            else
            {
                fG = m_afDiag[i2]-m_afDiag[i0]+m_afSubd[i0]/(fG+fR);
            }
            Real fSin = (Real)1.0, fCos = (Real)1.0, fP = (Real)0.0;
            for (int i3 = i2-1; i3 >= i0; i3--)
            {
                Real fF = fSin*m_afSubd[i3];
                Real fB = fCos*m_afSubd[i3];
                if (CXYMath<Real>::FAbs(fF) >= CXYMath<Real>::FAbs(fG))
                {
                    fCos = fG/fF;
                    fR = CXYMath<Real>::Sqrt(fCos*fCos+(Real)1.0);
                    m_afSubd[i3+1] = fF*fR;
                    fSin = ((Real)1.0)/fR;
                    fCos *= fSin;
                }
                else
                {
                    fSin = fF/fG;
                    fR = CXYMath<Real>::Sqrt(fSin*fSin+(Real)1.0);
                    m_afSubd[i3+1] = fG*fR;
                    fCos = ((Real)1.0)/fR;
                    fSin *= fCos;
                }
                fG = m_afDiag[i3+1]-fP;
                fR = (m_afDiag[i3]-fG)*fSin+((Real)2.0)*fB*fCos;
                fP = fSin*fR;
                m_afDiag[i3+1] = fG+fP;
                fG = fCos*fR-fB;

                for (int i4 = 0; i4 < m_iSize; i4++)
                {
                    fF = m_kMat[i4][i3+1];
                    m_kMat[i4][i3+1] = fSin*m_kMat[i4][i3]+fCos*fF;
                    m_kMat[i4][i3] = fCos*m_kMat[i4][i3]-fSin*fF;
                }
            }
            m_afDiag[i0] -= fP;
            m_afSubd[i0] = fG;
            m_afSubd[i2] = (Real)0.0;
        }
        if (i1 == iMaxIter)
        {
            return false;
        }
    }

    return true;
}
//----------------------------------------------------------------------------
template <class Real>
void CXYEigen<Real>::DecreasingSort ()
{
    // sort eigenvalues in decreasing order, e[0] >= ... >= e[iSize-1]
    for (int i0 = 0, i1; i0 <= m_iSize-2; i0++)
    {
        // locate maximum eigenvalue
        i1 = i0;
        Real fMax = m_afDiag[i1];
        int i2;
        for (i2 = i0+1; i2 < m_iSize; i2++)
        {
            if (m_afDiag[i2] > fMax)
            {
                i1 = i2;
                fMax = m_afDiag[i1];
            }
        }

        if (i1 != i0)
        {
            // swap eigenvalues
            m_afDiag[i1] = m_afDiag[i0];
            m_afDiag[i0] = fMax;

            // swap eigenvectors
            for (i2 = 0; i2 < m_iSize; i2++)
            {
                Real fTmp = m_kMat[i2][i0];
                m_kMat[i2][i0] = m_kMat[i2][i1];
                m_kMat[i2][i1] = fTmp;
                m_bIsRotation = !m_bIsRotation;
            }
        }
    }
}
//----------------------------------------------------------------------------
template <class Real>
void CXYEigen<Real>::IncreasingSort ()
{
    // sort eigenvalues in increasing order, e[0] <= ... <= e[iSize-1]
    for (int i0 = 0, i1; i0 <= m_iSize-2; i0++)
    {
        // locate minimum eigenvalue
        i1 = i0;
        Real fMin = m_afDiag[i1];
        int i2;
        for (i2 = i0+1; i2 < m_iSize; i2++)
        {
            if (m_afDiag[i2] < fMin)
            {
                i1 = i2;
                fMin = m_afDiag[i1];
            }
        }

        if (i1 != i0)
        {
            // swap eigenvalues
            m_afDiag[i1] = m_afDiag[i0];
            m_afDiag[i0] = fMin;

            // swap eigenvectors
            for (i2 = 0; i2 < m_iSize; i2++)
            {
                Real fTmp = m_kMat[i2][i0];
                m_kMat[i2][i0] = m_kMat[i2][i1];
                m_kMat[i2][i1] = fTmp;
                m_bIsRotation = !m_bIsRotation;
            }
        }
    }
}
//----------------------------------------------------------------------------
template <class Real>
void CXYEigen<Real>::GuaranteeRotation ()
{
    if (!m_bIsRotation)
    {
        // change sign on the first column
        for (int iRow = 0; iRow < m_iSize; iRow++)
        {
            m_kMat[iRow][0] = -m_kMat[iRow][0];
        }
    }
}
//----------------------------------------------------------------------------
template <class Real>
void CXYEigen<Real>::EigenStuffN ()
{
    TridiagonalN();
    QLAlgorithm();
    GuaranteeRotation();
}
//----------------------------------------------------------------------------
template <class Real>
void CXYEigen<Real>::DecrSortEigenStuffN ()
{
    TridiagonalN();
    QLAlgorithm();
    DecreasingSort();
    GuaranteeRotation();
}
//----------------------------------------------------------------------------
template <class Real>
void CXYEigen<Real>::IncrSortEigenStuffN ()
{
    TridiagonalN();
    QLAlgorithm();
    IncreasingSort();
    GuaranteeRotation();
}

//----------------------------------------------------------------------------
// explicit instantiation
//----------------------------------------------------------------------------
template class CXYEigen<float>;
//
template class CXYEigen<double>;
//----------------------------------------------------------------------------
