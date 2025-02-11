#include "XYLapack.h"

#ifdef WIN32
	extern "C"  int SSYEV(char *jobz, char *uplo, int *n, float *a, 
        int *lda, float *w, float *work, int *lwork, int *info);
#else
	extern "C"  int ssyev_(char *jobz, char *uplo, int *n, float *a, 
        int *lda, float *w, float *work, int *lwork, int *info);
#endif

#ifdef WIN32
	extern "C" int SPOTRI(char *uplo, int *n, float *a, int *lda, int *info);
#else
	extern "C" int spotri_(char *uplo, int *n, float *a, int *lda, int *info);
#endif

//----------------------------------------------------------------------------
template <class Real>
CXYLapackEigen<Real>::CXYLapackEigen (int iSize)
    :
    m_kMat(iSize,iSize)
{
    assert(iSize >= 2);

	m_jobz = (char *)"V";
	m_uplo = (char *)"U";
	m_n = iSize;
	m_a = new Real[m_n*m_n];
	m_lda = max(1,m_n);
	m_w = new Real[m_n];
	m_lwork = -1;
	m_work = new Real[max(1,m_lwork)];
	m_info = 0;
}

//----------------------------------------------------------------------------
template <class Real>
CXYLapackEigen<Real>::CXYLapackEigen (const CXYMatrix<Real>& rkM)
    :
    m_kMat(rkM)
{
    m_n = rkM.GetRows();
    assert(m_n >= 2 && (rkM.GetColumns() == m_n));
	m_jobz = (char *)"V";
	m_uplo = (char *)"U";
	m_a = new Real[m_n*m_n];
	m_kMat.GetColumnMajor(m_a);

	m_lda = max(1,m_n);
	m_w = new Real[m_n];
	m_lwork = -1;
	m_work = new Real[max(1,m_lwork)];
	m_info = 0;
}
//----------------------------------------------------------------------------
template <class Real>
CXYLapackEigen<Real>::~CXYLapackEigen ()
{
	delete [] m_a;
	delete [] m_w;
	delete [] m_work;
}
//----------------------------------------------------------------------------
template <class Real>
Real& CXYLapackEigen<Real>::operator() (int iRow, int iCol)
{
    return m_kMat[iRow][iCol];
}
//----------------------------------------------------------------------------
template <class Real>
CXYLapackEigen<Real>& CXYLapackEigen<Real>::operator= (const CXYMatrix<Real>& rkM)
{
    m_kMat = rkM;
	m_kMat.GetColumnMajor(m_a);
    return *this;
}
////----------------------------------------------------------------------------
template <class Real>
void CXYLapackEigen<Real>::SetLwork()
{
  float *work = new float[max(1,m_lwork)];
  m_lwork = -1;
  m_info = 0;
#ifdef WIN32
  SSYEV(m_jobz, m_uplo, &m_n, m_a, &m_lda, m_w, work, &m_lwork, &m_info);
#else
  ssyev_(m_jobz, m_uplo, &m_n, m_a, &m_lda, m_w, work, &m_lwork, &m_info);
#endif
  m_lwork = int(work[0]);
  delete [] work;

  // re-assign space to m_work
  delete [] m_work;
  m_work = new float[max(1,m_lwork)];
}
//----------------------------------------------------------------------------
template <class Real>
void CXYLapackEigen<Real>::Lapack_SSYEV()
{
  SetLwork();
#ifdef WIN32
  SSYEV(m_jobz, m_uplo, &m_n, m_a, &m_lda, m_w, m_work, &m_lwork, &m_info);
#else
  ssyev_(m_jobz, m_uplo, &m_n, m_a, &m_lda, m_w, m_work, &m_lwork, &m_info);
#endif

  m_kMat = CXYMatrix<float>(m_n, m_n, m_a);
  m_kMat = m_kMat.Transpose();
}

//----------------------------------------------------------------------------

template <class Real>
Real CXYLapackEigen<Real>::GetEigenValue (int i) const
{
	return m_w[i];
}
//----------------------------------------------------------------------------
/*
template <class Real>
const Real* CXYLapackEigen<Real>::GetEigenvalues () const
{
	return m_w;
}
*/
//----------------------------------------------------------------------------
template <class Real>
CXYVector<Real> CXYLapackEigen<Real>::GetEigenValues() const
{
	CXYVector<Real> kVec(m_n,m_w);
	return kVec;
}
//----------------------------------------------------------------------------
template <class Real>
CXYVector<Real> CXYLapackEigen<Real>::GetEigenVector (int i) const
{
    return m_kMat.GetColumn(i);
}
//----------------------------------------------------------------------------
template <class Real>
const CXYMatrix<Real>& CXYLapackEigen<Real>::GetEigenVectors () const
{
    return m_kMat;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// explicit instantiation
//----------------------------------------------------------------------------
template class CXYLapackEigen<float>;
//
//template class CXYLapackEigen<double>;
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// m_info = 0 and m_lwork = -1 to get the value of m_work, this is optimization
// for calculation.
