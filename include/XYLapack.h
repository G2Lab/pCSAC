#ifndef XYLAPACK_H
#define XYLAPACK_H

#include <iostream>
#include "XYMath.h"
#include "XYVector.h"
#include "XYMatrix.h"

using namespace std;

template <class Real>
class CXYLapackEigen{
public:
	// constructor
	CXYLapackEigen(int iSize);
	CXYLapackEigen (const CXYMatrix<Real>& rkM);
	~CXYLapackEigen(void);

	// set the matrix for eigensolving
	Real& operator() (int iRow, int iCol);
	CXYLapackEigen& operator= (const CXYMatrix<Real>& rkM);

	// Get the eigenresults (eigenvectors are columns of eigenmatrix).
	Real GetEigenValue (int i) const;
	//const Real* GetEigenvalues () const;
	CXYVector<Real> GetEigenValues() const;
	CXYVector<Real> GetEigenVector (int i) const;
	const CXYMatrix<Real>& GetEigenVectors () const;

	// set lwork value to optmize
	void SetLwork();
	// solve eigensystem
	void Lapack_SSYEV();


private:
	CXYMatrix<Real> m_kMat;
	char* m_jobz;
	char* m_uplo;
	int m_n;
	Real* m_a;
	int m_lda;
	Real* m_w;
	Real* m_work;
	int m_lwork;
	int m_info;
};

typedef CXYLapackEigen<float> LapackEigenf;
//typedef CXYLapackEigen<double> LapackEigend;

#endif
