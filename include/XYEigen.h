#ifndef XYEIGEN_H
#define XYEIGEN_H

template <class Real>
class CXYEigen{
public:
	// constructor
	CXYEigen(int iSize);
	CXYEigen (const CXYMatrix<Real>& rkM);
	~CXYEigen(void);

    // set the matrix for eigensolving
    Real& operator() (int iRow, int iCol);
    CXYEigen& operator= (const CXYMatrix<Real>& rkM);

	// Get the eigenresults (eigenvectors are columns of eigenmatrix).
    Real GetEigenvalue (int i) const;
    const Real* GetEigenvalues () const;
    CXYVector<Real> GetEigenvector (int i) const;
    const CXYMatrix<Real>& GetEigenvectors () const;

	// solve eigensystem
    void EigenStuffN ();

	// solve eigensystem, use decreasing sort on eigenvalues
    void DecrSortEigenStuffN ();

	// solve eigensystem, use increasing sort on eigenvalues
    void IncrSortEigenStuffN ();


private:
    int m_iSize;
    CXYMatrix<Real> m_kMat;
    Real* m_afDiag;
    Real* m_afSubd;

	// For odd size matrices, the Householder reduction involves an odd
    // number of reflections.  The product of these is a reflection.  The
    // QL algorithm uses rotations for further reductions.  The final
    // orthogonal matrix whose columns are the eigenvectors is a reflection,
    // so its determinant is -1.  For even size matrices, the Householder
    // reduction involves an even number of reflections whose product is a
    // rotation.  The final orthogonal matrix has determinant +1.  Many
    // algorithms that need an eigendecomposition want a rotation matrix.
    // We want to guarantee this is the case, so m_bRotation keeps track of
    // this.  The DecrSort and IncrSort further complicate the issue since
    // they swap columns of the orthogonal matrix, causing the matrix to
    // toggle between rotation and reflection.  The value m_bRotation must
    // be toggled accordingly.
    bool m_bIsRotation;
    void GuaranteeRotation ();

	// Householder reduction to tridiagonal form
    void TridiagonalN ();

    // QL algorithm with implicit shifting, applies to tridiagonal matrices
    bool QLAlgorithm ();

    // sort eigenvalues from largest to smallest
    void DecreasingSort ();

    // sort eigenvalues from smallest to largest
    void IncreasingSort ();

};

typedef CXYEigen<float> Eigenf;
typedef CXYEigen<double> Eigend;

#endif
