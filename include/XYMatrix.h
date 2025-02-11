#ifndef XYMATRIX_H
#define XYMATRIX_H
#include "XYVector.h"

template <class Real>
class CXYMatrix
{
public:
	// construction
    CXYMatrix (int iRows = 0, int iCols = 0);
    CXYMatrix (int iRows, int iCols, const Real* afData);
    CXYMatrix (int iRows, int iCols, const Real** aafEntry);
    CXYMatrix (const CXYMatrix& rkM);
	~CXYMatrix(void);

	// member access
    void SetSize (int iRows, int iCols);
    void GetSize (int& riRows, int& riCols) const;

    int GetRows () const;
    int GetColumns () const;
    int GetQuantity () const;

	operator const Real* () const;
    operator Real* ();
    const Real* operator[] (int iRow) const;
    Real* operator[] (int iRow);
    void SwapRows (int iRow0, int iRow1);
    Real operator( ) (int iRow, int iCol) const;
    Real& operator( ) (int iRow, int iCol);
	void SetRow (int iRow, const CXYVector<Real>& rkV);
    CXYVector<Real> GetRow (int iRow) const;
    void SetColumn (int iCol, const CXYVector<Real>& rkV);
    CXYVector<Real> GetColumn (int iCol) const;

	void SetZeros();

	void SetMatrix (int iRows, int iCols, const Real* afEntry);
    void SetMatrix (int iRows, int iCols, const Real** aafMatrix);

	void GetColumnMajor (Real* afCMajor) const;

	// assignment
    CXYMatrix& operator= (const CXYMatrix& rkM);

	// comparison
    bool operator== (const CXYMatrix& rkM) const;
    bool operator!= (const CXYMatrix& rkM) const;
    bool operator<  (const CXYMatrix& rkM) const;
    bool operator<= (const CXYMatrix& rkM) const;
    bool operator>  (const CXYMatrix& rkM) const;
    bool operator>= (const CXYMatrix& rkM) const;

    // arithmetic operations
    CXYMatrix operator+ (const CXYMatrix& rkM) const;
    CXYMatrix operator- (const CXYMatrix& rkM) const;
    CXYMatrix operator* (const CXYMatrix& rkM) const;
    CXYMatrix operator* (Real fScalar) const;
    CXYMatrix operator/ (Real fScalar) const;
    CXYMatrix operator- () const;

    // arithmetic updates
    CXYMatrix& operator+= (const CXYMatrix& rkM);
    CXYMatrix& operator-= (const CXYMatrix& rkM);
    CXYMatrix& operator*= (Real fScalar);
    CXYMatrix& operator/= (Real fScalar);

	// matrix products
    CXYMatrix Transpose () const;  // M^T
    CXYMatrix TransposeTimes (const CXYMatrix& rkM) const;  // this^T * M
    CXYMatrix TimesTranspose (const CXYMatrix& rkM) const;  // this * M^T

    // matrix-vector operations
    CXYVector<Real> operator* (const CXYVector<Real>& rkV) const;  // M * v
    Real QForm (const CXYVector<Real>& rkU, const CXYVector<Real>& rkV)
        const;  // u^T*M*v

    // Inversion.  The matrix must be square.  The function returns true
    // whenever the matrix is square and invertible.
    bool GetInverse (CXYMatrix<Real>& rkInverse) const;

private:
    // Support for allocation and deallocation.  The allocation call requires
    // m_iRows, m_iCols, and m_iQuantity to have already been correctly
    // initialized.
    void Allocate (bool bSetToZero);
    void Deallocate ();

	// support for comparisons
    int CompareArrays (const CXYMatrix& rkM) const;

	int m_iRows, m_iCols, m_iQuantity;

    // the matrix is stored in row-major form as a 1-dimensional array
    Real* m_afData;

    // An array of pointers to the rows of the matrix.  The separation of
    // row pointers and actual data supports swapping of rows in linear
    // algebraic algorithms such as solving linear systems of equations.
    Real** m_aafEntry;
};

// c * M
template <class Real>
CXYMatrix<Real> operator* (Real fScalar, const CXYMatrix<Real>& rkM);
// v^T * M
template <class Real>
CXYVector<Real> operator* (const CXYVector<Real>& rkV, const CXYMatrix<Real>& rkM);
// x^T * y , x_{1xn}  y_{1xn}
template <class Real>
CXYMatrix<Real> XTTimesY (const CXYVector<Real>& rkV1, const CXYVector<Real>& rkV2);

#include "XYMatrix.inl"

typedef CXYMatrix<float>	CXYMatrixf;
typedef CXYMatrix<double>	CXYMatrixd;

#endif

