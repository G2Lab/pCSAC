/*
 * XYEnsemble.h
 * Author: Yun Xu
 * Email: yxu7@uic.edu
 * Date: Sep 20, 2011
 */
#ifndef XYENSEMBLE_H
#define XYENSEMBLE_H

#include "XYMatrix.h"
#include "tree.hh"
#include "XYMath.h"
#include "MersenneTwister.h"
#include "XYSO3Sequence.h"
#include <sstream>
#include "octree.h"
#include "libgen.h"

using namespace std;
using namespace spatialaggregate;

//template<class T1, class T2 >
//struct sort_pair_first_greater {
//	bool operator()(const pair<T1,T2>&left, const pair<T1,T2>&right) {
//		return left.first > right.first;
//	}
//};


class CXYEnsemble{
public:
	CXYEnsemble();
	CXYEnsemble(
		char* cOutPath,
		float fPersistenLength,
		float fCollisionLength,
		float fPackingDensity,
		float fNucleusSphereDiameter,
		int iNumNodes,
		int iNumSamplePoints,
		char* cStartEndFile,
    char* cSegLenFile,
		char* cContIndFile);
	
	~CXYEnsemble(void);
	
	// set and get bindingangle begin
	void SetBendingAngleBeg(float fAng);
	float GetBendingAngleBeg(void);
	
	// set and get bindingangle end
	void SetBendingAngleEnd(float fAng);
	float GetBendingAngleEnd(void);
	
	// set and get number of binding angles
	void SetNumBendingAngle(int iNum);
	int GetNumBendingAngle(void);
	
	// set and get number of torsion angles
	void SetNumTorsionAngle(int iNum);
	int GetNumTorsionAngle(void);
	
	// set and get torsionangle begin
	void SetTorsionAngleBeg(float fAng);
	float GetTorsionAngleBeg(void);
	
	// set and get torsionangle end
	void SetTorsionAngleEnd(float fAng);
	float GetTorsionAngleEnd(void);
	
	// set and get persistence length
	void SetPersistenceLength(float fPL);
	float GetPersistenceLength(void);
	
	// set and get collision length
	void SetCollisionLength(float fCollision);
	float GetCollisionLength(void);
	
	void SetPackingDensity(float fPackingDensity);
	float GetPackingDensity(void);
	
	// set and get number of nodes
	void SetNumNodes(int iNumNodes);
	int GetNumNodes(void);
	
	void SetSamplesOrg();
	CXYMatrix<float> GetSamplesOrg();

	// use V1->V2 as Z axis, reconstruct new xyz coordinate system
	CXYMatrix<float> NewXYZ(CXYVector<float> kV_1, CXYVector<float> kV_2);
	CXYMatrix<float> NewXYZ(CXYVector<float> kV_0, CXYVector<float> kV_1, CXYVector<float> kV_2);
	// get rotation matrix
	CXYMatrix<float> GetRotMatrix(CXYMatrix<float> &rkMXYZ);
	
	// get node samples based on current node
	CXYMatrix<float> GetNodeSamples(tree<CXYVector<float> >::iterator  &itNode, int iSegInd);
	
	
	// initialize chain (0, 0, 0) and (0, 0, PersistenceLength)
	void InitializeChain();
	tree< CXYVector<float> >* GetTree();
	CXYVector<float> RndSetStartPoint(void);
	

	void SetNucleusSphereDiameter(float fNucleusSphereDiameter);
	float GetNucleusSphereDiameter(void);
	
	bool IsInsideSphere(CXYVector<float> kV_point);
	bool IsInsideSphere(float* fCoord);
	bool IsCollision(CXYVector<float> kV_point);
  bool IsCollision(CXYMatrix<float> & MiddleEndPoints);
	// growth chain
	bool GrowOneChain();

	void WriteChain(char* fn);
	void WriteDistance(char* fn);	
  void WriteContactDistance(char* fn);

	void SetOutPath(char* cPathName);
	char* GetOutPath(void);
	// void SetSegLengths(char* cStartEndFile);
	void SetSegLengths(char* cStartEndFile,const char* cMethod);
	void SetSegLengths(void);
	vector<float> & GetSegLengths(void);
	float GetSegLength(int ind);
	
	void SetContIndex(char* cContIndFile);
  void SetContIndex(void);
  
  void GetMiddlePoints(CXYVector<float>& StartPoint, CXYVector<float>& EndPoint, CXYMatrix<float> & MiddleEndPoints);
  bool IsSatisfyCondition(CXYMatrix<float> & MiddleEndPoints);

  void GetGoodPoints( CXYVector<float>& prvnode, CXYMatrix<float>& kMSamplesPoints,vector<int>& GoodPointInd, vector<int>& NoCollisionPointInd);
  
  void WriteWeight(char* fn);
private:
	char* m_cOutPath;
	float m_fBendingAngleBeg; // binding angle begin
	float m_fBendingAngleEnd; // binding angle end
	float m_fTorsionAngleBeg;	// torsion angle begin
	float m_fTorsionAngleEnd;	// torsion angle end
	int   m_iNumBendingAngle; 	// number of binding angle
	int   m_iNumTorsionAngle; 	// number of torsion angle
	float m_fPersistenceLength; // persistence length
	float m_fCollisionLength;   // collision length
	float m_fPackingDensity;   // packing density
	float m_fNucleusSphereDiameter; 
	int   m_iNumNodes;
	int   m_iNumSamplePoints;
	CXYMatrix<float>* m_pMSamplesOrg;
	
	tree< CXYVector<float> >* m_pTChain; // tree stored coordinates
	
	vector<float> m_vfSegLength; // vector of segment length (may different)
	CXYMatrix<float> m_MContInd; // vector of contact index
  
  OcTree<float,int>* m_pOctree;
  Eigen::Matrix< float, 4, 1 > m_center;
  float m_minimumVolumeSize;
  float m_dr;
  int m_maxDepth;
  
  int m_iNumMiddleEndPoints;
  double m_LogWeight;
  vector<float> m_vLogWeight;
};

#endif
