/*
 * XYSIS.h
 * Author: Yun Xu
 * Email: yxu7@uic.edu
 * Date: Nov 27, 2011
 */
#ifndef XYENSEMBLE_H
#define XYENSEMBLE_H
#include <iostream>
#include <unordered_map>

#include "XYMatrix.h"
#include "tree.hh"
#include "tree_util.hh"
#include "XYMath.h"
#include "MersenneTwister.h"
#include "XYSO3Sequence.h"
#include "XYUtility.h"
#include <functional>
#include <algorithm>
#include <list>
#include <map>
#include <deque>
#include "octree.h"
#include <iterator>


using namespace std;
using namespace spatialaggregate;

typedef pair<float,int> FloatIntPair;
struct FloatIntPairCompare{
  bool operator () ( const FloatIntPair& left, const FloatIntPair& right)
  { return left.first < right.first; }
};

struct RefOrder {
  bool operator ()(pair<int, int> const& left, pair<int, int> const& right) {
    return left.first < right.first;
  }
};

typedef map<
  tree<CXYVector<float> >::iterator, 
  OcTree<float,int>*, 
  tree<CXYVector<float> >::iterator_base_less>  Iterator_OcTree_Map;

typedef map<  
  tree<CXYVector<float> >::iterator, 
  int, 
  tree<CXYVector<float> >::iterator_base_less>  Iterator_Int_Map;

typedef map<
  tree<CXYVector<float> >::iterator, 
  vector<OcTree<float,int>* >,
  tree<CXYVector<float> >::iterator_base_less>  TreeIterator_VecpOcTree_Map;

typedef vector< vector< float> > MultiArray;

typedef map< 
  tree<CXYVector<float> >::iterator,
  MultiArray ,
  tree<CXYVector<float> >::iterator_base_less>  TreeIterator_MultiArray_Map;
  
class CXYSIS{
public:
  CXYSIS();
  CXYSIS(
    char* cOutPath,
    float fPersistenLength,
    float fCollisionLength,
    float fCutoff,    
    float fPackingDensity,
    float fNucleusSphereDiameter,
    int iNumNodes,
    int iNumSamplePoints,
    char* cStartEndFile,
    int iMmax,
    float fRho_1,
    float fRho_2,
    float fRho_3,
    float fTau_t,
    float fAdjust,
    char* cPvalFile,
    char* cContIndFile,
    float fthreshold_A,
    float fthreshold_B
    );
    
  ~CXYSIS(void);
  
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
  bool IsCollision(tree<CXYVector<float> >::iterator  &ritNode, CXYVector<float> kV_point);

  // growth chain
  bool GrowthOneChain();

  void WriteChain(char* fn);
  void WriteDistance(char* fn); 
  void SetOutPath(char* cPathName);
  char* GetOutPath(void);
  void SetSegLengths(void);
//  void SetSegLengths(char* cStartEndFile);
  void SetSegLengths(char* cStartEndFile,const char* cMethod);

  vector<float> & GetSegLengths(void);
  float GetSegLength(int ind);
  
  // for SIS
  void SetMmax(int iMmax);
  int GetMmax(void);

  void SetRho_1(float fRho_1);
  float GetRho_1(void);
  void SetRho_2(float fRho_2);
  float GetRho_2(void);
  void SetRho_3(float fRho_3);
  float GetRho_3(void);
  void SetTau_t(float fTau_t);
  float GetTau_t(void);
  void SetAjust(float fAjust);
  float GetAjust(void);

  void SIS_Algorithm(void);
  void SetDistFileName(char* cFileName);
  char* GetDistFileName(void);
  void SetDistMatrix(void);
  CXYMatrix<float>& GetDistMatrix(void);
  vector<tree<CXYVector<float> >::iterator > GetNodeSet(int iLevel);  
  CXYMatrix<float> GrowthChain(tree< CXYVector<float> >::iterator &ritNode, 
                                CXYVector<float> &rkVPoint);
  CXYMatrix<float> GrowthChain_NoCon(tree< CXYVector<float> >::iterator &ritNode, 
                                     CXYVector<float> &rkVPoint);
  float CalBeta_t(CXYMatrix<float> &kM);
  float BinSearchConstC(vector<float>& vfBeta_t);
  float h1Function(CXYMatrix<float>& kM);
  float h2Function(CXYMatrix<float>& kM);
  int GetCountContactNode(int iNode);
  void GetTargetDistribution(void);
  int GetCountContactNode_All(int iNode);
  void WritePDB(char* cFN, CXYMatrix<float>& kM, int iInd);
  void WritePDBArr(void);
  void GetOneChain(CXYMatrix<float>& rkM,
                   tree<CXYVector<float > >::iterator itNode, int reading_way);
  void InitializeChain_SISRoot();
  void WritePtsArr(void);
  
  map<int,float> GetContactMap(int iNode);
  int GetCountContactNodesFromMap(int ind);
  void RandomPermutation(int ArrN[], int n, int ArrM[], int m); // N > M
  CXYVector<float> GetNodeNextPosition(tree<CXYVector<float> >::iterator  &ritNode, int iSegInd, int ind);
  
  vector<int> GetPrvContactsListFromMap(int ind);
  void GrowOneChain_ByVector(  tree< CXYVector<float> >::iterator &ritCurNode, 
  tree< CXYVector<float> >::iterator &ritSurfNode,  int SegInd, vector<tree<CXYVector<float> >::iterator>* pvCahin);
  void CalBeta_t_ByVector(vector< vector<tree<CXYVector<float> >::iterator >* >& rvvChains,  int iSegInd, vector<float>& rvfBeta_t);
  void h1Function_ByVector(vector<tree<CXYVector<float> >::iterator >& rvOneChain);
  void h2Function_ByVector(vector<tree<CXYVector<float> >::iterator >& rvOneChain, CXYVector<float>& rKV_Dist);
  void h3Function_ByVector(vector<tree<CXYVector<float> >::iterator >& rvOneChain, CXYVector<float>& rKV_Dist);
  CXYVector<float> GetPrvContactsDistFromMap(int ind);
  
  float GetNoConMaxDist(int iSegInd);
  CXYVector<float> GetPrvConPosition(tree<CXYVector<float> >::iterator  &ritNode, int iSegInd);
  
  
  //--------------------------------------------------------
  // Segment may contain many nodes
  //--------------------------------------------------------
  void SetSegSegPval(char* cPvalFile);
  void SetContIndex(char* cContIndFile);
  int FindSegIndFromNodeInd(int iNodeInd);
  int FindContactionTypeFromSegSegInd(int iSegInd1, int iSegInd2);
  vector<int> GetPrvContactSegIndsFromSegInd(int iSegInd);
  int GetPrvContactNumberFromSegInd(int iSegInd);
  void GrowChainSegment_InitRoot(void);
  void GrowChainSegment_SISAlgorithm(void);
  void GrowChainSegment_NoRestriction(int iSegInd, int iNodeInd, int iNodeStartInd, int iNodeEndInd);
  unordered_map<int,CXYVector<float> > GrowChainSegment_GeneratePotentialPoints(void);
  void GrowOneNode_Segment(tree< CXYVector<float> >::iterator &ritCurNode, 
                           tree< CXYVector<float> >::iterator &ritSurfNode,  
                           int SegInd, 
                           vector<tree<CXYVector<float> >::iterator>* pvChain);
  void GrowChainSegment_RandSelect(vector< vector<tree<CXYVector<float> >::iterator >* >& vvChains,
                                   vector< tree< CXYVector<float> >::iterator >& vitPrvNodes,
                                   int iSegInd,
                                   int iNodeInd,
                                   int iNodeStartInd,
                                   int iNodeEndInd);
  void GrowChainSegment_AllSelect(vector< vector<tree<CXYVector<float> >::iterator >* >& vvChains,
                                  vector< tree< CXYVector<float> >::iterator >& vitPrvNodes,
                                  int iSegInd,
                                  int iNodeInd,
                                  int iNodeStartInd,
                                  int iNodeEndInd);
  void UpdateOctreeMap_New(tree<CXYVector<float> >::iterator oldpos, 
                           tree<CXYVector<float> >::iterator newpos,
                           vector<tree<CXYVector<float> >::iterator>& v_octreemap_prvkey);
  void UpdateOctreeMap_Del(vector<tree<CXYVector<float> >::iterator>& v_octreemap_prvkey);
  //  void UpdateOctreeMap_One(
//                           tree<CXYVector<float> >::iterator oldpos, 
//                           tree<CXYVector<float> >::iterator newpos);
//  void UpdateOctreeMap_Batch(void);
  void UpdateMap_TreeIt_Vec_Append(tree<CXYVector<float> >::iterator oldpos, 
                                  tree<CXYVector<float> >::iterator newpos,
                                  int iNodeInd,
                                  int back=0);
  void UpdateMap_TreeIt_Vec_Clean(void); 
//  void UpdateMap_TreeIt_MultiArray_Append(tree<CXYVector<float> >::iterator oldpos, 
//                                          tree<CXYVector<float> >::iterator newpos,
//                                          int iSegInd,
//                                          int iNodeInd,
//                                          int iNodeStartInd,
//                                          int iNodeEndInd);
//  void UpdateMap_TreeIt_MultiArray_Clean();
    
  void FindNodeStartEndIndFromSegInd(int iSegInd, 
                                     int& iNodeStartInd, 
                                     int& iNodeEndInd);
  
  void GrowChainSegment_CalBeta_t(vector< vector<tree<CXYVector<float> >::iterator >* >& vvChains,
                                  vector< tree< CXYVector<float> >::iterator >& vitPrvNodes,
                                  int iSegInd,
                                  int iNodeInd,
                                  vector<FloatIntPair>& vBeta_t);

  void GrowChainSegment_h1Function(vector<tree<CXYVector<float> >::iterator >& vCurChain,
                                   tree< CXYVector<float> >::iterator itPrvNode,
                                   int iSegInd,
                                   int iNodeInd,
                                   int iPrvNodeSegInd);
  void GrowChainSegment_h2Function(vector<tree<CXYVector<float> >::iterator >& vCurChain,
                                   tree< CXYVector<float> >::iterator itPrvNode,
                                   int iSegInd,
                                   int iNodeInd,
                                   int iPrvNodeSegInd);
  void GrowChainSegment_h3Function(vector<tree<CXYVector<float> >::iterator >& vCurChain,
                                   tree< CXYVector<float> >::iterator itPrvNode,
                                   int iSegInd,
                                   int iNodeInd,
                                   int iPrvNodeSegInd,
                                   vector<vector<int > > & vStartEndInd);
    
  float GrowChainSegment_BinSearchConstC(vector<FloatIntPair>& vfBeta_t);
  bool IsCollision(OcTree<float,int>* pOctree_, CXYVector<float>& kV_point);
  bool IsCollision(OcTree<float,int>* pOctree_, float x_, float y_, float z_);
  bool IsOverlap(OcTree<float,int>* pOctree_, CXYVector<float>& kV_point, float dr_);
  bool IsOverlap(OcTree<float,int>* pOctree_, float x_, float y_, float z_, float dr_);
  void FindContactStartEndSegIndFromSegInd(int iSegInd, vector<vector<int> > & vStartEndInd);
  
  void SetContIndex(void);
  void GetMiddlePoints(CXYVector<float>& StartPoint, CXYVector<float>& EndPoint, CXYMatrix<float> & MiddleEndPoints);
  
  void GrowChainSegment_SavePosInd(void);
  void GrowChainSegment_GetErr(void);
  void GrowChainSegment_GetWeight(void);
  void GrowChainSegment_WritePtsArr(tree< CXYVector<float> >::iterator lastNode, int chainInd, int reading_way);
  
  
  //--------------------------------------------------------
  // Map propensity to distance
  //--------------------------------------------------------
  void GrowChainSegment_SISAlgorithm_2(int reading_way);
  void SetSegSegPval_2(char* cPvalFile);
  vector<int> GetPrvContactSegIndsFromSegInd_2(int iSegInd);
  int GetPrvContactNumberFromSegInd_2(int iSegInd);

  void GrowChainSegment_CalBeta_t_2(vector< vector<tree<CXYVector<float> >::iterator >* >& vvChains,
                                  vector< tree< CXYVector<float> >::iterator >& vitPrvNodes,
                                  int iSegInd,
                                  int iNodeInd,
                                  vector<FloatIntPair>& vBeta_t, /* growth function */
                                  vector<FloatIntPair>& vGamma_t); /* target function */ 
  void FindContactStartEndSegIndFromSegInd_2(int iSegInd, vector<vector<int> > & vStartEndInd);
  void GrowChainSegment_h1Function_2(vector<tree<CXYVector<float> >::iterator >& vCurChain,
                                   tree< CXYVector<float> >::iterator itPrvNode,
                                   int iSegInd,
                                   int iNodeInd,
                                   int iPrvNodeSegInd);
  void GrowChainSegment_h2Function_2(vector<tree<CXYVector<float> >::iterator >& vCurChain,
                                   tree< CXYVector<float> >::iterator itPrvNode,
                                   int iSegInd,
                                   int iNodeInd,
                                   int iPrvNodeSegInd,
                                   vector<int> & rvConInd,
                                   vector<float> & rvConDist,
                                   CXYMatrix<float> & rKM_PrvNodePos);
  void GrowChainSegment_h3Function_2(vector<tree<CXYVector<float> >::iterator >& vCurChain,
                                   tree< CXYVector<float> >::iterator itPrvNode,
                                   int iSegInd,
                                   int iNodeInd,
                                   int iPrvNodeSegInd,
                                   vector<vector<int > > & vStartEndInd,
                                   CXYMatrix<float>& rKM_PrvNodePos);
  CXYVector<float> GetPrvContactsDist_2(int iSegInd) ;
  void GrowChainSegment_GetErr_2(void);
  void GrowChainSegment_GetNodes_2(tree<CXYVector<float > >::iterator itNode, CXYMatrix<float>& kM );
  void FindContactStartEndNodeIndFromNodeInd_2(int iNodeInd, vector<vector<int> > & vStartEndInd, int iChain);
  float GrowChainSegment_GetTwoNodesDistFromProfile_2 ( int iLeftNodeInd, int iRightNodeInd, int iChain);
  float GrowChainSegment_GetTotalSegmentLengthFromTwoNodeInd_2 ( int iStartNodeInd, int iEndNodeInd);
  float GrowChainSegment_GetTotalSegmentLengthFromSegInd_2 ( int iStartSegInd, int iEndSegInd);
  vector<int> GetPrvNodeIndsFromNodeInd_2 (int iNodeInd);
  CXYVector<float> GetPrvNodeDistsFromNodeInd_2 (int iNodeInd);
  void GetNodeIndAndDistFromNodeInd_2 (int iNodeInd, vector<int> & rvConInd, vector<float> & rvConDist, int iChain );
  int GrowChainSegment_GetNumMiddlePoints_2 ( int iSegInd );
  void GetMiddlePoints(CXYVector<float>& StartPoint, CXYVector<float>& EndPoint, int iNumMiddlePoints, CXYMatrix<float> & MiddleEndPoints);
  void GrowChainSegment_NormalizeBeta( vector<FloatIntPair>& rvb );
  void GrowChainSegment_WriteContactDistanceArr(void);
  float GrowChainSegment_SumVectorFloatIntPair ( vector<FloatIntPair>& rvFloatIntPair );
  void GrowChainSegment_WriteAllNodesDistanceArr ( void );
  void GrowChainSegment_RearrangeVectorFloatIntPair ( vector<FloatIntPair>& rvFIP_source, vector<FloatIntPair>& rvFIP_target);
  float GrowChainSegment_BinSearchConstC_2 ( vector<FloatIntPair>& vfBeta_t, vector<FloatIntPair>& vfGamma_t );
private:
  char* m_cOutPath;
  float m_fPersistenceLength; // persistence length
  float m_fCollisionLength;   // collision length
  float m_fPackingDensity;   // packing density
  float m_fNucleusSphereDiameter; 
  int   m_iNumSegments;
  int   m_iNumSamplePoints;
  CXYMatrix<float>* m_pMSamplesOrg;
  
  tree< CXYVector<float> >* m_pTChain; // tree stored coordinates
  
  vector<float> m_vfSegLength; // vector of segment length (may different)

  // for SIS
  int m_iMmax; // number of samples selected
  float m_fRho_1, m_fRho_2, m_fRho_3, m_fTau_t; // function parameter
  float m_fAdjust; // adjust gamma matrix
  char* m_cDistFileName;
  CXYMatrix<float>* m_pMDist;
  CXYVector<float>* m_pVErr;
  CXYVector<float>* m_pVEColl;
  vector<CXYTriple<int,int,float> > m_vtriple;

  vector<int> m_vConNodeInd; // unique sorted connection node index
  deque<tree<CXYVector<float> >::iterator> m_qPrvNodes; //stores the last node iterator of each mMax chain

  //--------------------------------------------------------
  // Segment may contain many nodes
  //--------------------------------------------------------
  vector<CXYTriple<int,int,int> > m_vPvalTriple;
  // SegInd, SegInd, Attraction or Repulsion index
  CXYMatrix<int> m_MSegSegPval;
  // SegInd, NodeStartInd, NodeEndInd
  CXYMatrix<int> m_MContInd;
  Iterator_OcTree_Map m_MapOctree;

  // Octree parameter
  Eigen::Matrix<float, 4, 1> m_center;
  float m_minimumVolumeSize;
  float m_dr, m_Dr;
  int m_maxDepth;
//  Iterator_Int_Map m_Map_PrvKey_Count;
  TreeIterator_VecpOcTree_Map m_Treeit_VecpOctree_Map;
  vector<tree<CXYVector<float> >::iterator> m_PrvTreeit_Vec;
//  TreeIterator_MultiArray_Map m_Treeit_MultiArray_Map;
  vector<tree<CXYVector<float> >::iterator> m_PrvSegmentLastTreeit_Vec;
//  int m_iNumMiddleEndPoints;
  
  vector<pair<float,int> >  m_Err_Ind;
  map<tree<CXYVector<float> >::iterator, float,tree<CXYVector<float> >::iterator_base_less> m_Err_map;
  vector<pair<float,int> >  m_Weight_Ind;
  map<tree<CXYVector<float>>::iterator, float, tree<CXYVector<float> >::iterator_base_less> m_Weight_map;
  vector<tree<CXYVector<float> >::iterator> m_posInd_Vec;
  float m_fInvNumber;

  //--------------------------------------------------------
  // Map propensity to distance
  //--------------------------------------------------------
  CXYMatrix<float> m_MSegSegPval_2;
  CXYMatrix<float> m_MSegSegPval_2_reversed;

  //--------------------------------------------------------
  deque<tree<CXYVector<float> >::iterator> m_qPrvNodes_2; //stores the last node iterator of each mMax chain
  tree< CXYVector<float> >* m_pTChain_2; // tree stored coordinates
  map<tree<CXYVector<float> >::iterator, int, tree<CXYVector<float> >::iterator_base_less> m_repeatNode;
  map<tree<CXYVector<float> >::iterator, int, tree<CXYVector<float> >::iterator_base_less> m_reminiscentRepeatNode;
  int m_numInteraction = 0;
  CXYMatrix<float> GetAbsoluteNodeSamples(
    tree<CXYVector<float> >::iterator  &ritNode,
    CXYMatrix<float> kMSamplesOrg);
  void GrowChainSegment_RandSelect_2(
    vector<CXYMatrix<float> >& vkMSamples,
    vector<CXYMatrix<float> >& vkMWeightedSamples,
    vector<int>& vfSampleSelInd,
    vector<float>& vfSampleSelBeta_t,
    vector<float>& vfSampleSelGamma_t,
    float& fSumBeta_t,
    float& fSumGamma_t);
  void GrowChainSegment_Single_RandSelect_2(
    CXYMatrix<float> kMSamples,
    vector<CXYMatrix<float> >& vkMWeightedSamples,
    int iChain,
    tree<CXYVector<float>>::iterator& itPrvNode,
    vector<int>& vfSampleSelInd,
    vector<float>& vfSampleSelBeta_t,
    vector<float>& vfSampleSelGamma_t,
    float& fSumBeta_t,
    float& fSumGamma_t);
  void GrowChainSegment_CalBeta_t_2_Single(
    CXYMatrix<float>& kMSamples,
    tree<CXYVector<float> >::iterator& itPrvNode,
    int iNodeInd,
    vector<FloatIntPair>& vBeta_t,
    vector<FloatIntPair>& vGamma_t,
    int iChain);
  void GrowChainSegment_h1Function_2_mod(
    CXYVector<float>& kV_Vector_cur,
    tree< CXYVector<float> >::iterator itPrvNode,
    int iSegInd);
  void GrowChainSegment_h2Function_2_mod(
    CXYVector<float>& kV_Vector_cur,
    vector<int>& rvConInd,
    vector<float>& rvConDist,
    CXYMatrix<float>& rKM_PrvNodePos);
  void GrowChainSegment_h3Function_2_mod(
    CXYVector<float>& kV_Vector_cur,
    int iSegInd,
    vector<vector<int > >& vStartEndInd,
    CXYMatrix<float>& rKM_PrvNodePos,
    int iChain);
  void AddNewNodes(
    unordered_map<int,
    CXYVector<float> > &potNodes, int reading_way);
  void DistributeInteractions(int n_interactions);
  void DeleteDistribution(int iChainInd);
  void GrowChainSegment_Reset(void);
  float m_fthreshold_A, m_fthreshold_B;
  int m_writtenFiles=0;
  vector<vector<int>> m_MSegSegPval_2_per_chain;
};

#endif
