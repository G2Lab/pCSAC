#include "XYSIS.h"
#include <unordered_map>
#include <algorithm>
//------------------------------------------------------------------------
CXYSIS::CXYSIS()
{
}

CXYSIS::CXYSIS(
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
    float fthreshold_B,
    int interaction_distance,
    float fmc_coeff
    )
{
  //Set all global variables to construtor-given values
  m_cOutPath = new char[1024];
  m_cDistFileName = new char[1024];
  CXYFile::MakeDirectory(cOutPath,0755);
  strcpy(m_cOutPath,cOutPath);
  m_fPersistenceLength = fPersistenLength;
  m_fPackingDensity = fPackingDensity;
  m_fNucleusSphereDiameter = fNucleusSphereDiameter;
  m_iNumSamplePoints = iNumSamplePoints;
  if (strcmp(cStartEndFile, "") == 0){
    // coarse version
    // given number of nodes
    m_iNumSegments = iNumNodes;
    // set each segment length equal to persistence length
    SetSegLengths();
  } else {
    // fine version
    // read start end position file
    // dynamic setting each segment length determined by mass density
    // SetSegLengths(cStartEndFile,"DIF");
    SetSegLengths(cStartEndFile, "LEN");
    m_iNumSegments = m_vfSegLength.size();
  }
  SetContIndex();
  SetCollisionLength( fCollisionLength );
  m_pTChain_2 = new tree< CXYVector<float> > ();
  m_pMSamplesOrg = new CXYMatrix<float>;
  SetSamplesOrg();

  // For SIS
  m_fInvNumber = iMmax * iNumNodes;
  m_iMmax = iMmax;
  m_fRho_1 = fRho_1;
  m_fRho_2 = fRho_2;
  m_fRho_3 = fRho_3;
  m_fTau_t = fTau_t;
  m_fAdjust = fAdjust;
  m_fthreshold_A = fthreshold_A;
  m_fthreshold_B = fthreshold_B;
  m_interaction_distance = interaction_distance;
  mfmc_coeff = fmc_coeff;
  cout<<"INTERACTION DISTANCE = "<<m_interaction_distance<<endl;

  
  // For interaction matrix generation
  SetSegSegPval_2(cPvalFile);

  m_pVErr = new CXYVector<float>(iMmax);
  m_pVEColl = new CXYVector<float>(iMmax);

  // Set octree parameter
  m_center = Eigen::Matrix< float, 4, 1 > (0.0f,0.0f,0.0f,1);
  m_minimumVolumeSize = m_fCollisionLength/4;
  m_dr = m_fCollisionLength*2;
  boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator  =
  boost::make_shared< OcTreeNodeAllocator< float , int > >();
  OcTree<float,int> octree_tmp(m_center,m_minimumVolumeSize,m_fNucleusSphereDiameter,allocator);
  m_maxDepth = (int) ceil(octree_tmp.depthForVolumeSize(m_minimumVolumeSize));
}

//------------------------------------------------------------------------
CXYSIS::~CXYSIS(void)
{
  // Record which octree has been deleted so that no need do again
  vector<size_t> vec_deleted;

  vector<size_t>::iterator it;
  for (TreeIterator_VecpOcTree_Map::iterator map_it=m_Treeit_VecpOctree_Map.begin(); map_it !=m_Treeit_VecpOctree_Map.end(); map_it++) {
    for (vector<OcTree<float,int>* >::iterator vec_it=(map_it->second).begin(); vec_it!=(map_it->second).end(); vec_it++) {
      if (find(vec_deleted.begin(), vec_deleted.end(), (size_t)(*vec_it)) == vec_deleted.end()) {
        vec_deleted.push_back((size_t) *vec_it);
        delete *vec_it;
      }
    }
  }

  delete m_pVErr;
  delete m_pVEColl;
//  delete m_pMDist;
  delete m_cDistFileName;
  delete m_cOutPath;
  delete m_pMSamplesOrg;
  delete m_pTChain_2;

}

// For setting and getting global variables
//------------------------------------------------------------------------
void CXYSIS::SetPersistenceLength(float fPL)
{
  m_fPersistenceLength = fPL;
}
//------------------------------------------------------------------------
float CXYSIS::GetPersistenceLength(void)
{
  return m_fPersistenceLength;
}
//------------------------------------------------------------------------
void CXYSIS::SetCollisionLength(float fCollision)
{
  vector<float>& vfSegment = GetSegLengths();
  float min_diameter = *(min_element(vfSegment.begin(), vfSegment.end()));
  m_fCollisionLength = min(fCollision,min_diameter);
}
//------------------------------------------------------------------------
float CXYSIS::GetCollisionLength(void)
{
  return m_fCollisionLength;
}
//------------------------------------------------------------------------
float CXYSIS::GetPackingDensity(void)
{
  return m_fPackingDensity;
}
//------------------------------------------------------------------------
float CXYSIS::GetNucleusSphereDiameter(void)
{
  return m_fNucleusSphereDiameter;
}
//------------------------------------------------------------------------
int CXYSIS::GetNumNodes(void)
{
  return m_iNumSegments;
}
//------------------------------------------------------------------------
char* CXYSIS::GetOutPath(void)
{
  return m_cOutPath;
}
//------------------------------------------------------------------------
void CXYSIS::SetSegLengths(char* cStartEndFile,const char* cMethod){
  CXYMatrix<float> kMStartEnd = CXYFile::ReadMatrix(cStartEndFile);

  if (strcmp(cMethod, "AVG") == 0) {
    // average different node length
    float len = 0L;
    for (int i= 0; i<kMStartEnd.GetRows(); i++) {
      len = len + ( kMStartEnd[i][1] - kMStartEnd[i][0] +1 ) * GetPackingDensity();
    }
    float avglen = len/(float)kMStartEnd.GetRows();
    // set fixed length
    SetPersistenceLength( avglen );

    for (int i= 0; i<kMStartEnd.GetRows(); i++) {
      m_vfSegLength.push_back(GetPersistenceLength());
    }

  } else if (strcmp(cMethod, "DIF") == 0) {
    // different node length
    for (int i= 0; i<kMStartEnd.GetRows(); i++) {
      int baselen = (int) ( kMStartEnd[i][1] - kMStartEnd[i][0] );
      m_vfSegLength.push_back(((float)baselen * GetPackingDensity()));
    }
  } else if (strcmp(cMethod, "LEN") == 0) {
    for (int i= 0; i<kMStartEnd.GetRows(); i++) {
      m_vfSegLength.push_back(kMStartEnd[i][0]);
    }
  }
  else {
    cerr  << "Please type Method" << endl;
    exit(-1);
  }
}
//------------------------------------------------------------------------
vector<float> & CXYSIS::GetSegLengths(void){
  return m_vfSegLength;
}
//------------------------------------------------------------------------
float CXYSIS::GetSegLength(int ind){
  return m_vfSegLength[ind];
}
//------------------------------------------------------------------------
void CXYSIS::SetContIndex(void)
{
  int iNumNodes = m_iNumSegments + 1;
  m_MContInd.SetSize(iNumNodes,3);
  for (int i=0; i<iNumNodes; i++) {
    m_MContInd[i][0] = i;
    m_MContInd[i][1] = i;
    m_MContInd[i][2] = i;
  }
}
//------------------------------------------------------------------------
void CXYSIS::SetSegLengths(void){
  for (int i= 0; i<GetNumNodes(); i++) {
    m_vfSegLength.push_back(GetPersistenceLength());
  }
}
//------------------------------------------------------------------------

// Generate sphere sample points with radius persistence length.
void CXYSIS::SetSamplesOrg(void){
  CXYSO3Sequence sO3sequence(m_iNumSamplePoints);
  sO3sequence.SetSO3Sequence();
  (*m_pMSamplesOrg) = sO3sequence.GetSO3Sequence();
}
//------------------------------------------------------------------------
CXYMatrix<float> CXYSIS::GetSamplesOrg(){
  return (*m_pMSamplesOrg);
}
//------------------------------------------------------------------------
int CXYSIS::GetMmax(void){
  return m_iMmax;
}
//------------------------------------------------------------------------
float CXYSIS::GetRho_1(void){
  return  m_fRho_1;
}
//------------------------------------------------------------------------
float CXYSIS::GetRho_2(void){
  return m_fRho_2;
}
//------------------------------------------------------------------------
float CXYSIS::GetRho_3(void){
  return m_fRho_3;
}
//------------------------------------------------------------------------
float CXYSIS::GetTau_t(void){
  return m_fTau_t;
}
//------------------------------------------------------------------------
float CXYSIS::GetAjust(void){
  return  m_fAdjust;
}
//------------------------------------------------------------------------
tree< CXYVector<float> >* CXYSIS::GetTree(){
  return m_pTChain_2;
}
//------------------------------------------------------------------------
// Checks if two nodes are colliding or not
bool CXYSIS::IsCollision(tree<CXYVector<float> >::iterator  &ritNode, CXYVector<float> kV_point)
{
  float fCollisionLength = GetCollisionLength();
  // to speed up the caculation, we calculate power(2)
  float fCollisionLength_Squar = fCollisionLength*fCollisionLength;
  // enumerate previous nodes
  tree< CXYVector<float> >* ptr = GetTree();
  //node from where we start comparing the positions/last node added to the chain
  tree<CXYVector<float> >::iterator pos = ritNode;
  pos = ptr->parent(pos);
  while (ptr->is_valid(pos))
  {
    if( (kV_point-(*pos)).SquaredLength() < fCollisionLength_Squar ){
      return true;
    }
    pos = ptr->parent(pos);
  }
  return false;

}
//------------------------------------------------------------------------
float LRInteractionValues(char* cPvalFile){
  CXYMatrix<float> MSegSegPval = CXYFile::ReadMatrix(cPvalFile);
  int iRowCount = MSegSegPval.GetRows();
  int iColCount = MSegSegPval.GetColumns();
  //FLAG: this was updated to 0.95 from 0.97, change it back if it doesn't work
  int cut_off = static_cast<int>(0.95 * iRowCount); //index for the 97 percentile

  std::vector<float> LRValues;

  for (int i=0; i<iRowCount; i++){
    LRValues.push_back(MSegSegPval[i][4]);
  }
  
  std::sort(LRValues.begin(), LRValues.end());
  
  // Guard against out-of-bounds access
  if (cut_off >= LRValues.size()) {
      cut_off = LRValues.size() - 1;
  }
  
  float threshold = LRValues[cut_off]; //threshold value for the interactions that need more details
  cout << "Threshold value = " << threshold << endl;
  return threshold;
}

float superFineStructures(char* cPvalFile){
  CXYMatrix<float> MSegSegPval = CXYFile::ReadMatrix(cPvalFile);
  int iRowCount = MSegSegPval.GetRows();
  int iColCount = MSegSegPval.GetColumns();
  int sf_cut_off = static_cast<int>(0.975 * iRowCount); //index for the 98.5 percentile

  std::vector<float> SFLRValues;

  for (int i=0; i<iRowCount; i++){
    SFLRValues.push_back(MSegSegPval[i][4]);
  }
  
  std::sort(SFLRValues.begin(), SFLRValues.end());
  
  // Guard against out-of-bounds access
  if (sf_cut_off >= SFLRValues.size()) {
      sf_cut_off = SFLRValues.size() - 1;
  }
  
  float sf_threshold = SFLRValues[sf_cut_off]; //threshold value for the interactions that need more details
  cout << "Threshold value for finest structures = " << sf_threshold << endl;
  return sf_threshold;
}

/*
 * ===  FUNCTION  ======================================================================
 *  Name:  SetSegSegPval_2
 *  Description: Fills m_MSegSegPval_2 matrix with [node i, node j, number of interactions to distribute] using 
 *  input MSegSegPval interaction frequencies 
 * =====================================================================================
 */
void CXYSIS::SetSegSegPval_2(char* cPvalFile){
  CXYMatrix<float> MSegSegPval = CXYFile::ReadMatrix(cPvalFile);
  int iRowCount = MSegSegPval.GetRows();
  int iColCount = MSegSegPval.GetColumns();
  m_nb_nodes = MSegSegPval[iRowCount-1][3];

  m_MSegSegPval_2.SetSize(iRowCount*iRowCount, MSegSegPval.GetColumns());
  int i_count=0;

  float threshold = LRInteractionValues(cPvalFile); 
  int fine=0;

  float tsuperFine = superFineStructures(cPvalFile);
  int superFine=0;

  MTRand drand;
  MTRand xrand;

  for (int i=0; i<iRowCount; i++) {
      int start1 = MSegSegPval[i][0];
      int end1 = MSegSegPval[i][1];
      int start2 = MSegSegPval[i][2];
      int end2 = MSegSegPval[i][3];

    if(MSegSegPval[i][4] >= threshold && abs(start2 - end1) != 1){ //if the interaction is important enough to be considered in the fine version. 
      fine++;

      for (int f1 = start1; f1 <= end1; f1++) {
        for (int f2 = start2; f2 <= end2; f2++) {

          if(abs(f1-f2) == 1) continue;

          m_MSegSegPval_2[i_count][0] = f1;
          m_MSegSegPval_2[i_count][1] = f2;
          
          if (MSegSegPval[i][4] >= tsuperFine) {
            if (round(MSegSegPval[i][4] * GetMmax() + (GetMmax()/50)) > GetMmax()) {
              //cout << "ERROR : There are more interactions than maximum number of chains" << endl;
              m_MSegSegPval_2[i_count][2] = round(MSegSegPval[i][4] * GetMmax());
            }
            else {
              m_MSegSegPval_2[i_count][2] = round(MSegSegPval[i][4] * GetMmax() + (GetMmax()/50));
              superFine++;
            }
          }
          else {
            m_MSegSegPval_2[i_count][2] = round(MSegSegPval[i][4] * GetMmax());
          }
          cout << m_MSegSegPval_2[i_count][0] << "\t" << m_MSegSegPval_2[i_count][1] << "\t" << m_MSegSegPval_2[i_count][2] << endl;
          i_count++;
        }
      }
    } else {

        if(abs(start1-start2) <= 4) continue;
        
        int cc = MSegSegPval[i][1] - MSegSegPval[i][0];
        int vv = drand.randInt(cc);
        int f1 = MSegSegPval[i][0] + vv;
        int cc2 = MSegSegPval[i][3] - MSegSegPval[i][2];
        int vv2 = xrand.randInt(cc2);
        int f2 = MSegSegPval[i][2] + vv2;

        if(abs(f1-f2) == 0) continue;

        m_MSegSegPval_2[i_count][0] = f1;
        m_MSegSegPval_2[i_count][1] = f2;
        m_MSegSegPval_2[i_count][2] = round(MSegSegPval[i][4] * GetMmax());
        cout << m_MSegSegPval_2[i_count][0] << "\t" << m_MSegSegPval_2[i_count][1] << "\t" << m_MSegSegPval_2[i_count][2] << endl;
        i_count++;
      
    }
  }
  m_numInteraction = i_count;
  cout<<"Number of interactions to satisfy : "<<m_numInteraction<<endl;
  cout<<"Number of fine interactions : "<<fine-superFine<<endl;
  cout<<"Number of super fine interactions : "<<superFine<<endl;
}

//------------------------------------------------------------------------
 /*
 * ===  FUNCTION  ======================================================================
 *  Name:  DistributeInteractions
 *  Description: Generates the different interaction distributions that altogether 
 *  respect the frequency given in input
 * =====================================================================================
 */
void CXYSIS::DistributeInteractions(int iChainInd){
  m_MSegSegPval_2_distributed.clear();
  necessary_distribution.clear();
  int RowNumber = m_MSegSegPval_2.GetRows();
  int m = GetMmax();
  int n_interactions;
  int rint;
  for(int iRow=0; iRow<RowNumber; iRow++){
    n_interactions = m_MSegSegPval_2[iRow][2];
    if(0<n_interactions && n_interactions<(m-iChainInd)){
      //Monte-carlo factor
      float fr = rand() / static_cast<float>(RAND_MAX); //random float between 0 and 1
      if(fr<mfmc_coeff){
        m_MSegSegPval_2_distributed.push_back(m_fCollisionLength); 
      }
      else{
        //choose randomly if there will be an interaction
        rint = (rand()%(m)); //1) rint = (rand()%(m)) & rint<n_interactions 2) rint = (rand()%(m-iChainInd)) & rint<n_interactions 3) rint=(rand()%2) & rint==1
        m_MSegSegPval_2_distributed.push_back((rint<n_interactions?m_fCollisionLength:m_fCollisionLength*5)); // setting up alwas as m_fCollisionLength
        // but m_fCollisionLength is alwas half of the distance between two nodes ... why would it be 1?
      }
      necessary_distribution.push_back(0); //the interaction was chosen randomly and not forced
    }
    else{
      //when the interaction has to happend or when this interaction already happened enough times
      m_MSegSegPval_2_distributed.push_back(n_interactions<=0?m_fCollisionLength*5:m_fCollisionLength);
      necessary_distribution.push_back(1);
    }
  }
  //"returns" m_MSegSegPval_2_distributed with new distribution of interaction to follow
}
//------------------------------------------------------------------------
 /*
 * ===  FUNCTION  ======================================================================
 *  Name:  ClearDistribution
 *  Description: Resets every container to start a new chain
 * =====================================================================================
 */
void CXYSIS::ClearDistribution(){
    m_MSegSegPval_2_distributed.clear();
    tree< CXYVector<float> >* ptr = GetTree();
    ptr->clear();
    m_qPrvNodes_2.clear();
    m_repeatNode.clear();
    m_reminiscentRepeatNode.clear();
    m_Treeit_VecpOctree_Map.clear();
    m_PrvTreeit_Vec.clear();

    //clears every global variable that relates to the latest chain
}
//------------------------------------------------------------------------
/*void CXYSIS::ConfirmDistribution(){
    tree<CXYVector<float>>::iterator itNode;
    tree< CXYVector<float> >* ptr = GetTree();
    int iRowNum = ptr->max_depth()+1;
    tree<CXYVector<float>>::iterator lastItNode = ptr->end();
    CXYMatrix<float> kM_PrvNodePos(iRowNum, m_nb_features);
    GrowChainSegment_GetNodes_2(lastItNode, kM_PrvNodePos);
    int index_node_A, index_node_B;
    float fPoint_A[3], fPoint_B[3], fD_AB;
    int erase;


    for(int i=0; i<m_numInteraction; i++){
      index_node_A = m_MSegSegPval_2[i][0];
      fPoint_A = {kM_PrvNodePos[index_node_A][0],kM_PrvNodePos[index_node_A][1],kM_PrvNodePos[index_node_A][2]};
      CXYVector<float> kV_Point_A(3,fPoint_A);

      index_node_B = m_MSegSegPval_2[i][1];
      fPoint_B = {kM_PrvNodePos[index_node_B][0],kM_PrvNodePos[index_node_B][1],kM_PrvNodePos[index_node_B][2]};
      CXYVector<float> kV_Point_B(3,fPoint_B);

      CXYVector<float> kV_Point_diff = kV_Point_A - kV_Point_B;
      fD_AB = kVDiff.Length();

      erase = fD_AB<800/3?1:0; //mon probleme c'est que des fois meme quand c'est 0 il faut l'enlever
      m_MSegSegPval_2[i][2]=min(0, m_MSegSegPval_2[i][2]-erase)

    }
}*/
//------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *  Name:  GrowChainSegment_ConfirmDistribution_2
 *  Description: Minimizes, when possible, the error between the target distribution and
 *  the final one
 * =====================================================================================
 */
void CXYSIS::GrowChainSegment_ConfirmDistribution_2(tree<CXYVector<float>>::iterator lastNode){
  vector<int> new_dis; //closest distribution possible to the generated
  bool current_state, target;
  int index_node_A, index_node_B, count;
  float fD_AB;
  tree< CXYVector<float> >* ptr = GetTree();
  int iRowNum = ptr->max_depth()+1;
  tree<CXYVector<float>>::iterator itNode;
  CXYMatrix<float> kM_PrvNodePos(iRowNum, m_nb_features);
  GrowChainSegment_GetNodes_2(lastNode, kM_PrvNodePos);
  float fBeta_t, h1, h2;
  int mindex;
  cout<<"Errors in original distribution : "<<(*lastNode)[8]<<endl;

  /*for(int i=0; i<m_numInteraction; i++){
      cout<<"Before:"<<"\t"<<m_MSegSegPval_2[i][0]<<"\t"<<m_MSegSegPval_2[i][1]<<"\t"<<m_MSegSegPval_2[i][2]<<endl;
  }*/


  for(int i=0; i<m_numInteraction; i++){
    index_node_A = m_MSegSegPval_2[i][0];
    float fPoint_A[3] = {kM_PrvNodePos[index_node_A][0],kM_PrvNodePos[index_node_A][1],kM_PrvNodePos[index_node_A][2]};
    CXYVector<float> kV_Point_A(3,fPoint_A);
    index_node_B = m_MSegSegPval_2[i][1];
    float fPoint_B[3] = {kM_PrvNodePos[index_node_B][0],kM_PrvNodePos[index_node_B][1],kM_PrvNodePos[index_node_B][2]};
    CXYVector<float> kV_Point_B(3,fPoint_B);
    CXYVector<float> kV_Point_diff = kV_Point_A - kV_Point_B;
    fD_AB = kV_Point_diff.Length();

    
    current_state = fD_AB<=(m_interaction_distance+1)?1:0;
    target = m_MSegSegPval_2_distributed[i]==m_fCollisionLength?1:0; // if collision length is equal to m_MSegSegPval_2_distributed[i] then target is 1


    //cout<<target<<" ; "<<current_state<<" ; "<<necessary_distribution[i]<<" ; "<<(current_state!=target)<<" ; "<<(necessary_distribution[i]==0)<<" ; "<<m_MSegSegPval_2[i][2]<<" ; ";

   // cout<<"before"<<m_MSegSegPval_2[i][0]<<"\t"<<m_MSegSegPval_2[i][1]<<"\t"<<m_MSegSegPval_2[i][2]<<endl;  
  
    if(current_state==target || necessary_distribution[i]==1){
      m_MSegSegPval_2[i][2]=m_MSegSegPval_2[i][2]-target;
    }
    else{
      count++;
      m_MSegSegPval_2[i][2]=m_MSegSegPval_2[i][2]-current_state;
      mindex = max(index_node_A, index_node_B);
      h2 = --kM_PrvNodePos[mindex][4];
      h1 = kM_PrvNodePos[mindex][3];
      //fBeta_t  = exp( - (m_fRho_1 * h1 + 2 * m_fRho_2 * h2)/ (m_fTau_t) );
      kM_PrvNodePos[mindex][7] = - (m_fRho_1 * h1 + m_fRho_2 * h2)/ (m_fTau_t);
      itNode = (ptr->begin_fixed(ptr->begin(), mindex));
      *itNode = kM_PrvNodePos.GetRow(mindex);
    }
    //cout<<m_MSegSegPval_2[i][2]<<endl;
  }
  //cout<<endl;
  cout<<"Number of corrections : "<<count<<endl;
  
  /*for(int i=0; i<m_numInteraction; i++){
      cout<<"After:"<<"\t"<<m_MSegSegPval_2[i][0]<<"\t"<<m_MSegSegPval_2[i][1]<<"\t"<<m_MSegSegPval_2[i][2]<<endl;
  }*/
 

 /*for(int i=0; i<m_numInteraction; i++){
    cout<<m_MSegSegPval_2[i][2]<<" ; ";
  }
  cout<<endl;*/
  float fh2t = kM_PrvNodePos.GetColumn(4).Sum();
  float wt = kM_PrvNodePos.GetColumn(7).Sum();
  (*(lastNode))[7] = wt-1; //total weight after modifications
  (*(lastNode))[8] = fh2t; //total h2 errors after modifications
}
//------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *  Name:  GrowChainSegment_SISAlgorithm_2
 *  Description:  Main algorithm for chain generation :
 * - Initialize root
 * - Choose new candidate node based on weight calculation and adds it to the chain (tree
 *  and m_qPrvNodes_2) if error doesn't go above a certain threshold
 * - Tries again or generates another target distribution if the current chain failed too 
 *  many times
 * =====================================================================================
 */
void CXYSIS::GrowChainSegment_SISAlgorithm_2(void){
    unordered_map<int,CXYVector<float> > potNodes;
    tree< CXYVector<float> >* ptr = GetTree();

    for(int iChain=0; iChain<m_iMmax; iChain++){
      cout<<"Starting chain "<<iChain<<endl;
      // Initialization of first segment
      StartNewChain(iChain);
      m_iteration = 0;
      while(!m_qPrvNodes_2.empty()){
          m_iteration++;
          //cout<<ptr->depth(m_qPrvNodes_2[0])+1<<" : ";
          potNodes = GrowChainSegment_GeneratePotentialPoints(); //optimized random potential new node at the end of each parallel chain
          AddNewNodes(potNodes); //evalue the candidates and add the valid nodes
          if(m_nb_full_reset>=10){
            cout<<"Changing distribution for chain "<<iChain<<endl;
            StartNewChain(iChain);
            }
      }
    }
    

  //---------------
  // output
  float err=0;
  float iterations=0;
  for(map<tree<CXYVector<float>>::iterator, float,tree<CXYVector<float> >::iterator_base_less>::iterator it = m_Err_map.begin();  it!=m_Err_map.end(); it++){
    cout << "# Err= " << (*it).second << " # Weight= "<<m_Weight_map[(*it).first]<< endl;
    err+=(*it).second;
    iterations+=m_Iteration_map[(*it).first];
  }
  cout << "AVERAGE ERROR "<<err/GetMmax() << endl;
  cout << "AVERAGE ITERATIONS "<< iterations/GetMmax() << endl;

}
//------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *  Name:  StartNewChain
 *  Description:  Resets everything, creates a new target distribution and initialize the
 *  current chain
 * =====================================================================================
 */
void CXYSIS::StartNewChain(int iChain){
  ClearDistribution();
  DistributeInteractions(iChain);
  GrowChainSegment_InitRoot();
  m_nb_full_reset=0;
}
//------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *  Name:  GrowChainSegment_InitRoot
 *  Description:  Creates the first nodes of the chain in a two layer tree of nodes 
 *  ((0,0,0) and (0,0,300)) and updates the VecOctreeMap
 * =====================================================================================
 */
void CXYSIS::GrowChainSegment_InitRoot(void){
  // clear the queue
  m_qPrvNodes_2.clear();

  tree< CXYVector<float> >* ptr = GetTree();
  tree< CXYVector<float> >::iterator oldpos, newpos;

  // first node
  float fV0[9] = {0,0,0,0,0,0,0,1,0};  // x, y, z, h1, h2, h3, total target interactions, w, h2 cumulated
  oldpos = ptr->insert(ptr->begin(), CXYVector<float>(m_nb_features,fV0));
  // second node
  float fLen = GetSegLength(0);
  float fV1[9] = {0,0,fLen,0,0,0,0,1,0};
  newpos = ptr->append_child(oldpos, CXYVector<float>(m_nb_features,fV1));

  m_qPrvNodes_2.clear();
  m_qPrvNodes_2.push_back(newpos);

  // Only one segment, so need update one time
  // put last point in queue
  m_PrvTreeit_Vec.clear();
  m_PrvTreeit_Vec.push_back(oldpos);

  UpdateMap_TreeIt_Vec_Append(oldpos, newpos, 0, 0);
  UpdateMap_TreeIt_Vec_Clean();

  m_PrvTreeit_Vec.push_back(newpos);

}
//------------------------------------------------------------------------
/* GrowChainSegment_GeneratePotentialPoints doit :
- comprendre l'avancee des chaines avec depth mqprvnodes2
- generer les 640 points spheriques dans un vecteur                                   //(function)
- pour i in range size mqprvnodes2
    si la chaine n'est pas deja terminee :
    - 640 pts + coordonnees de mqprvnodes2[i]
    - while is collision :                                                            //(function)
        - choose random                                                        //(pq function ici?)
    else :
    - choose -1,-1,-1 SI CHAINE TERMINEE
- return potnodes de la taille de mqprvnodes2 sans toucher a mqprvnodes2*/
/*
 * ===  FUNCTION  ======================================================================
 *  Name:  GrowChainSegment_GeneratePotentialPoints
 *  Description:  Returns a map associating the last iterator of each chain and the next randomly selected
 *  sample point (now that it works seauentially and not simultaneously, having a map is not useful
 *  but I just didn't change the function to be more flexible)
 * =====================================================================================
 */
unordered_map<int,CXYVector<float> > CXYSIS::GrowChainSegment_GeneratePotentialPoints(void){
    tree<CXYVector<float> >::iterator itPrvNode; 
    CXYMatrix<float> kMSamplesOrg = GetSamplesOrg(); //only relative normalized x y z
    vector<CXYMatrix<float> > vkMSamples; //vector with matrix of potential points for each chain
    CXYVector<float> kVector(m_nb_features); //vector with selected potential point
    tree<CXYVector<float> >* ptr = GetTree();
    vector<int> vfSampleSelInd; //vector with selected indexes in vkMSamples
    vector<CXYMatrix<float> > vkMWeightedSamples; //vector with updated matrices with the weight calculated
    vector<float> vfSampleSelBeta_t; //vector with all the Beta each selected candidate of each chain
    vector<float> vfSampleSelGamma_t; //vector with all the Gamma each selected candidate of each chain
    int iChainInd, iSampleInd, iNodeInd;
    float fSumBeta_t = 0.f, fSumGamma_t = 0.f;
    unordered_map<int,CXYVector<float> > randNodes; //keeping the randomly selected node for the corresponding iterator in mptchain2

    //for each parallel chain get the sphere of sample nodes
    for(int i=0; i<m_qPrvNodes_2.size(); i++){ 
        itPrvNode = m_qPrvNodes_2[i];
        iNodeInd = ptr->depth(itPrvNode);
        CXYMatrix<float> kMSamples = GetAbsoluteNodeSamples(itPrvNode,kMSamplesOrg * GetSegLength(iNodeInd)); 
        vkMSamples.push_back(kMSamples);
    }

    //select randomly potential points for each chain and puts their indexes in vfSampleSelInd
    GrowChainSegment_RandSelect_2(vkMSamples, vkMWeightedSamples, vfSampleSelInd, vfSampleSelBeta_t, vfSampleSelGamma_t, fSumBeta_t, fSumGamma_t);

    int im_t = vfSampleSelInd.size();
    //add the selected node in the final map if the chain is not finished
    for (int iSample = 0; iSample < im_t; iSample++) {
        iChainInd = min(int(m_qPrvNodes_2.size()-1), iSample); //vfSampleSelInd[iSample]/numSamples; //what chain in mqprvnode
        iSampleInd = vfSampleSelInd[iSample]; //%numSamples; //what sample in corresponding kMSamples
        itPrvNode = m_qPrvNodes_2[iChainInd];
        iNodeInd = ptr->depth(itPrvNode);
        if(iNodeInd==m_iNumSegments){
            float endOfChain[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1}; //indicates that the chain has actually ended and the selected node is useless
            kVector = CXYVector<float>(m_nb_features, endOfChain);
        }
        else{
            kVector = vkMWeightedSamples[iChainInd].GetRow(iSampleInd);
        }
        randNodes.insert(std::make_pair<int,CXYVector<float> >(int(iSample), CXYVector<float>(kVector)));
    }
    return randNodes;
}
//------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *  Name:  GetAbsoluteNodeSamples
 *  Description:  Adds spherical relative node samples to previous node and returns them 
 *  with previous weights and h functions (n, m_nb_features)
 * =====================================================================================
 */
CXYMatrix<float> CXYSIS::GetAbsoluteNodeSamples(tree<CXYVector<float> >::iterator  &ritNode, CXYMatrix<float> kMSamplesOrg)
{
  tree<CXYVector<float> >* ptr = GetTree();
  tree<CXYVector<float> >::iterator pos = ritNode;
  //tree<CXYVector<float> >::iterator root = ptr->begin();

  // X, Y, Z, h1, h2, h3 of previous node
  float XYZ_n2[3] = {(*pos)[0], (*pos)[1], (*pos)[2]};
  float fh1 = (*pos)[3]; //records the number of collisions for this node
  float fh2 = (*pos)[4]; // records the missed interactions at this node
  float fh3 = (*pos)[5]; // records the triangular inequality error for this node
  float interaction = (*pos)[6]; // records total interactions until this node
  float fWeight = (*pos)[7]; // records the exponential weight that keeps in track collision and interaction
  float fh2c = (*pos)[8]; // records total h2 errors


  CXYVector<float> kV_2(3,XYZ_n2); // previous node coordinates

  int iNumRow = kMSamplesOrg.GetRows();
  CXYMatrix<float> kMSamplesVector(iNumRow,m_nb_features);
  for (int iRowInd = 0; iRowInd<iNumRow; iRowInd++)
  {
    // translation
    CXYVector<float> kVpoint = kMSamplesOrg.GetRow(iRowInd) + kV_2; //new node = previous coordinates + sample coordinates (works well)
    float fVector[9] = {kVpoint[0],kVpoint[1],kVpoint[2],fh1,fh2,fh3,interaction, fWeight, fh2c};
    CXYVector<float> kVector(m_nb_features,fVector);
    kMSamplesVector.SetRow(iRowInd, kVector);
  }

  return kMSamplesVector; //x y z + previous coordinates and previous weight, h1, h2, h3, fdistdiff
}
//------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *  Name:  GrowChainSegment_Single_RandSelect_2
 *  Description:  The function that goes through every candidate for one chain and chooses
 *  one based on growth function and fills vfSampleSelInd with the selected candidate
 * =====================================================================================
 */
void CXYSIS::GrowChainSegment_Single_RandSelect_2(
    CXYMatrix<float> kMSamples,
    vector<CXYMatrix<float> >& vkMWeightedSamples,
    int iChain,
    tree<CXYVector<float>>::iterator& itPrvNode,
    vector<int>& vfSampleSelInd,
    vector<float>& vfSampleSelBeta_t,
    vector<float>& vfSampleSelGamma_t,
    float& fSumBeta_t,
    float& fSumGamma_t){
    
    MTRand mtrand;
    vector<FloatIntPair> vfBeta_t;                                     /* growth function */
    vector<FloatIntPair> vfGamma_t;                                    /* target function */
    tree< CXYVector<float> >* ptr = GetTree();
    int iNodeInd = ptr->depth(itPrvNode);
    vfBeta_t.clear();
    vfGamma_t.clear();

    //computes all weights and errors for each potential point
    GrowChainSegment_CalBeta_t_2_Single(kMSamples, itPrvNode, iNodeInd, vfBeta_t, vfGamma_t, iChain);
    vkMWeightedSamples.push_back(kMSamples);
    fSumBeta_t = GrowChainSegment_SumVectorFloatIntPair(vfBeta_t);
    fSumGamma_t = GrowChainSegment_SumVectorFloatIntPair(vfGamma_t);
    float fConst_C = GrowChainSegment_BinSearchConstC_2(vfBeta_t,vfGamma_t); //what is happening here ?
    //cout<<"fConst_C : "<<fConst_C<<endl;
    
    //old way of selecting a candidate
    /*if(m_repeatNode[itPrvNode]!=0 || m_reminiscentRepeatNode[itPrvNode]!=0){
      vector<float> vfCumSum, vfMinCBeta1;
      float fCumSum = 0, fMinCBeta1;
      cout<<"fConst_C : "<<fConst_C<<endl;
      cout<<vfBeta_t[vfBeta_t.size() - 1].first<<endl;
      for (int iSample = 0; iSample < (int) vfBeta_t.size(); iSample ++) {
        fMinCBeta1 = min(fConst_C * (vfBeta_t[iSample].first),1.0f);
        fCumSum += fMinCBeta1;
        vfMinCBeta1.push_back(fMinCBeta1);
        vfCumSum.push_back(fCumSum);
      }

      // We can get a real number in the range 0 to 1, excluding 1
      float fRand = 0;
      unsigned int iCumSumInd = 0;
      fRand = 1 - mtrand.randExc();
      
      while (fRand >  vfCumSum[iCumSumInd]) {
        if (iCumSumInd == vfCumSum.size()-1) { //finit toujours ici donc pqs si random
          break;
        }
        ++iCumSumInd;
      }

      cout<<"fRand : "<<fRand<<" and vfCumSum : "<<iCumSumInd<<endl;
      vfSampleSelInd.push_back(vfBeta_t[iCumSumInd].second); 
      vfSampleSelBeta_t.push_back(vfBeta_t[iCumSumInd].first);
      vfSampleSelGamma_t.push_back(vfGamma_t[iCumSumInd].first);
      kMSamples[vfBeta_t[iCumSumInd].second][3] += (vfGamma_t[iCumSumInd].first/fSumGamma_t) / (vfBeta_t[iCumSumInd].first/fSumBeta_t);
      cout<<vfBeta_t[iCumSumInd].second<<" Random"<<endl;
    }*/

    //picks randomly a node between all the candidates that have the best beta_t value
    int nb_max = 0;
    float max = vfBeta_t[vfBeta_t.size()-1].first;
    while(vfBeta_t[vfBeta_t.size()-1-(++nb_max)].first==max){}; //counts the number of nodes that have the best beta_t value
    int index = vfBeta_t.size() - 1 - mtrand.randInt(nb_max-1); //picks an index between the nb_max first candidates which all have the best beta_t
    //cout<<"Range : "<<vfBeta_t[0].first<<" - "<<vfBeta_t[vfBeta_t.size()-1].first<<endl;
    //cout<<"Number of possibilities : "<<nb_max<<endl;
    vfSampleSelInd.push_back(vfBeta_t[index].second);
    vfSampleSelBeta_t.push_back(vfBeta_t[index].first);
    vfSampleSelGamma_t.push_back(vfGamma_t[index].first);

    /*for(int j=0; j<3; j++){
      cout<<(vkMWeightedSamples.back())[vfSampleSelInd.back()][j]<<" ";
    }
    cout<<"- ";*/
    //(vkMWeightedSamples.back())[vfBeta_t[index].second][7] += (vfGamma_t[index].first/fSumGamma_t) / (vfBeta_t[index].first/fSumBeta_t);
    // cout<<max<<" "<<vfBeta_t[index].first<<" "<<vfBeta_t[index].second<<" Random best"<<endl;
    //cout<<"weight calculation : "<<(vkMWeightedSamples.back())[vfBeta_t[index].second][3]<<endl;
    //cout<<"h1 : "<<(vkMWeightedSamples.back())[vfBeta_t[index].second][3]<<" h2 : "<<(vkMWeightedSamples.back())[vfBeta_t[index].second][4] - (*itPrvNode)[4]<<endl;
    //}

  //"returns" vfSampleSelInd, vfSampleSelBeta_t, vfSampleSelGamma_t with chosen candidate and respective errors
}
//------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *  Name:  GrowChainSegment_RandSelect_2
 *  Description:  Calls GrowChainSegment_Single_RandSelect_2 for every chain concerned
 *  (was designed for the parallel version otherwise, not that useful)
 * =====================================================================================
 */
void CXYSIS::GrowChainSegment_RandSelect_2(
  vector<CXYMatrix<float> >& vkMSamples,
  vector<CXYMatrix<float> >& vkMWeightedSamples,
  vector<int>& vfSampleSelInd,
  vector<float>& vfSampleSelBeta_t,
  vector<float>& vfSampleSelGamma_t,
  float& fSumBeta_t,
  float& fSumGamma_t
){
    tree< CXYVector<float> >* ptr = GetTree();
    tree< CXYVector<float> >::iterator pos, itPrvNode;
    //vector<CXYMatrix<float>> vkMWeightedSamples;


    vfSampleSelInd.clear(); //vector of randomly selected indexes
    vfSampleSelBeta_t.clear();
    vfSampleSelGamma_t.clear();
    //vkMWeightedSamples.clear();

    vector<float> vfSumBeta_t, vfSumGamma_t;
    
    itPrvNode = m_qPrvNodes_2[0];
    GrowChainSegment_Single_RandSelect_2(vkMSamples[0], vkMWeightedSamples, 0, itPrvNode, vfSampleSelInd, vfSampleSelBeta_t, vfSampleSelGamma_t, fSumBeta_t, fSumGamma_t);
    vfSumBeta_t.push_back(fSumBeta_t);
    vfSumGamma_t.push_back(fSumGamma_t);
    
    //"returns" vfSampleSelInd, vfSampleSelBeta_t, vfSampleSelGamma_t, vkMWeightedSamples
}
//-------------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  GrowChainSegment_CalBeta_t_2_Single
 *  Description:  Compute weight for all potential nodes when compared with itPrvNode chain
 * and stores it in kMSamples. Stores Gamma and Beta in associated vector for each sample
 * =====================================================================================
 */
void CXYSIS::GrowChainSegment_CalBeta_t_2_Single(
    CXYMatrix<float>& kMSamples,
    tree<CXYVector<float> >::iterator& itPrvNode,
    int iNodeInd,
    vector<FloatIntPair>& vBeta_t,
    vector<FloatIntPair>& vGamma_t,
    int iChain)
{
  //cout<<kMSamples.GetRows()<<endl;

  vector<vector<int> > vStartEndInd; //to sotre loop indexes for h3
  int iCurrentNodeInd = iNodeInd+1; //why not = ?
  FindContactStartEndNodeIndFromNodeInd_2(iCurrentNodeInd, vStartEndInd, iChain); //vector with loop paired indexes to check h3
  int iAttractionNum = vStartEndInd.size();

  vector<int> vConInd_1;
  vector<float> vConDist_1;
  //cout<<"iCurrentNodeInd : "<<iCurrentNodeInd<<endl;
  GetNodeIndAndDistFromNodeInd_2(iCurrentNodeInd, vConInd_1, vConDist_1, iChain); //gets the interaction constraints concerning this node
  float h1,h2,h3;
  float fBeta_t;                                                     /* growth function */
  float fGamma_t;                                                    /* target function */
  int iRowNum = iNodeInd + 1;
  CXYMatrix<float> kM_PrvNodePos(iRowNum, m_nb_features);                      /* first 3 element is x, y, z */
  GrowChainSegment_GetNodes_2(itPrvNode, kM_PrvNodePos); //gets the chain nodes in a matrix
  // compute weight for all potential nodes
  for (int iSample=0; iSample<kMSamples.GetRows(); iSample++) {
    CXYVector<float> vSample = kMSamples.GetRow(iSample);
    GrowChainSegment_h1Function_2_mod( //collision weight
        vSample, itPrvNode, iNodeInd);
    GrowChainSegment_h2Function_2_mod( //interaction weight
        vSample,vConInd_1, vConDist_1, kM_PrvNodePos);
    if (iAttractionNum != 0) {
      GrowChainSegment_h3Function_2_mod( //loop triangular inequality weight
        vSample, iNodeInd, vStartEndInd, kM_PrvNodePos, iChain);
    }

    h1 = vSample[3] ; //collisions at this node
    h2 = vSample[4] ; //missed interaction at this node 
    h3 = vSample[5] ; //triangular error

    //check vSample and kMSample row are the same, if not :
    //cout<<"vSample "<<iSample<<" : "<<vSample[0]<<" "<<vSample[1]<<" "<<vSample[2]<<" "<<vSample[3]<<" "<<vSample[4]<<" "<<vSample[5]<<" "<<vSample[6]<<endl;

    float fTemp = 0.5f; //*pow(2,iNodeInd);
    fGamma_t   = exp( - (m_fRho_1 * h1 +  m_fRho_2 * h2)/ fTemp*m_fTau_t ) ;   //gamze 3/18/2014
    fBeta_t  = exp( - (m_fRho_1 * h1 + m_fRho_2 * h2)/ (m_fTau_t) );
    vSample[7] = log(fBeta_t); //vSample[7]+log(fBeta_t);
    kMSamples.SetRow(iSample,vSample);

    
    //cout<<fBeta_t<<endl;
    vBeta_t.push_back(FloatIntPair(fBeta_t,iSample));
    vGamma_t.push_back(FloatIntPair(fGamma_t,iSample));
  }
  
  //"returns" vBeta_t, vGamma_t and kMSamples all for 1 chain (one itPrvNode)
}
//-------------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  GrowChainSegment_GetNodes_2
 *  Description:  Get all the nodes in the itNode's chain in gettree
 * =====================================================================================
 */
void CXYSIS::GrowChainSegment_GetNodes_2(tree<CXYVector<float > >::iterator itNode, CXYMatrix<float>& kM )
{
  tree<CXYVector<float > >* ptr = GetTree();
  int iRowInd = ptr->depth(itNode);
//  cout << iRowInd << endl;
  while (ptr->is_valid(itNode)){
    kM.SetRow(iRowInd, (*itNode));
//    cout << (*itNode)[0] <<"\t" <<(*itNode)[1] <<"\t"<<(*itNode)[2]<<endl;
    itNode = ptr->parent(itNode);
    iRowInd --;
  }
}
//------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *  Name:  GrowChainSegment_SumVectorFloatIntPair
 *  Description:  summation of vector floatintpair first element
 * =====================================================================================
 */
float CXYSIS::GrowChainSegment_SumVectorFloatIntPair ( vector<FloatIntPair>& rvFloatIntPair )
{
  float fSum = 0.0f;
  for (vector<FloatIntPair>::iterator it = rvFloatIntPair.begin(); it != rvFloatIntPair.end();
      it ++){
    fSum += (*it).first;
  }
  return fSum;
}		/* -----  end of function GrowChainSegment_SumVectorFloatIntPair  ----- */
//------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *  Name:  GrowChainSegment_BinSearchConstC_2
 *  Description: Function computes a certain constant C used in the previous version to 
 *  select a candidate node (and sorts vfbeta_t)
 * =====================================================================================
 */
float CXYSIS::GrowChainSegment_BinSearchConstC_2 ( vector<FloatIntPair>& vfBeta_t, vector<FloatIntPair>& vfGamma_t ){
  // vfBeta_t is changed
  sort(vfBeta_t.begin(), vfBeta_t.end(),FloatIntPairCompare());
  // vfGamma_t is changed
  GrowChainSegment_RearrangeVectorFloatIntPair(vfBeta_t,vfGamma_t);
  float l_max = 1; //m_iMmax;

  if ( vfBeta_t[vfBeta_t.size()-1].first * float(l_max) < 1.0f){
    float fSumBeta_t = GrowChainSegment_SumVectorFloatIntPair(vfBeta_t);
    float fConstC = l_max * fSumBeta_t;
    return fConstC;
  }

  CXYMath<float> kMath;
  int iSampleSize = (int) vfBeta_t.size();
  float fminnonzero = iSampleSize*l_max;
  int iRest = 1;
  for (int i=0; i<(int)vfBeta_t.size(); i++) {
    if (vfBeta_t[i].first >= kMath.ZERO_TOLERANCE) {
      fminnonzero = vfBeta_t[i].first;
      iRest = (int)vfBeta_t.size()-i;
      break;
    }
  }

  float fCons_Min = 0.0f;
  float fCons_Max = iRest/fminnonzero;

  float ftmp = 0;
  float fSumofMin = (float) (l_max+2);
  float fCons_Cur = fCons_Max;
  float fCons_Prv = fCons_Cur;

  int iMaxCellInd = iSampleSize-1;
  int iPrvCellInd = iMaxCellInd;
  int iCurCellInd = 0;

  bool bGreatThanOne = false;
  int iCount = 0;
  vector<float>::iterator it_begin, it_end;
  float iMmax = l_max ;

  while (fabs(fSumofMin-iMmax)>kMath.ZERO_TOLERANCE && (fCons_Max - fCons_Min)/2 > kMath.ZERO_TOLERANCE) {

    //  while (iPrvCellInd != iCurCellInd) {
    iPrvCellInd = iCurCellInd;
    //    if ( fSumofMin-m_iMmax > CXYMath<float>::ZERO_TOLERANCE )
    if ( fSumofMin-iMmax > 0)
    {
      fCons_Max = fCons_Cur;
    }else{
      fCons_Min = fCons_Cur;
    }
    fCons_Cur = fCons_Min + (fCons_Max - fCons_Min)/2.0f;


    fSumofMin = 0;
    bGreatThanOne = false;

    for (int i = 0 ; i < iSampleSize ; i++)
    {
      //    fSumofMin = fSumofMin + min(ftmp,1.0f);
      ftmp = fCons_Cur * vfBeta_t[i].first;
      if (bGreatThanOne) {
        fSumofMin += 1.0f;
      }else {
        if (ftmp < 1.0f) {
          fSumofMin += ftmp;
        }else {
          fSumofMin += 1.0f;
          bGreatThanOne = true;
          iCurCellInd = i;
        }
      }
    }

    if (kMath.AlmostEqual2sComplement(fCons_Prv,fCons_Cur,1) || iCount > min(max(100,(int)iMmax), 100) ) {
      break;
    }

    iCount ++;
    fCons_Prv = fCons_Cur;
  }
  return fCons_Cur;
} 		/* -----  end of function GrowChainSegment_BinSearchConstC_2  ----- */
//------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  GrowChainSegment_RearrangeVectorFloatIntPair
 *  Description:  rearrange vector FloatIntPair according to input FloatIntPair second index
 * =====================================================================================
 */
void CXYSIS::GrowChainSegment_RearrangeVectorFloatIntPair ( vector<FloatIntPair>& rvFIP_source, vector<FloatIntPair>& rvFIP_target)
{
  vector<FloatIntPair> vFIP_copy = rvFIP_target;
  int iInd_Tgt;
  for (unsigned int iInd_Src = 0; iInd_Src < rvFIP_source.size(); iInd_Src++){
    iInd_Tgt = rvFIP_source[iInd_Src].second;
    rvFIP_target[iInd_Src] =  vFIP_copy[iInd_Tgt];
  }
}		/* -----  end of function GrowChainSegment_RearrangeVectorFloatIntPair  ----- */
//------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  UpdateMap_TreeIt_Vec_Append
 *  Description:  Fills a map with the last nodes and their associated octree (chain they are
 *  attached to)
 * =====================================================================================
 */
void CXYSIS::UpdateMap_TreeIt_Vec_Append(
    tree<CXYVector<float> >::iterator oldpos,
    tree<CXYVector<float> >::iterator newpos,
    int iNodeInd,
    int back)
{
//	if (iSegInd == 0) {return;}

  tree<CXYVector<float> >* ptr = GetTree();
  tree<CXYVector<float> >::iterator PrvPos = oldpos;
  
  if (m_Treeit_VecpOctree_Map.find(oldpos) == m_Treeit_VecpOctree_Map.end()) { //if oldpos is not in the tree, adds the chain to the map
      boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator = boost::make_shared< OcTreeNodeAllocator< float , int > >();
      OcTree<float,int>* pOctree = new OcTree<float,int>(m_center,m_minimumVolumeSize,m_fNucleusSphereDiameter,allocator);
      pOctree->addPoint((*newpos)[0],(*newpos)[1],(*newpos)[2],1,m_maxDepth);
      while (ptr->is_valid(PrvPos))
      {
        pOctree->addPoint((*PrvPos)[0],(*PrvPos)[1],(*PrvPos)[2],1,m_maxDepth);
        PrvPos = ptr->parent(PrvPos);
      } //add the chain in the Octree 
      m_Treeit_VecpOctree_Map[newpos].push_back(pOctree);
      
  }
  else {
    // first one node of the segment
    m_Treeit_VecpOctree_Map[newpos] = m_Treeit_VecpOctree_Map[oldpos];
    //create a new entry with newpos as the key and rest of the chain as the mapped value
    OcTree<float,int>* pOctree;
    if(back==0){ //if we don't have to go back in the chain, just create a new entry [newpos : new octree with the new node]
      boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator = boost::make_shared< OcTreeNodeAllocator< float , int > >();
      pOctree = new OcTree<float,int>(m_center,m_minimumVolumeSize,m_fNucleusSphereDiameter,allocator);
      pOctree->addPoint((*newpos)[0],(*newpos)[1],(*newpos)[2],1,m_maxDepth);
      m_Treeit_VecpOctree_Map[newpos].push_back(pOctree);
    }
    else{ //else, just remove what should be removed in the entry [oldpos] and assign it to the new entry [newpos]
      for(int i=back-1; i>0; i--){
        m_Treeit_VecpOctree_Map[oldpos].pop_back();
      }
      m_Treeit_VecpOctree_Map[newpos]=m_Treeit_VecpOctree_Map[oldpos];
    }
  }
  //cout<<"OCTREE : "<<m_Treeit_VecpOctree_Map[newpos].size()<<endl;
}
//------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  UpdateMap_TreeIt_Vec_Clean
 *  Description:  Cleans the m_Treeit_VecpOctree_Map by removing the old nodes entries
 * =====================================================================================
 */
void CXYSIS::UpdateMap_TreeIt_Vec_Clean(void)
{
  tree<CXYVector<float>>* ptr = GetTree();
  for (vector<tree<CXYVector<float> >::iterator>::iterator
	  it=m_PrvTreeit_Vec.begin(); it != m_PrvTreeit_Vec.end(); it++) {
      TreeIterator_VecpOcTree_Map::iterator tvmit = m_Treeit_VecpOctree_Map.find(*it);
      if ( tvmit != m_Treeit_VecpOctree_Map.end() && ptr->depth(*it)>1) {
        m_Treeit_VecpOctree_Map.erase(tvmit);
      }
  }
  m_PrvTreeit_Vec.clear();
}
//------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  GrowChainSegment_h1Function_2_mod
 *  Description:  Computes the collision error for itPrvNode with its associated chain
 *  (using the octree stored in m_Treeit_VecpOctree_Map)
 * =====================================================================================
 */
void CXYSIS::GrowChainSegment_h1Function_2_mod(
    CXYVector<float>& kV_Vector_cur,
    tree< CXYVector<float> >::iterator itPrvNode,
    int iSegInd){
    float fh = 0;
    //cout<<iSegInd<<" - "<<m_Treeit_VecpOctree_Map[itPrvNode].size()<<endl;
	  for (int i=0; i<m_Treeit_VecpOctree_Map[itPrvNode].size(); i++) {
        fh += (IsCollision(m_Treeit_VecpOctree_Map[itPrvNode][i],kV_Vector_cur) ? 1 : 0);
    }
    kV_Vector_cur[3] = fh;
    //fh /= (iSegInd+1);
    //cout<<fh<<endl;
    //"returns" kV_Vector_cur with 3 modified //entre 0 et 1 selon les collisions
}
//------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  GrowChainSegment_h2Function_2_mod
 *  Description:  Computes the interaction error for itPrvNode using the constraints in
 *  rvConDist. Counts the number of times when the constraint is not respected and stores it
 *  in rKM_PrvNodePos
 * 
 * =====================================================================================
 */
void CXYSIS::GrowChainSegment_h2Function_2_mod(
    CXYVector<float>& kV_Vector_cur,
    vector<int>& rvConInd,
    vector<float>& rvConDist,
    CXYMatrix<float>& rKM_PrvNodePos)
{
  //cout<<" rvConDist length : "<<rvConDist.size()<<endl;
  CXYVector<float> kV_Vector_B = kV_Vector_cur;                              /* current node */
  float fPoint_B[3] = {kV_Vector_B[0], kV_Vector_B[1], kV_Vector_B[2]};
  CXYVector<float> kV_Point_B (3,fPoint_B);
  float fh = 0.0f;

  for (unsigned int iRow=0; iRow < rvConInd.size(); iRow++){
    int iNodeInd = rvConInd[iRow];
    float fPoint_A[3] = {rKM_PrvNodePos[iNodeInd][0],rKM_PrvNodePos[iNodeInd][1],rKM_PrvNodePos[iNodeInd][2]};
    CXYVector<float> kV_Point_A(3,fPoint_A);
    CXYVector<float> kVDiff = kV_Point_A - kV_Point_B;

    float fD_AB = kVDiff.Length();
    if (rvConDist[iRow]==m_fCollisionLength){ 
        fh+=(fD_AB > m_interaction_distance ? 1 : 0); 
    }
    else{ //if not interacting
        fh+=(fD_AB <= m_interaction_distance ? 1 : 0);
    }
  }

  /*if (rvConInd.size() > 0){
    fh /= rvConInd.size();
  }*/

  kV_Vector_cur[6] = rvConInd.size()+kV_Vector_cur[6]; //number of interactions until this node in the chain 
  kV_Vector_cur[8] = fh+kV_Vector_cur[8]; //number of missed interaction in total chain
  kV_Vector_cur[4] = fh; //number of missed interactions for this node only
  
  //"returns" kV_Vector_cur with 4, 6, 8 modified 
}
//------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  GrowChainSegment_h3Function_2_mod
 *  Description:  Computes the triangular error based on the loops the node should be in
 *  (stored in vStartEndInd)
 * 
 * =====================================================================================
 */
void CXYSIS::GrowChainSegment_h3Function_2_mod(
    CXYVector<float>& kV_Vector_cur,
    int iSegInd,
    vector<vector<int > >& vStartEndInd,
    CXYMatrix<float>& rKM_PrvNodePos,
    int iChain){
  /*-----------------------------------------------------------------------------
   *  Use triangle property
   *  | D_AB - D_AC | < D_B..C
   *  where B is current node, A is closest left node to B, A is already known
   *  C is the closest right node to B, C is not known
   *  D_AB is the distance of point A and point B
   *  D_AC is the distance of the condition from input
   *  D_B..C is is the summation of segments distance between point B and point C
   *-----------------------------------------------------------------------------*/

  int iLeftNodeInd, iRightNodeInd;
  vector<vector<int > >::iterator it;
  float fPoint_B[3] = {kV_Vector_cur[0],kV_Vector_cur[1],kV_Vector_cur[2]};

  CXYVector<float> kV_Point_B(3, fPoint_B);

  int iCurNodeInd = iSegInd + 1; //why not = ????
  float fh = 0.0f;
  for ( it=vStartEndInd.begin(); it < vStartEndInd.end(); it++){ /* how many loop the current node in */
    iLeftNodeInd = vStartEndInd[0][0];
    iRightNodeInd = vStartEndInd[0][1];

    float fPoint_A[3] = {rKM_PrvNodePos[iLeftNodeInd][0],rKM_PrvNodePos[iLeftNodeInd][1],rKM_PrvNodePos[iLeftNodeInd][2]};
    CXYVector<float> kV_Point_A(3,fPoint_A);
    CXYVector<float> kVDiff = kV_Point_A - kV_Point_B;
    float fD_AB = kVDiff.Length();
    float fD_BC = GrowChainSegment_GetTotalSegmentLengthFromTwoNodeInd_2(iCurNodeInd, iRightNodeInd);
    float fD_AC = GrowChainSegment_GetTwoNodesDistFromProfile_2(iLeftNodeInd, iRightNodeInd, iChain);

    if ( abs(fD_AB - fD_AC) > fD_BC ) {
      fh += 1;
    }
  }

  /*if (vStartEndInd.size() != 0){
    fh /= vStartEndInd.size();
  }*/

  kV_Vector_cur[5] = fh;

  //"returns" kV_Vector_cur with 5 modified
}
//------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  IsCollision
 *  Description:  Checks if the kV_point collides with any nodes of the octree
 * 
 * =====================================================================================
 */
bool CXYSIS::IsCollision(OcTree<float,int>* octree_, CXYVector<float>& kV_point)
{

  float x_ = kV_point[0];
  float y_ = kV_point[1];
  float z_ = kV_point[2];
//  cout << "x_ = " << x_ << ", y_ = " << y_ << ", z_ = " << z_ << endl ;
  return IsCollision(octree_, x_, y_, z_);

}
//------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  IsCollision
 *  Description:  Checks if the kV_point collides with any nodes of the octree
 * 
 * =====================================================================================
 */
bool CXYSIS::IsCollision(OcTree<float,int>* pOctree_, float x_, float y_, float z_)
{

//  Eigen::Matrix<float,4,1> currentpos_(x_,y_,z_,1);
  float min_x_ = max(x_ - m_dr, -m_fNucleusSphereDiameter);
  float min_y_ = max(y_ - m_dr, -m_fNucleusSphereDiameter);
  float min_z_ = max(z_ - m_dr, -m_fNucleusSphereDiameter);

  Eigen::Matrix<float,4,1> minPosition_ (min_x_,min_y_,min_z_,1);
  Eigen::Matrix<float,4,1> maxPosition_ (x_+m_dr,y_+m_dr,z_+m_dr,1);
  vector< OcTreeNode< float, int >* > nodes_;
  pOctree_->getAllNodesInVolumeOnDepth(nodes_,minPosition_,maxPosition_,m_maxDepth,true);

  if (nodes_.size() > 0) {
    bool bFlag = false;
    for (int i = 0; i< (int)nodes_.size(); i++) {
      Eigen::Matrix<float,4,1> point_ = (nodes_[i])->getPosition();
      float point_x_ = point_[0];
      float point_y_ = point_[1];
      float point_z_ = point_[2];
      float dist_ = sqrt((point_x_ - x_)*(point_x_ - x_) + (point_y_ - y_ )*(point_y_ - y_) + (point_z_ - z_)*(point_z_ - z_));
      if (dist_ < m_fCollisionLength) {
        bFlag = true;
        break;
        //        cout << x_ << " " << y_ << " " << z_ <<endl;
        //        cout << point_x_ << " " << point_y_ << " " << point_z_ <<endl;
        //        cout << dist_ << endl;
        //        cout << "collision happend";
      }
    }
    return bFlag;
  }else {
    return false;
  }

}
//-------------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  FindContactStartEndNodeIndFromNodeInd_2
 *  Description:  Stores in vStartEndInd the indexes of the loop in which the node iNodeInd
 *  is inside
 * 
 * =====================================================================================
 */
void CXYSIS::FindContactStartEndNodeIndFromNodeInd_2(int iNodeInd, vector<vector<int> > & vStartEndInd, int iChain)
{

  /*-----------------------------------------------------------------------------
   *  Find the loops which the point inside
   *  Include every loop starting and ending point
   *-----------------------------------------------------------------------------*/
  for (int i=0; i<m_MSegSegPval_2.GetRows(); i++) {
    if (m_MSegSegPval_2_distributed[i] != m_fCollisionLength*5 ) { //pas 0 la condition !!!
      if ((int)m_MSegSegPval_2[i][0] < iNodeInd && (int)m_MSegSegPval_2[i][1] > iNodeInd) {
        vector<int> vInd;
        vInd.push_back((int)m_MSegSegPval_2[i][0]);
        vInd.push_back((int)m_MSegSegPval_2[i][1]);
        vStartEndInd.push_back(vInd);
      }
    }
  }
}
//------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  GetNodeIndAndDistFromNodeInd_2
 *  Description:  Get node index and distance restriction from node index
 * =====================================================================================
 */
void CXYSIS::GetNodeIndAndDistFromNodeInd_2 (
    int iNodeInd, vector<int> & rvConInd, vector<float> & rvConDist, int iChain )
{
  int interval = 3; //since there has been a random selection at beginning in sets of 4 indexes each time, we have to take this in consideration when navigating through the matrix
  for (int i=0; i<m_MSegSegPval_2.GetRows(); i++) {
    if ((int)m_MSegSegPval_2[i][1] > iNodeInd+interval) {
      break;
    }
    if ((int)m_MSegSegPval_2[i][1] == iNodeInd) {
      rvConInd.push_back((int)m_MSegSegPval_2[i][0]);
      rvConDist.push_back(m_MSegSegPval_2_distributed[i]);
    }
  }
}		/* -----  end of function GetNodeIndAndDistFromNodeInd_2  ----- */
//------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  AddNewNodes
 *  Description:  Add the new node to the chain if the threshold are not overpassed, otherwise
 *  stays at old node to try again or goes back to retry the chain differently
 * =====================================================================================
 */
void CXYSIS::AddNewNodes(unordered_map<int,CXYVector<float> > &potNodes){
    tree<CXYVector<float> >* ptr = GetTree();
    tree< CXYVector<float> >::iterator itPrvNode, itmp;
    CXYVector<float> kVector(m_nb_features);
    deque<tree<CXYVector<float> >::iterator> vitmp;


    //for the beginning
    if(ptr->depth(m_qPrvNodes_2[0])<=1){
      itPrvNode = m_qPrvNodes_2[0];
      int iNodeInd = ptr->depth(itPrvNode);
      for(int iChain = 0; iChain<potNodes.size(); iChain++){    
        kVector = potNodes[iChain];
        //cout<<"nh1 : "<<kVector[3]/(ptr->depth(itPrvNode))<<" - "<<"nth2 : "<<kVector[8]<<" - "<<"sh2 : "<<kVector[4]<<" - "<<"h3 : "<<kVector[5]<<" - "<<"int : "<<kVector[6]<<" - "<<"weight : "<<kVector[7]<<" - ";
        if(kVector[3]!=-1){ //if not end of chain
            itmp = ptr->append_child(itPrvNode,kVector);
            vitmp.push_back(itmp);
            m_repeatNode[itmp]=0;
            m_reminiscentRepeatNode[itmp]=0;  
            UpdateMap_TreeIt_Vec_Append(itPrvNode, itmp, iNodeInd, 0);
        }
        else{
            GrowChainSegment_WritePtsArr(itPrvNode,++m_writtenFiles);
          }
      }
    }

    //for the rest of the chain
    else{
      int s = 0; //stay count
      for(int iChain = 0; iChain<potNodes.size(); iChain++){
          itPrvNode = m_qPrvNodes_2[iChain];
          kVector = potNodes[iChain];
          if(kVector[3]!=-1){ //if not end of chain
              int i = 0;
              //int den = (kVector[6]>0?kVector[6]:m_numInteraction);
              //cout<<"nh1 : "<<kVector[3]/(ptr->depth(itPrvNode))<<" - "<<"nth2 : "<<kVector[8]<<" - "<<"sh2 : "<<kVector[4]<<" - "<<"h3 : "<<kVector[5]<<" - "<<"int : "<<kVector[6]<<" - "<<"weight : "<<kVector[7]<<" - ";
              //cout<<<<" - "<<kVector[8]<<"/"<<(m_numInteraction)<<"="<<kVector[8]/m_numInteraction<<" - "<<m_repeatNode[itPrvNode]<<" - ";
              bool validWeight = (kVector[3]/(ptr->depth(itPrvNode))<m_fthreshold_A && kVector[8]/(m_numInteraction)<=m_fthreshold_B); //kVector[3]>threshold;
              if(validWeight){ //add node
                  itmp = ptr->append_child(itPrvNode,kVector);
                  vitmp.push_back(itmp);
                  m_repeatNode.erase(m_repeatNode.find(itPrvNode)); //erase last node count
                  m_repeatNode[itmp]=0;
                  m_reminiscentRepeatNode[itmp]=0;  
                  //cout<<ptr->depth(itmp)<<" - NODE ADDED";
                  int iNodeInd = ptr->depth(itPrvNode);
                  UpdateMap_TreeIt_Vec_Append(itPrvNode, itmp, iNodeInd, 0);
              }
              else if (!validWeight && m_repeatNode[itPrvNode]<10){ //stay at node
                /* does nothing except erasing the iterator m_PrvTreeit_Vec so that it won't be cleaned from m_Treeit_VecpOctree_Map */
                m_PrvTreeit_Vec.erase(m_PrvTreeit_Vec.begin()+iChain-s); //delete it so it won't be cleaned out of m_Treeit_VecpOctree_Map //ah check if ca marche bien maintenant aue le s est aussi ++ dans going back
                m_repeatNode[itPrvNode]=m_repeatNode[itPrvNode]+1;
                vitmp.push_back(itPrvNode);
                //cout<<ptr->depth(itPrvNode)<<" - STAYING";
                s++;
              }
              else{ //go back
                tree<CXYVector<float>>::iterator child; //used to erase it from the tree so there won't be any memory 
                itmp = itPrvNode;
                while(i<5 && ptr->depth(itmp)>1){
                  child = itmp;
                  itmp = ptr->parent(itmp); //parent of bad position
                  i++;
                  m_reminiscentRepeatNode.erase(m_reminiscentRepeatNode.find(child));
                }
                int n_backToThisNode = m_reminiscentRepeatNode[itmp];
                if(n_backToThisNode>10){
                  int reset = ptr->depth(itmp); //resets all the chain
                  while((i<reset+5 || m_reminiscentRepeatNode[itmp]>10) && ptr->depth(itmp)>1){ //&& (*itmp)[7]>0 but is it too restrictive ?
                  child = itmp;
                  itmp = ptr->parent(itmp); //parent of bad position
                  i++;
                  }
                }
                if(ptr->depth(itmp)==1){
                  m_nb_full_reset++;
                  cout<<"WENT BACK TO NODE AT DEPTH "<<ptr->depth(itmp)<<endl;
                }
                
                //at this point, itmp is the "new" node and child is its child in the current chain
                vitmp.push_back(itmp);
                m_repeatNode[itmp]=0; //reset it
                int iNodeInd = ptr->depth(ptr->parent(itmp));
                UpdateMap_TreeIt_Vec_Append(itPrvNode, itmp, iNodeInd, i+1); 
                m_reminiscentRepeatNode[itmp]=m_reminiscentRepeatNode[itmp]+1;
                m_repeatNode.erase(m_repeatNode.find(itPrvNode));

                //because we're going to erase itprvnode from the tree so we don't want it to be referenced in UpdateMap_TreeIt_Vec_Clean
                m_Treeit_VecpOctree_Map.erase(m_Treeit_VecpOctree_Map.find(itPrvNode));
                m_PrvTreeit_Vec.erase(m_PrvTreeit_Vec.begin()+iChain-s);
                s++;
                ptr->erase_children(child);
                ptr->erase(child);
              }
          }
          else{
            GrowChainSegment_WritePtsArr(itPrvNode,++m_writtenFiles); //if the chain is ended, save the file
          }
      }
    }

    UpdateMap_TreeIt_Vec_Clean();
    m_qPrvNodes_2.clear();
    potNodes.clear();

    while(!vitmp.empty()){
        itmp = vitmp.front();
        m_qPrvNodes_2.push_back(itmp);
        m_PrvTreeit_Vec.push_back(itmp);
        vitmp.pop_front();
    }
    //cout<<endl;

}

//------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  AddNewNodes
 *  Description:  Writes the files with all coordinates and h1 h2 h3 weight of each node 
 *  constituting the chain
 * =====================================================================================
 */
void CXYSIS::GrowChainSegment_WritePtsArr(tree< CXYVector<float> >::iterator lastNode, int chainInd)
{
  GrowChainSegment_ConfirmDistribution_2(lastNode); //fixes the target distribution if needed and confirms it by deleting it from m_MSegSegPval_2
  tree<CXYVector<float> > * ptr = GetTree();

  char buffer[1024];
  CXYMatrix<float> kM = CXYMatrix<float> (ptr->max_depth()+1, m_nb_features);
  GetOneChain(kM,lastNode);
  float err = (*lastNode)[8]/m_numInteraction; //number of missed interactions/number of interactions in chain
  m_Err_map[lastNode] = err;
  float weight = (*lastNode)[7];
  m_Weight_map[lastNode] = weight;
  m_Iteration_map[lastNode] = m_iteration/m_nb_nodes;

  //cout << "# Err= " << err << " # Weight= "<<weight<< endl;
  cout<< "Writing chain...";
  sprintf(buffer, "%s/%05d.pts",GetOutPath(),chainInd); //originally %04d.pts
  ostream* fp;
  fp = new std::ofstream(buffer);
  *fp << "# LogWeight= " << weight << endl;
  *fp << "# Err= " << err << endl;
  *fp << "# Iterations/NbNodes= " << m_iteration/m_nb_nodes << endl;
  *fp << "# NbChains= " << GetMmax() << endl;
  *fp << "# Threshold used= " << m_fthreshold_B << endl;
  char itembuff[1024];
  int iNumColumn = kM.GetColumns()-2; //we don't want to write the last column which probably doesn't have much sense after potentially fixing the target distribution, h3 is written but is actually not good either
  for (int iRow = 0; iRow < kM.GetRows(); iRow++){
    for ( int iCol=0; iCol < iNumColumn; iCol++){
      sprintf(itembuff,"%10.12f\t",kM[iRow][iCol]);
      *fp << itembuff;
      }
    *fp << endl;
  }
  delete fp;
  cout<<" finished."<<endl;
  cout << "# Err= " << err << endl;
  cout << "# LogWeight= " << weight << endl;
  cout << "# Iterations/nb_nodes= " << m_iteration/m_nb_nodes << endl;
  cout << "______________________________________________________________________________" << endl << endl;


}
//------------------------------------------------------------------------
// Gets the chain from itNode to the top of the tree
void CXYSIS::GetOneChain(CXYMatrix<float>& rkM, tree<CXYVector<float > >::iterator itNode)
{
  tree<CXYVector<float> >* ptr = GetTree();
  tree<CXYVector<float> >::iterator pos = itNode;

  int iNumRow = ptr->max_depth() + 1;

  int iRowInd = iNumRow - 1;
  while (ptr->is_valid(pos))
  {
    rkM.SetRow(iRowInd, (*pos));
    --iRowInd;
    pos = ptr->parent(pos);
  }
}
//------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  GrowChainSegment_GetTotalSegmentLengthFromTwoNodeInd_2
 *  Description:  Get total segment length between two node index
 * =====================================================================================
 */
float CXYSIS::GrowChainSegment_GetTotalSegmentLengthFromTwoNodeInd_2 ( int iStartNodeInd, int iEndNodeInd)
{
  float fdist;
  int iStartSegInd = iStartNodeInd;
  int iEndSegInd = iEndNodeInd - 1;
  fdist = GrowChainSegment_GetTotalSegmentLengthFromSegInd_2 (iStartSegInd, iEndSegInd);
  return (fdist);
}		/* -----  end of function GrowChainSegment_GetTotalSegmentLengthFromTwoNodeInd_2  ----- */
//------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  GrowChainSegment_GetTwoNodesDistFromProfile
 *  Description:  Get 2 nodes distance restriction from restriction profile
 * =====================================================================================
 */
float CXYSIS::GrowChainSegment_GetTwoNodesDistFromProfile_2 ( int iLeftNodeInd, int iRightNodeInd, int iChain)
{
  float fdist;

  for (int i=0; i<m_MSegSegPval_2.GetRows(); i++) {
    if ((int)m_MSegSegPval_2[i][0] == iLeftNodeInd && (int)m_MSegSegPval_2[i][1] == iRightNodeInd) {
      fdist = m_MSegSegPval_2[iChain][i];
      break;
    }
  }
  return (fdist);
}		/* -----  end of function GrowChainSegment_GetTwoNodesDistFromProfile  ----- */
//------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  GrowChainSegment_GetNumMiddlePoints_2
 *  Description:  Get the number of middle points from SegInd
 * =====================================================================================
 */
int CXYSIS::GrowChainSegment_GetNumMiddlePoints_2 ( int iSegInd )
{
  float fSegLen = m_vfSegLength[iSegInd];
  int iNum = (int) ( fSegLen / (m_fCollisionLength) );
  return (max(iNum,1));
}		/* -----  end of function GetNumMiddlePoint  ----- */
//------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  GetMiddlePoints
 *  Description:  Get Middle Points
 * =====================================================================================
 */
void CXYSIS::GetMiddlePoints(CXYVector<float>& StartPoint, CXYVector<float>& EndPoint, int iNumMiddlePoints, CXYMatrix<float> & MiddleEndPoints)
{

  float dx = (EndPoint[0] - StartPoint[0])/(float)(iNumMiddlePoints);
  float dy = (EndPoint[1] - StartPoint[1])/(float)(iNumMiddlePoints);
  float dz = (EndPoint[2] - StartPoint[2])/(float)(iNumMiddlePoints);

  for (int i =0; i< iNumMiddlePoints; i++) {
    MiddleEndPoints[i][0] = StartPoint[0] + (i+1) * dx;
    MiddleEndPoints[i][1] = StartPoint[1] + (i+1) * dy;
    MiddleEndPoints[i][2] = StartPoint[2] + (i+1) * dz;
  }

}		/* -----  end of function GetMiddlePoints  ----- */
//------------------------------------------------------------------------
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  GrowChainSegment_GetTotalSegmentLengthFromTwoSegInd
 *  Description:  Get total segment length between two segment index
 * =====================================================================================
 */
float CXYSIS::GrowChainSegment_GetTotalSegmentLengthFromSegInd_2 ( int iStartSegInd, int iEndSegInd)
{
  assert((unsigned int) iEndSegInd <= m_vfSegLength.size());
  float fdist = 0.0f;
  for (int i= iStartSegInd; i<= iEndSegInd; i++){
    fdist += m_vfSegLength[i];
  }
  return(fdist);
}		/* -----  end of function GrowChainSegment_GetTotalSegmentLengthFromTwoSegInd_2  ----- */
//------------------------------------------------------------------------
