#include "XYEnsemble.h"
using namespace spatialaggregate;

//----------------------------------------------------------------------------
CXYEnsemble::CXYEnsemble()
{
}

CXYEnsemble::CXYEnsemble(
	char* cOutPath,
	float fPersistenLength,
	float fCollisionLength,
	float fPackingDensity,
	float fNucleusSphereDiameter,
	int iNumNodes,
	int iNumSamplePoints,
	char* cStartEndFile,
    char* cSegLenFile,
	char* cContIndFile)
{
	m_cOutPath = new char[1024];
	SetOutPath(cOutPath);
	SetPersistenceLength( fPersistenLength );
	SetPackingDensity(fPackingDensity);
	SetNucleusSphereDiameter(fNucleusSphereDiameter);

	m_iNumSamplePoints = iNumSamplePoints;
//	if (strcmp(cStartEndFile, "") == 0){
//		// coarse version
//		// given number of nodes
//		SetNumNodes(iNumNodes); 
//		// set each segment length equal to persistence length
//		SetSegLengths();
//    SetContIndex();
//	} else {
//		// fine version
//		// read start end position file
//		// dynamic setting each segment length determined by mass density
//		SetSegLengths(cStartEndFile, "AVG");
//		SetNumNodes(m_vfSegLength.size()); // set node numbers
//    SetContIndex(cContIndFile);
//	}

  SetSegLengths(cSegLenFile,"Len");
  SetNumNodes(m_vfSegLength.size()+1); // set node numbers
  SetContIndex(cContIndFile);
  
	SetCollisionLength( fCollisionLength );
	m_pTChain = new tree< CXYVector<float> > ();
	m_pMSamplesOrg = new CXYMatrix<float>;
	SetSamplesOrg();
  
//  m_iNumMiddleEndPoints = max(int(m_fPersistenceLength / m_fCollisionLength)-1 ,1);
  m_iNumMiddleEndPoints = 0;
  
  cout << "persistence length = " << m_fPersistenceLength << endl;
  cout << "collision length = " << m_fCollisionLength << endl;
  cout << "number nodes = " << m_iNumNodes << endl;
//  cout << "m_iNumMiddleEndPoints = " << m_iNumMiddleEndPoints <<endl;
  

  // Construct octree
  m_center = Eigen::Matrix< float, 4, 1 > (0.0f,0.0f,0.0f,1);
  m_minimumVolumeSize = m_fCollisionLength/2;
  m_dr = m_fCollisionLength*2;
//  m_dr = m_fNucleusSphereDiameter/2;
  boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator  = 
    boost::make_shared< OcTreeNodeAllocator< float , int > >();
	m_pOctree = new OcTree<float,int>(m_center,m_minimumVolumeSize,m_fNucleusSphereDiameter,allocator);
  m_maxDepth = ceil(m_pOctree->depthForVolumeSize(m_minimumVolumeSize));
  cout << "minimumVolumeSize =" << m_minimumVolumeSize << endl;
  cout << "maxDistance =" << m_fNucleusSphereDiameter << endl;
  cout << "maxDepth = " <<  m_maxDepth << endl;
  
}

//----------------------------------------------------------------------------
CXYEnsemble::~CXYEnsemble(void)
{
  delete m_pOctree;
	delete m_cOutPath;
	delete m_pMSamplesOrg;
	delete m_pTChain;
}

//----------------------------------------------------------------------------
void CXYEnsemble::SetPersistenceLength(float fPL)
{
	m_fPersistenceLength = fPL;
}
//----------------------------------------------------------------------------
float CXYEnsemble::GetPersistenceLength(void)
{
	return m_fPersistenceLength;
}
//----------------------------------------------------------------------------
void CXYEnsemble::SetCollisionLength(float fCollision)
{
//	vector<float>& vfSegment = GetSegLengths();
//	float min_diameter = *(min_element(vfSegment.begin(), vfSegment.end()));
//	m_fCollisionLength = min(fCollision,min_diameter);
  m_fCollisionLength = fCollision;
}
//----------------------------------------------------------------------------
float CXYEnsemble::GetCollisionLength(void)
{
	return m_fCollisionLength;
}
//----------------------------------------------------------------------------
void CXYEnsemble::SetPackingDensity(float fPackingDensity)
{
	m_fPackingDensity = fPackingDensity;
}
//----------------------------------------------------------------------------
float CXYEnsemble::GetPackingDensity(void)
{
	return m_fPackingDensity;
}
//----------------------------------------------------------------------------
void CXYEnsemble::SetNucleusSphereDiameter(float fNucleusSphereDiameter)
{
	m_fNucleusSphereDiameter = fNucleusSphereDiameter;
}
//----------------------------------------------------------------------------
float CXYEnsemble::GetNucleusSphereDiameter(void)
{
	return m_fNucleusSphereDiameter;
}
//----------------------------------------------------------------------------
void CXYEnsemble::SetNumNodes(int iNumNodes){
	m_iNumNodes = iNumNodes;
}
//----------------------------------------------------------------------------
int CXYEnsemble::GetNumNodes(void){
	return m_iNumNodes;
}

//----------------------------------------------------------------------------
// Generate sphere sample points with radius persistence length.
void CXYEnsemble::SetSamplesOrg(void)
{
	CXYSO3Sequence sO3sequence(m_iNumSamplePoints);
	sO3sequence.SetSO3Sequence();
	(*m_pMSamplesOrg) = sO3sequence.GetSO3Sequence();
}
//----------------------------------------------------------------------------
CXYMatrix<float> CXYEnsemble::GetSamplesOrg()
{
	return (*m_pMSamplesOrg);
}

//----------------------------------------------------------------------------
//CXYMatrix<float> CXYEnsemble::NewXYZ(CXYVector<float> kV_0, CXYVector<float> kV_1, CXYVector<float> kV_2)
//{
//	CXYMatrix<float> kMXYZ(3,3);
//	CXYVector<float> kV_Z = kV_2 - kV_1;
//	CXYVector<float> kV_01 = kV_1 - kV_0;
//	kV_01.Normalize();
//	if (fabs(kV_Z.Dot(kV_01)) < CXYMath<float>::ZERO_TOLERANCE)
//	{ // two line parallel
//		kMXYZ = NewXYZ(kV_1,kV_2);
//	} else {
//		float a = kV_Z[0],  b = kV_Z[1],  c = kV_Z[2];
//		float x = kV_01[0], y = kV_01[1], z = kV_01[2];
//		
//		float fX[3] = {	y*c - z*b, z*a - x*c, x*b - y*a };
//		CXYVector<float> kV_X(3,fX);
//		kV_X.Normalize();
//		
//		x = kV_Z[0], y = kV_Z[1], z = kV_Z[2];
//		a = kV_X[0], b = kV_X[1], c = kV_X[2];
//		float fY[3] = {	y*c - z*b, z*a - x*c, x*b - y*a };
//		CXYVector<float> kV_Y(3,fY);
//		kV_Y.Normalize();
//		kMXYZ.SetRow(0, kV_X);
//		kMXYZ.SetRow(1, kV_Y);
//		kMXYZ.SetRow(2, kV_Z);
//	}
//	
//	return kMXYZ;
//}
////----------------------------------------------------------------------------
//CXYMatrix<float> CXYEnsemble::NewXYZ(CXYVector<float> kV_1, CXYVector<float> kV_2)
//{
//	CXYMatrix<float> kMXYZ(3,3);
//	CXYVector<float> kV_Z = kV_2 - kV_1;
//	kV_Z.Normalize();
//	float a = kV_Z[0], b= kV_Z[1], c = kV_Z[2];
//	
//	float x, y, z;
//	if (c > CXYMath<float>::ZERO_TOLERANCE ) 
//	{
//		x = 1; y = 0; z = -(a*x+b*y)/c;
//	} else if (b > CXYMath<float>::ZERO_TOLERANCE ) 
//	{
//		x = 0; z = 1; y = -(a*x+c*z)/b;
//	} else {
//		y = 1; z = 0; x = -(b*y+c*z)/a;
//	}
//	float fX[3] = {x,y,z};
//	CXYVector<float> kV_X(3,fX);
//	float fY[3] = {b*z-c*y, c*x-a*z, a*y-b*x};
//	CXYVector<float> kV_Y(3,fY);
//	
//	kMXYZ.SetRow(0, kV_X);
//	kMXYZ.SetRow(1, kV_Y);
//	kMXYZ.SetRow(2, kV_Z);
//	
//	
//	return kMXYZ;
//}

//----------------------------------------------------------------------------
//CXYMatrix<float> CXYEnsemble::GetRotMatrix(CXYMatrix<float> &rkMXYZ)
//{
//	CXYVector<float> kV_X = rkMXYZ.GetRow(0);
//	CXYVector<float> kV_Y = rkMXYZ.GetRow(1);
//	CXYVector<float> kV_Z = rkMXYZ.GetRow(2);
//	
//	CXYMatrix<float> kM_A(3,3);
//	float alpha, beta, gamma;
//	float Z3diff = fabs(1-kV_Z[2]*kV_Z[2]);
//	if ( Z3diff < CXYMath<float>::ZERO_TOLERANCE) 
//	{
//		alpha = acos( min(float(1.0), max(float(-1.0), kV_X[0])));
//		beta  = 0;
//		gamma = 0;
//	} else {		// http://www.macosxguru.net/article.php?story=20040210124637626
//		float Denorm = sqrt(Z3diff);
//		alpha = acos(min(float(1.0), max(float(-1.0), -kV_Z[1]/Denorm)));
//		beta  = acos(min(float(1.0), max(float(-1.0), -kV_Z[2])));
//		gamma = acos(min(float(1.0), max(float(-1.0), -kV_Y[2]/Denorm)));
//	}
//	
//	kM_A[0][0] = cos(gamma) * cos(alpha) - cos(beta) * sin(alpha) * sin(gamma);
//	kM_A[0][1] = cos(gamma) * sin(alpha) + cos(beta) * cos(alpha) * sin(gamma);
//	kM_A[0][2] = sin(gamma) * sin(beta);
//	kM_A[1][0] = -sin(gamma) * cos(alpha) - cos(beta) * sin(alpha) * cos(gamma);
//	kM_A[1][1] = -sin(gamma) * sin(alpha) + cos(beta) * cos(alpha) * cos(gamma);
//	kM_A[1][2] = cos(gamma) * sin(beta);
//	kM_A[2][0] = sin(beta) * sin(alpha);
//	kM_A[2][1] = -sin(beta) * cos(alpha);
//	kM_A[2][2] = cos(beta);
//	
//	kM_A.GetInverse(kM_A);
//	
//	
//	return kM_A;
//}

//----------------------------------------------------------------------------
void CXYEnsemble::InitializeChain()
{
  // renew octree
  if (m_pOctree->root_) {
    delete m_pOctree;
  }

  boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator  = 
    boost::make_shared< OcTreeNodeAllocator< float , int > >();
  m_pOctree = new OcTree<float,int>(m_center,m_minimumVolumeSize,m_fNucleusSphereDiameter,allocator);
  
	tree< CXYVector<float> >* ptr = GetTree();
	// erase all nodes if not empty
	if (!ptr->empty()) {
		ptr->clear();
	}
	
	tree< CXYVector<float> >::iterator top, root;
	
	top = ptr->begin();
	
	// first node position
//  float fCoord[3] = {
//    0.0f, 0.0f, 0.0f
//  };
//  CXYVector<float> kV_startpoint (3,fCoord);
  CXYVector<float> kV_startpoint = RndSetStartPoint();

	root = ptr->insert(top, kV_startpoint);
  
  m_pOctree->addPoint(kV_startpoint[0],kV_startpoint[1],kV_startpoint[2],1,m_maxDepth);
  
  m_LogWeight = 0;
  m_vLogWeight.empty();
  m_vLogWeight.push_back(m_LogWeight);
  
//  cout << "StartingPoint: " << kV_startpoint[0] << " " << kV_startpoint[1] << " " << kV_startpoint[2] << endl;
	
//	CXYMatrix<float> kMSamplesOrg = GetSamplesOrg();
//	kMSamplesOrg *= GetSegLength(0); // index is 0, first one segment
//	int iSampleSize = kMSamplesOrg.GetRows();
//  
//	int iRnd;
//	CXYVector<float> kV_Point(3);
//	// second node position
//	MTRand mtrand;
//	while (1) {
//		iRnd = mtrand.randInt(iSampleSize-1);  // integer in [0,iSampleSize-1]
//		kV_Point = kMSamplesOrg.GetRow(iRnd) + kV_startpoint; // second point
//    
//    //----------------------------
//    // for debug
////    fCoord[0] = -129.225;
////    fCoord[1] = 129.225;
////    fCoord[2] = 881.25 ;
////    kV_Point = CXYVector<float> (3,fCoord);
//    //----------------------------
//    
//    CXYMatrix<float> kM_MiddlePoints(m_iNumMiddleEndPoints,3);
//    GetMiddlePoints(kV_startpoint, kV_Point, kM_MiddlePoints);
//    
//		if (IsInsideSphere(kV_Point) ) {
//			ptr->append_child(root, kV_Point);
////      m_pOctree->addPoint(kV_Point[0], kV_Point[1],kV_Point[2],1,m_maxDepth);
//      for (int iNumPoint= m_iNumMiddleEndPoints-1; iNumPoint >=0; iNumPoint--) {
//        m_pOctree->addPoint(kM_MiddlePoints[iNumPoint][0], kM_MiddlePoints[iNumPoint][1],kM_MiddlePoints[iNumPoint][2],1,m_maxDepth);
//      }
////      cout << kV_Point[0] << " " << kV_Point[1] << " " << kV_Point[2] << endl;
//			break;
//		}
//	}
}
//----------------------------------------------------------------------------
tree< CXYVector<float> >* CXYEnsemble::GetTree()
{
	return m_pTChain;
}
//----------------------------------------------------------------------------
CXYVector<float> CXYEnsemble::RndSetStartPoint(void)
{
	float fdiameter = GetNucleusSphereDiameter();
	float fradius = fdiameter/2;
	MTRand mtrand;
	
	// Check if the random point located in the neuclus sphere.
	// When one point satisfy the condition, we accept it.
	float fCoord[3];
	while (1) {
		fCoord[0] = mtrand.randExc(fradius);
		fCoord[1] = mtrand.randExc(fradius);
		fCoord[2] = mtrand.randExc(fradius);
		if (IsInsideSphere(fCoord)){
			break;
		}
	}
//	cout << x << ";" << y << ";" << z<<";" <<endl;
//	cout << x*x+y*y+z*z << endl;
//	cout << fradius2 << endl;
	CXYVector<float> kV (3,fCoord);
	return kV;
}
//----------------------------------------------------------------------------
bool CXYEnsemble::IsInsideSphere(CXYVector<float> kV_point)
{
	float fradius = GetNucleusSphereDiameter()/2;
	float fradius2 = fradius*fradius;
	bool flag;
	if (kV_point.SquaredLength() < fradius2) {
		flag = true;
	}else {
		flag = false;
	}
	return flag;
	
}
//----------------------------------------------------------------------------
bool CXYEnsemble::IsInsideSphere(float* fCoord)
{
	float fradius = GetNucleusSphereDiameter()/2;
	float fradius2 = fradius*fradius;
	bool flag;
	if (fCoord[0]*fCoord[0] + fCoord[1]*fCoord[1] + fCoord[2]*fCoord[2]
			< fradius2) {
		flag = true;
	}else {
		flag = false;
	}
	return flag;
}

//----------------------------------------------------------------------------
bool CXYEnsemble::IsCollision(CXYVector<float> kV_point)
{
  
  float x_ = kV_point[0];
  float y_ = kV_point[1];
  float z_ = kV_point[2];

  float min_x_ = max(x_ - m_dr, -m_fNucleusSphereDiameter);
  float min_y_ = max(y_ - m_dr, -m_fNucleusSphereDiameter);
  float min_z_ = max(z_ - m_dr, -m_fNucleusSphereDiameter);
  Eigen::Matrix<float,4,1> minPosition_ (min_x_,min_y_,min_z_,1);
  Eigen::Matrix<float,4,1> maxPosition_ (x_+m_dr,y_+m_dr,z_+m_dr,1);
  vector< OcTreeNode< float, int >* > nodes_;
  m_pOctree->getAllNodesInVolumeOnDepth(nodes_,minPosition_,maxPosition_,m_maxDepth,true);

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

//----------------------------------------------------------------------------
bool CXYEnsemble::GrowOneChain()
{
	bool Flag = true;
	InitializeChain();
	MTRand mtrand;

	int iNumNodes = GetNumNodes();
	tree< CXYVector<float> >* ptr = GetTree();
	tree< CXYVector<float> >::pre_order_iterator node;


//  iNumNodes = 3;
	for (int i =1 ; i < iNumNodes; i++) {
//		cout << i << endl;
		
		node = ptr->end(); 
		node --; // go to last node

		CXYMatrix<float> kMSamplesPoints = GetSamplesOrg();
		kMSamplesPoints *= GetSegLength(i-1); // because initialize chain already done
    m_iNumMiddleEndPoints = max((int)(GetSegLength(i-1)/m_fCollisionLength) - 1,1);
		
		//int iSampleSize = kMSamplesPoints.GetRows();
    
    vector<int> GoodPointInd;
    vector<int> NoCollisionPointInd;
    GetGoodPoints( (*node), kMSamplesPoints, GoodPointInd, NoCollisionPointInd);
    int GoodPointSize = GoodPointInd.size();
    int NoCollisionPointSize = NoCollisionPointInd.size();
    if (GoodPointSize == 0 || NoCollisionPointSize == 0) {
      Flag = false;
      break;
    }
    
    
//    cout << "i = " << i << ",GoodPointSize = " << GoodPointSize << endl;
    m_LogWeight += log(NoCollisionPointSize);
    m_vLogWeight.push_back(m_LogWeight);
    
		int iRnd = mtrand.randInt(GoodPointSize-1); // integer in [0, GoodPointSize-1]
		CXYVector<float> kV_Point(3);


    CXYMatrix<float> kM_MiddlePoints(m_iNumMiddleEndPoints,3);
    
    kV_Point = kMSamplesPoints.GetRow(GoodPointInd[iRnd]) + (*node);
    GetMiddlePoints((*node), kV_Point, kM_MiddlePoints);

    ptr->append_child(node, kV_Point);
    for (int iNumPoint=m_iNumMiddleEndPoints-1; iNumPoint >= 0; iNumPoint --) {
      m_pOctree->addPoint(kM_MiddlePoints[iNumPoint][0],kM_MiddlePoints[iNumPoint][1],kM_MiddlePoints[iNumPoint][2],1,m_maxDepth);
    }

    // 
//    cout << "---------------\n";
//    for (int iPoint = 0; iPoint < GoodPointSize; iPoint ++) {
//      kV_Point = kMSamplesPoints.GetRow(GoodPointInd[iPoint]) + (*node);
//      cout << kV_Point[0] << " " << kV_Point[1] << " " << kV_Point[2] << endl;
//    }
//    cout << "---------------all\n";
//    for (int iPoint = 0; iPoint < kMSamplesPoints.GetRows(); iPoint ++) {
//      kV_Point = kMSamplesPoints.GetRow(iPoint) + (*node);
//      cout << kV_Point[0] << " " << kV_Point[1] << " " << kV_Point[2] << endl;
//    }
    
    
//		int n_trial = max(1000,iNumNodes);
//		while (n_trial > 0) {
//			iRnd = mtrand.randInt(iSampleSize-1);  // integer in [0,iSampleSize-1]
//			kV_Point = kMSamplesPoints.GetRow(iRnd) + (*node); // next point position
//      
//      // Get Middle and End points
//      GetMiddlePoints((*node), kV_Point, kM_MiddlePoints);
//      
//      if(IsSatisfyCondition(kM_MiddlePoints)){
//        ptr->append_child(node, kV_Point);
//        for (int iNumPoint=m_iNumMiddleEndPoints-1; iNumPoint >= 0; iNumPoint --) {
//          m_pOctree->addPoint(kM_MiddlePoints[iNumPoint][0],kM_MiddlePoints[iNumPoint][1],kM_MiddlePoints[iNumPoint][2],1,m_maxDepth);
//        }
//        break;
//      }
//      
//			// check if inside neuclus sphere and there is no collision
////			if (IsInsideSphere(kV_Point) && (! IsCollision(kV_Point))) {
//////        cout << kV_Point[0] << " " << kV_Point[1] << " " << kV_Point[2] << endl;
////				ptr->append_child(node, kV_Point);
////        m_pOctree->addPoint(kV_Point[0],kV_Point[1],kV_Point[2],1,m_maxDepth);
////				break;
////			}
//			n_trial --;
//		}
//
//
//		// do not continue our growing chain
//		if (n_trial == 0){
//			Flag = false;
//			break;
//		}
	}
  cout << m_LogWeight << endl;
	return Flag;
}
//----------------------------------------------------------------------------
void CXYEnsemble::WriteChain(char* fn)
{
//http://stdcxx.apache.org/doc/stdlibug/34-2.html	
	
	tree<CXYVector<float> > * ptr = GetTree();
	tree< CXYVector<float> >::iterator pos = ptr->begin();
	char buffer[1024];
	
	ostream* fp;
	if (strcmp(fn,"")) {
    sprintf(buffer, "%s/%s", GetOutPath(),fn);
		fp = new std::ofstream(buffer);
	}else {
		fp = &std::cout;
	}

  *fp << "# LogWeight= " << m_LogWeight << endl;
	while (ptr->is_valid(pos)) {
		*fp << (*pos)[0] <<"\t" <<(*pos)[1] <<"\t"<<(*pos)[2]<<endl;
		++pos;
	}

	if (fp != &std::cout)
		delete fp;	
}
//----------------------------------------------------------------------------
void CXYEnsemble::WriteDistance(char* fn)
{
	//http://stdcxx.apache.org/doc/stdlibug/34-2.html	
	
	tree<CXYVector<float> > * ptr = GetTree();
	tree< CXYVector<float> >::iterator root = ptr->begin();
//	tree< CXYVector<float> >::fixed_depth_iterator pos1, pos2;
//	float deltaX, deltaY, deltaZ, length;
	char buffer[1024];
	
	ostream* fp;
	if (strcmp(fn,"")) {
    sprintf(buffer, "%s/%s", GetOutPath(),fn);
    fp = new std::ofstream(buffer);
	}else {
		fp = &std::cout;
	}

  // Calculate mean and variance, 
  tree< CXYVector<float> >::fixed_depth_iterator pos;
  int SegmentInd, StartNodeInd, EndNodeInd;
  float mean_x_, mean_y_, mean_z_;
  float var_x_, var_y_, var_z_;
  
  for (int irow = 0; irow < m_MContInd.GetRows(); irow++) {

    SegmentInd = m_MContInd[irow][0];
    StartNodeInd = m_MContInd[irow][1] - 1;
    EndNodeInd = m_MContInd[irow][2] - 1;
    
    mean_x_ = mean_y_ = mean_z_ = 0.0f;
    var_x_ = var_y_ = var_z_ = 0.0f;
    for (int iInd = StartNodeInd; iInd <= EndNodeInd; iInd++) {
      pos = ptr->begin_fixed(root, iInd);
      mean_x_ += (*pos)[0];
      mean_y_ += (*pos)[1];
      mean_z_ += (*pos)[2];
      var_x_  += (*pos)[0]*(*pos)[0];
      var_y_  += (*pos)[1]*(*pos)[1];
      var_z_  += (*pos)[2]*(*pos)[2];
    }
    
    int iNumNodes = EndNodeInd - StartNodeInd + 1;
    mean_x_ = mean_x_ / iNumNodes;
    mean_y_ = mean_y_ / iNumNodes;
    mean_z_ = mean_z_ / iNumNodes;

    var_x_  = var_x_  / iNumNodes;
    var_y_  = var_y_  / iNumNodes;
    var_z_  = var_z_  / iNumNodes;
    var_x_  = var_x_ - mean_x_ * mean_x_;
    var_y_  = var_y_ - mean_y_ * mean_y_;
    var_z_  = var_z_ - mean_z_ * mean_z_;
    
    float r_x_, r_y_, r_z_;
    r_x_ = sqrt(var_x_) * 3;
    r_y_ = sqrt(var_y_) * 3;
    r_z_ = sqrt(var_z_) * 3;
    
    if (iNumNodes == 1) {
      r_x_ = r_y_ = r_z_ = m_fPersistenceLength / 2;
    }
    
    *fp << SegmentInd << " " \
      << mean_x_ << " " \
      << mean_y_ << " " \
      << mean_z_ << " " \
      << r_x_ << " "  \
      << r_y_ << " "  \
      << r_z_ << " "  \
      << endl;
  }
  
	if (fp != &std::cout)
		delete fp;
}
//----------------------------------------------------------------------------
void CXYEnsemble::SetOutPath(char* cPathName)
{
	CXYFile::MakeDirectory(cPathName,0755);
  cout << "make directory " << cPathName << endl;
	strcpy(m_cOutPath,cPathName);
}

//----------------------------------------------------------------------------
char* CXYEnsemble::GetOutPath(void)
{
	return m_cOutPath;
}
//----------------------------------------------------------------------------
// void CXYEnsemble::SetSegLengths(char* cStartEndFile){
// 	CXYMatrix<float> kMStartEnd = CXYFile::ReadMatrix(cStartEndFile);
// 	for (int i= 0; i<kMStartEnd.GetRows(); i++) {
// 		int baselen = (int) ( kMStartEnd[i][1] - kMStartEnd[i][0] );
// 		m_vfSegLength.push_back(((float)baselen * GetPackingDensity()));
// 	}
// }
//----------------------------------------------------------------------------
void CXYEnsemble::SetSegLengths(char* cStartEndFile,const char* cMethod){
	CXYMatrix<float> kMStartEnd = CXYFile::ReadMatrix(cStartEndFile);
	m_vfSegLength.empty();
	if (strcmp(cMethod, "AVG") == 0) {
		// average different node length
		float len = 0L;
		for (int i= 0; i<kMStartEnd.GetRows(); i++) {
			len = len + ( kMStartEnd[i][1] - kMStartEnd[i][0] + 1 ) * GetPackingDensity();
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
	} else if (strcmp(cMethod, "Len") == 0) {
		for (int i= 0; i<kMStartEnd.GetRows(); i++) {
			m_vfSegLength.push_back((float)kMStartEnd[i][0]);
		}
  } else {
		cerr	<< "Please type Method" << endl;
		exit(-1);
	}
}
//----------------------------------------------------------------------------

vector<float> & CXYEnsemble::GetSegLengths(void){
	return m_vfSegLength;
}
//----------------------------------------------------------------------------
float CXYEnsemble::GetSegLength(int ind){
	return m_vfSegLength[ind];
}
//----------------------------------------------------------------------------
void CXYEnsemble::SetSegLengths(void){
  m_vfSegLength.empty();
	for (int i= 0; i<GetNumNodes(); i++) {
		m_vfSegLength.push_back(GetPersistenceLength());
	}
}
//----------------------------------------------------------------------------
void CXYEnsemble::SetContIndex(char* cContIndFile){
//	CXYVector<int> vtmp = CXYFile::ReadVectorInt(cContIndFile);
//	for (int i =0; i<vtmp.GetSize(); i++) {
//		// index start from 0
//		m_MContInd.push_back(vtmp[i]-1);
//	}
  m_MContInd = CXYFile::ReadMatrix(cContIndFile);
}
//----------------------------------------------------------------------------
void CXYEnsemble::SetContIndex(void){
  m_MContInd.SetSize(m_iNumNodes,3);
  for (int i=0; i<m_iNumNodes; i++) {
    m_MContInd[i][0] = i+1;
    m_MContInd[i][1] = i+1;
    m_MContInd[i][2] = i+1;
  }
}
//----------------------------------------------------------------------------
void CXYEnsemble::GetMiddlePoints(CXYVector<float>& StartPoint, CXYVector<float>& EndPoint, CXYMatrix<float> & MiddleEndPoints){

  float dx = EndPoint[0] - StartPoint[0];
  float dy = EndPoint[1] - StartPoint[1];
  float dz = EndPoint[2] - StartPoint[2];
  
  for (int i =0; i< m_iNumMiddleEndPoints; i++) {
    MiddleEndPoints[i][0] = StartPoint[0] + ((i+1)/(float)m_iNumMiddleEndPoints) * dx;
    MiddleEndPoints[i][1] = StartPoint[1] + ((i+1)/(float)m_iNumMiddleEndPoints) * dy;
    MiddleEndPoints[i][2] = StartPoint[2] + ((i+1)/(float)m_iNumMiddleEndPoints) * dz;
  }
  
}
//----------------------------------------------------------------------------
bool CXYEnsemble::IsSatisfyCondition(CXYMatrix<float> & MiddleEndPoints){
  bool flag = true;
  for (int i=0; i< m_iNumMiddleEndPoints; i++) {
    if (! IsInsideSphere(MiddleEndPoints.GetRow(i))  ) {
      flag = false;
      break;
    }
    if (IsCollision(MiddleEndPoints.GetRow(i))){
//      cout << MiddleEndPoints[m_iNumMiddleEndPoints-1][0] << " " << MiddleEndPoints[m_iNumMiddleEndPoints-1][1] << " " << MiddleEndPoints[m_iNumMiddleEndPoints-1][2] << endl;
      flag = false;
      break;
    }
  }
  return flag;
}
//----------------------------------------------------------------------------
bool CXYEnsemble::IsCollision(CXYMatrix<float> & MiddleEndPoints){
  bool flag = false;
  for (int i=0; i< m_iNumMiddleEndPoints; i++) {
    if ( IsCollision(MiddleEndPoints.GetRow(i))){
      //      cout << MiddleEndPoints[m_iNumMiddleEndPoints-1][0] << " " << MiddleEndPoints[m_iNumMiddleEndPoints-1][1] << " " << MiddleEndPoints[m_iNumMiddleEndPoints-1][2] << endl;
      flag = true;
      break;
    }
  }
  return flag;
}

//----------------------------------------------------------------------------
void CXYEnsemble::GetGoodPoints( CXYVector<float>& prvnode, CXYMatrix<float>& kMSamplesPoints,vector<int>& GoodPointInd, vector<int>& NoCollisionPointInd){

 // cout << prvnode[0] << " " << prvnode[1] << " " << prvnode[2] << endl;
  
  CXYVector<float> kV_Point(3);
  CXYMatrix<float> kM_MiddlePoints(m_iNumMiddleEndPoints,3);
  
  int iSampleSize = kMSamplesPoints.GetRows();
//  cout << "===================="<< iSampleSize << endl;
  for (int i = 0; i < iSampleSize; i++) {
    kV_Point = kMSamplesPoints.GetRow(i) + prvnode;
    
//    cout << kV_Point[0] << " " << kV_Point[1] << " " << kV_Point[2] << endl;
    
    GetMiddlePoints(prvnode, kV_Point, kM_MiddlePoints);
//    cout << "i = " << i <<endl;
    if (IsSatisfyCondition(kM_MiddlePoints)){
      GoodPointInd.push_back(i);
    }
    // No collision
    if (! IsCollision(kM_MiddlePoints)){
      NoCollisionPointInd.push_back(i);
    }
  }
//  cout << "===================="<< iSampleSize << endl;
  
}

//----------------------------------------------------------------------------
void CXYEnsemble::WriteWeight(char* fn)
{
	char buffer[1024];
	
	ostream* fp;
	if (strcmp(fn,"")) {
    sprintf(buffer, "%s/%s", GetOutPath(),fn);
		fp = new std::ofstream(buffer);
	}else {
		fp = &std::cout;
	}
  
  *fp << "# ChainLogWeight= " << m_LogWeight << endl;
  for (int i=0; i < (int)m_vLogWeight.size(); i++) {
    *fp << i+1 << " " << m_vLogWeight[i] << endl;
  }
  
	if (fp != &std::cout)
		delete fp;	
}


//----------------------------------------------------------------------------
void CXYEnsemble::WriteContactDistance(char* fn)
{
	//http://stdcxx.apache.org/doc/stdlibug/34-2.html	
	
	tree<CXYVector<float> > * ptr = GetTree();
	tree< CXYVector<float> >::iterator root = ptr->begin();
  //	tree< CXYVector<float> >::fixed_depth_iterator pos1, pos2;
  //	float deltaX, deltaY, deltaZ, length;
	char buffer[1024];
	
	ostream* fp;
	if (strcmp(fn,"")) {
    sprintf(buffer, "%s/%s", GetOutPath(),fn);
    fp = new std::ofstream(buffer);
	}else {
		fp = &std::cout;
	}
  
  *fp << "# LogWeight= " << m_LogWeight << endl;
  // Calculate mean and variance, 
  tree< CXYVector<float> >::fixed_depth_iterator firstpos,secondpos;
  int FirstInd, SecondInd;
  float x1_,y1_,z1_;
  float x2_,y2_,z2_;
  float xdiff_,ydiff_,zdiff_;
  float length;
  
  for (int irow = 0; irow < m_MContInd.GetRows(); irow++) {
    FirstInd = m_MContInd[irow][0];
    SecondInd = m_MContInd[irow][1];
    
    firstpos = ptr->begin_fixed(root, FirstInd-1);
    secondpos = ptr->begin_fixed(root, SecondInd-1);
    
    x1_ = (*firstpos)[0]; y1_ = (*firstpos)[1]; z1_ = (*firstpos)[2];
    x2_ = (*secondpos)[0]; y2_= (*secondpos)[1];z2_ = (*secondpos)[2];
    
    xdiff_ = (x2_ - x1_);
    ydiff_ = (y2_ - y1_);
    zdiff_ = (z2_ - z1_);
    
    length = sqrt(xdiff_*xdiff_+ydiff_*ydiff_+zdiff_*zdiff_);
    
    *fp << FirstInd << "\t" \
    << SecondInd << "\t " \
    << length << " " \
    << endl;
  }
  
	if (fp != &std::cout)
		delete fp;
}


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
