/* Test octree 
*/

#include <iostream>
#include "tree.hh"
#include "XYSO3Sequence.h"
#include "XYVector.h"
#include "octree.h"
#include <vector>

using namespace spatialaggregate;
int main (int argc, char *argv[]){
  Eigen::Matrix< float, 4, 1 > center(0.0f,0.0f,0.0f,1);
//  Eigen::Matrix< float, 4, 1 > center(0.1f,0.2f,0.3f,1);
  float minimumVolumeSize = 1; 
  float maxDistance = pow(2.0,8); 
  boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator  = boost::make_shared< OcTreeNodeAllocator< float , int > >();
  OcTree<float,int> octree_(center,minimumVolumeSize,maxDistance,allocator);
  int maxDepth = ceil(octree_.depthForVolumeSize(1));
  cout << "maxDepth = " <<  maxDepth << endl;
  
  
//  for (int i=0; i<4; i++) {
//    boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator  =  boost::make_shared< OcTreeNodeAllocator< float , int > >();        
//    OcTree<float,int>* pOctree = new OcTree<float,int> (center,minimumVolumeSize,maxDistance,allocator) ;
//    cout << "ooo= " << allocator << endl;        
//    cout << "aaa= " << pOctree << endl;
//       
//  }
//  return(0);
  
  vector<int *> vpint;
  int * a;
  a = new int;
  cout << "before " << a << endl;
  vpint.push_back(a);
  cout << "after " << vpint[0] << endl;
  
    vector<OcTree<float,int>* > vpoctree;
    for (int i=0; i<4; i++){
        boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator_test  =  boost::make_shared< OcTreeNodeAllocator< float , int > >();        
        OcTree<float,int>* pOctree = new OcTree<float,int> (center,minimumVolumeSize,maxDistance,allocator_test) ;
        cout << "BEFORE= " << pOctree << endl;
        vpoctree.push_back(pOctree);
        cout <<"AFTER= " << vpoctree[i] << endl;
        vpoctree[i]->addPoint(1256.19,  1478.53, -1194.18,1 , maxDepth);
    }
    exit(0);

    // Add point
  OcTreeNode<float,int>* poctreenode_;

  poctreenode_ = octree_.addPoint(0, 0, 0, 1, maxDepth);
  poctreenode_ = octree_.addPoint(1, 0, 0, 1, maxDepth);
  poctreenode_ = octree_.addPoint(1, 1, 0, 1, maxDepth);
  poctreenode_ = octree_.addPoint(0, 1, 0, 1, maxDepth);
 
  poctreenode_ = octree_.addPoint(0, 0, 1, 1, maxDepth);
  poctreenode_ = octree_.addPoint(1, 0, 1, 1, maxDepth);
  poctreenode_ = octree_.addPoint(1, 1, 1, 1, maxDepth);
  poctreenode_ = octree_.addPoint(1, 1, 1.01, 1, maxDepth);
  poctreenode_ = octree_.addPoint(0, 1, 1, 1, maxDepth);


  float x_ = 1.0f;
  float y_ = 1.0f;
  float z_ = 1.05f;
  float dr = pow(2.0,8);
  
  Eigen::Matrix<float,4,1> currentpos_(x_,y_,z_,1);
  cout << "Test position" << endl;
  cout << currentpos_[0] <<"," << currentpos_[1] << "," << currentpos_[2] << endl;
  cout << "dr = " << dr << endl;
  
  vector< OcTreeNode< float, int >* > nodes_;
  Eigen::Matrix<float,4,1> minPosition_ (x_-dr,y_-dr,z_-dr,1);
  Eigen::Matrix<float,4,1> maxPosition_ (x_+dr,y_+dr,z_+dr,1);
//  prepresentativetarget_->getAllNodesInVolumeOnDepth(nodes_,minPosition_,maxPosition_,maxDepth,true);
  //-------------------------------
  OcTreeKey<float,int> octreekey_;
  Eigen::Matrix<float,4,1> pos_;
  octreekey_ = octree_.getKey(minPosition_ );
  pos_ = octreekey_.getPosition(&octree_);
  cout << "minPosition_ " << pos_[0] << "," << pos_[1] << "," << pos_[2]<< endl;
  octreekey_ = octree_.getKey(maxPosition_ );
  pos_ = octreekey_.getPosition(&octree_);
  cout << "maxPosition_ " << pos_[0] << "," << pos_[1] << "," << pos_[2]<< endl;
  //-------------------------------
  
  octree_.getAllNodesInVolumeOnDepth(nodes_,minPosition_,maxPosition_,maxDepth,true);
  cout << "collision numbers = " << nodes_.size() << endl;
  
  if (nodes_.size() > 0) {
    cout << "-------this is collision-------" << endl;
  }else {
    cout << "-------no collision, we need add point.------" << endl;
  }

  for (unsigned int i=0; i<nodes_.size(); i++) {
    Eigen::Matrix<float,4,1> neighborpos_ = (nodes_[i])->getPosition();
    cout << nodes_[i] << "\t";
    cout << neighborpos_[0] << "," << neighborpos_[1] << "," << neighborpos_[2]  << endl;
  }


  //---------------debug---------------
    
  Eigen::Matrix< float, 4, 1 > m_center;
  m_center = Eigen::Matrix< float, 4, 1 >(0.0f,0.0f,0.0f,1);
  m_center = Eigen::Matrix< float, 4, 1 >(236.693, 995.11, 482.206, 1);
  
  float m_fCollisionLength = 100.0f;
  float m_fNucleusSphereDiameter = 6000.0f;
  float m_minimumVolumeSize = m_fCollisionLength; 
  float m_dr = m_minimumVolumeSize;
  allocator = boost::make_shared< OcTreeNodeAllocator< float , int > >();
  OcTree<float,int>* m_pOctree;
  m_pOctree = new OcTree<float,int>(m_center,m_minimumVolumeSize,m_fNucleusSphereDiameter,allocator);
  int m_maxDepth = ceil(m_pOctree->depthForVolumeSize(m_minimumVolumeSize));
  cout << "-------------------" << endl;
  cout << "minimumVolumeSize =" << m_minimumVolumeSize << endl;
  cout << "maxDistance =" << m_fNucleusSphereDiameter << endl;
  cout << "maxDepth = " <<  m_maxDepth << endl;
  cout << "mdr =" << m_dr << endl;
  
  m_pOctree->addPoint(1256.19,  1478.53, -1194.18,1 , m_maxDepth);
  x_ = 1214.87;
  y_ = 1459.49;
  z_ = -1198.86;

   minPosition_ = Eigen::Matrix<float, 4, 1> ( x_ - m_dr, y_ - m_dr, z_ - m_dr, 1);
   maxPosition_ = Eigen::Matrix<float, 4, 1> ( x_ + m_dr, y_ + m_dr, z_ + m_dr, 1);
  nodes_.clear();
  m_pOctree->getAllNodesInVolumeOnDepth(nodes_, minPosition_, maxPosition_, m_maxDepth, true);
  
  if (nodes_.size() > 0){
    bool bFlag = false;
    for (int i = 0; i< (int)nodes_.size(); i++) {
      Eigen::Matrix<float,4,1> point_ = (nodes_[i])->getPosition();
      float point_x_ = point_[0];
      float point_y_ = point_[1];
      float point_z_ = point_[2];
      float dist_ = sqrt((point_x_ - x_)*(point_x_ - x_) + (point_y_ - y_ )*(point_y_ - y_) + (point_z_ - z_)*(point_z_ - z_));
      if (dist_ < m_fCollisionLength) {
        bFlag = true;
        cout << x_ << " " << y_ << " " << z_ <<endl;
        cout << point_x_ << " " << point_y_ << " " << point_z_ <<endl;
        cout << dist_ << endl;
        cout << "collision happend";
        break;
      }
    }
    return bFlag;
  }else{
    cout << "aa" << endl;
    return false;
  }

  return 0;
  
  
}
