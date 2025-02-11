/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2011, University of Bonn, Computer Science Institute VI
 *  Author: Joerg Stueckler, 4/2011
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of University of Bonn, Computer Science Institute
 *     VI nor the names of its contributors may be used to endorse or
 *     promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include <cmath>



// insert a 0 bit after each of the 16 low bits of x
template< typename CoordType, typename ValueType >
inline uint64_t spatialaggregate::OcTreeKey< CoordType, ValueType >::dilate1By2( uint64_t x ) {
  x &= 0x0000ffffLL;
  x = (x ^ (x << 16)) & 0x0000ff0000ffLL;
  x = (x ^ (x << 8))  & 0x00f00f00f00fLL;
  x = (x ^ (x << 4))  & 0x0c30c30c30c3LL;
  x = (x ^ (x << 2))  & 0x249249249249LL;
  return x;
}


// Inverse of Part1By2 - "delete" all bits not at positions divisible by 3
template< typename CoordType, typename ValueType >
inline uint64_t spatialaggregate::OcTreeKey< CoordType, ValueType >::reduce1By2( uint64_t x ) {
  x &= 0x249249249249LL;
  x = (x ^ (x >>  2)) & 0x0c30c30c30c3LL;
  x = (x ^ (x >>  4)) & 0x00f00f00f00fLL;
  x = (x ^ (x >>  8)) & 0x0000ff0000ffLL;
  x = (x ^ (x >> 16)) & 0x0000ffffLL;
  return x;
}


// generates morton code for up to 16 bit in each dimension
template< typename CoordType, typename ValueType >
inline uint64_t spatialaggregate::OcTreeKey< CoordType, ValueType >::encodeMorton48( uint64_t x, uint64_t y, uint64_t z ) {

  return (dilate1By2(x) << 2) | (dilate1By2(y) << 1) | dilate1By2(z);
}

// decodes morton code for up to 16 bit in each dimension
template< typename CoordType, typename ValueType >
inline void spatialaggregate::OcTreeKey< CoordType, ValueType >::decodeMorton48( uint64_t key, uint64_t& x, uint64_t& y, uint64_t& z ) {

  x = reduce1By2( key >> 2 );
  y = reduce1By2( key >> 1 );
  z = reduce1By2( key );

}



template< typename CoordType, typename ValueType >
inline Eigen::Matrix< CoordType, 4, 1 > spatialaggregate::OcTreeKey< CoordType, ValueType >::getPosition( OcTree< CoordType, ValueType >* tree ) const {

  // normalize to min max position
  Eigen::Matrix< CoordType, 4, 1 > pos( x_, y_, z_, 1 );

  // tree position normalizer ensures that we only use 16 bit in each dimension (implies max depth of 16)
  pos = pos.cwiseProduct( tree->inv_position_normalizer_ ).eval();
  pos += tree->min_position_;

  return pos;

}


template< typename CoordType, typename ValueType >
inline void spatialaggregate::OcTreeKey< CoordType, ValueType >::setKey( const Eigen::Matrix< CoordType, 4, 1 >& position, OcTree< CoordType, ValueType >* tree ) {

  // normalize to min max position
  Eigen::Matrix< CoordType, 4, 1 > pos = position;
  pos -= tree->min_position_;

  // tree position normalizer ensures that we only use 16 bit in each dimension (implies max depth of 16)
  pos = pos.cwiseProduct( tree->position_normalizer_ ).eval();

  x_ = pos(0) + 0.5;
  y_ = pos(1) + 0.5;
  z_ = pos(2) + 0.5;

}


template< typename CoordType, typename ValueType >
inline unsigned int spatialaggregate::OcTreeNode< CoordType, ValueType >::getOctant( const spatialaggregate::OcTreeKey< CoordType, ValueType >& query ) {

  // just the xyz shuffle at the depth of the node

  const uint32_t mask = tree_->depth_masks_[depth_];

  return ((!!(query.x_ & mask)) << 2) |
       ((!!(query.y_ & mask)) << 1) |
       (!!(query.z_ & mask));

}


template< typename CoordType, typename ValueType >
inline bool spatialaggregate::OcTreeNode< CoordType, ValueType >::inRegion( const OcTreeKey< CoordType, ValueType >& minPosition, const OcTreeKey< CoordType, ValueType >& maxPosition ) {

  return pos_key_.x_ >= minPosition.x_ && pos_key_.x_ <= maxPosition.x_ && pos_key_.y_ >= minPosition.y_ && pos_key_.y_ <= maxPosition.y_ && pos_key_.z_ >= minPosition.z_ && pos_key_.z_ <= maxPosition.z_;

//  if( pos_key_.x_ >= minPosition.x_ && pos_key_.x_ <= maxPosition.x_ && pos_key_.y_ >= minPosition.y_ && pos_key_.y_ <= maxPosition.y_ && pos_key_.z_ >= minPosition.z_ && pos_key_.z_ <= maxPosition.z_ ) {
//    return true;
//  }
//  else
//    return false;

}

template< typename CoordType, typename ValueType >
inline bool spatialaggregate::OcTreeNode< CoordType, ValueType >::inRegion( const OcTreeKey< CoordType, ValueType >& position ) {

  return position.x_ >= min_key_.x_ && position.x_ <= max_key_.x_ && position.y_ >= min_key_.y_ && position.y_ <= max_key_.y_ && position.z_ >= min_key_.z_ && position.z_ <= max_key_.z_;

//  if( position.x_ >= min_key_.x_ && position.x_ <= max_key_.x_ && position.y_ >= min_key_.y_ && position.y_ <= max_key_.y_ && position.z_ >= min_key_.z_ && position.z_ <= max_key_.z_ ) {
//    return true;
//  }
//  else
//    return false;

}

template< typename CoordType, typename ValueType >
inline bool spatialaggregate::OcTreeNode< CoordType, ValueType >::overlap( const OcTreeKey< CoordType, ValueType >& minPosition, const OcTreeKey< CoordType, ValueType >& maxPosition ) {

  return (max_key_.x_ >= minPosition.x_ && max_key_.y_ >= minPosition.y_ && max_key_.z_ >= minPosition.z_ && min_key_.x_ <= maxPosition.x_ && min_key_.y_ <= maxPosition.y_ && min_key_.z_ <= maxPosition.z_);

//  if( max_key_.x_ < minPosition.x_ || max_key_.y_ < minPosition.y_ || max_key_.z_ < minPosition.z_ || min_key_.x_ > maxPosition.x_ || min_key_.y_ > maxPosition.y_ || min_key_.z_ > maxPosition.z_ ) {
//    return false;
//  }
//  else {
//    return true;
//  }

}


template< typename CoordType, typename ValueType >
inline bool spatialaggregate::OcTreeNode< CoordType, ValueType >::containedInRegion( const OcTreeKey< CoordType, ValueType >& minPosition, const OcTreeKey< CoordType, ValueType >& maxPosition ) {

  return (overlap( minPosition, maxPosition ) && min_key_.x_ >= minPosition.x_ && min_key_.y_ >= minPosition.y_ && min_key_.z_ >= minPosition.z_ && max_key_.x_ <= maxPosition.x_ && max_key_.y_ <= maxPosition.y_ && max_key_.z_ <= maxPosition.z_);

//  if( !overlap( minPosition, maxPosition ) )
//    return false;
//
//  if( min_key_.x_ < minPosition.x_ || min_key_.y_ < minPosition.y_ || min_key_.z_ < minPosition.z_ || max_key_.x_ > maxPosition.x_ || max_key_.y_ > maxPosition.y_ || max_key_.z_ > maxPosition.z_ )
//    return false;
//
//  return true;
}


template< typename CoordType, typename ValueType >
inline bool spatialaggregate::OcTreeNode< CoordType, ValueType >::regionContained( const OcTreeKey< CoordType, ValueType >& minPosition, const OcTreeKey< CoordType, ValueType >& maxPosition ) {

  return (overlap( minPosition, maxPosition ) && min_key_.x_ >= minPosition.x_ && min_key_.y_ >= minPosition.y_ && min_key_.z_ >= minPosition.z_ && max_key_.x_ <= maxPosition.x_ && max_key_.y_ <= maxPosition.y_ && max_key_.z_ <= maxPosition.z_);

//  if( !overlap( minPosition, maxPosition ) )
//    return false;
//
//  if( min_key_.x_ < minPosition.x_ || min_key_.y_ < minPosition.y_ || min_key_.z_ < minPosition.z_ || max_key_.x_ > maxPosition.x_ || max_key_.y_ > maxPosition.y_ || max_key_.z_ > maxPosition.z_ )
//    return false;

  return true;
}

template< typename CoordType, typename ValueType >
inline void spatialaggregate::OcTreeNode< CoordType, ValueType >::sweepUp( void* data, void (*f)( spatialaggregate::OcTreeNode< CoordType, ValueType >* current, spatialaggregate::OcTreeNode< CoordType, ValueType >* next, void* data ) ) {

  f( this, parent_, data );

  if( parent_ ) {
    parent_->sweepUp( data, f );
  }

}


template< typename CoordType, typename ValueType >
inline void spatialaggregate::OcTreeNode< CoordType, ValueType >::getAllLeaves( std::list< spatialaggregate::OcTreeNode< CoordType, ValueType >* >& points ) {

  if( type_ == OCTREE_LEAF_NODE ) {

    points.push_back( this );

  }
  else {

    // for all children_
    // - if regions overlap: call add points
    for( unsigned int i = 0; i < 8; i++ ) {
      if( !children_[i] )
        continue;

      children_[i]->getAllLeaves( points );
    }

  }

}


template< typename CoordType, typename ValueType >
inline void spatialaggregate::OcTreeNode< CoordType, ValueType >::getAllLeavesInVolume( std::list< spatialaggregate::OcTreeNode< CoordType, ValueType >* >& points, const spatialaggregate::OcTreeKey< CoordType, ValueType >& minPosition, const spatialaggregate::OcTreeKey< CoordType, ValueType >& maxPosition, int maxDepth ) {

  if( type_ == OCTREE_LEAF_NODE ) {

    // check if point in leaf is within region
    if( inRegion( minPosition, maxPosition ) ) {
      points.push_back( this );
    }

  }
  else {

    if( depth_ > maxDepth ) {
      return;
    }

    // for all children_
    // - if regions overlap: call add points
    for( unsigned int i = 0; i < 8; i++ ) {
      if( !children_[i] )
        continue;

      if( children_[i]->overlap( minPosition, maxPosition ) )
        children_[i]->getAllLeavesInVolume( points, minPosition, maxPosition, maxDepth );
    }

  }

}


template< typename CoordType, typename ValueType >
inline void spatialaggregate::OcTreeNode< CoordType, ValueType >::getAllNodesInVolumeOnDepth( std::list< OcTreeNode< CoordType, ValueType >* >& points, const OcTreeKey< CoordType, ValueType >& minPosition, const OcTreeKey< CoordType, ValueType >& maxPosition, int maxDepth, bool higherDepthLeaves ) {

  if( depth_ > maxDepth )
    return;

  if( type_ == OCTREE_LEAF_NODE ) {

    if( !higherDepthLeaves && depth_ != maxDepth )
      return;

    // check if point in leaf is within region
    if( inRegion( minPosition, maxPosition ) ) {
      points.push_back( this );
      return;
    }

  }
  else {

    if( depth_ == maxDepth ) {
      if( inRegion( minPosition, maxPosition ) )
        points.push_back( this );
      return;
    }

    // for all children
    // - if regions overlap: call function for the child
    for( unsigned int i = 0; i < 8; i++ ) {
      if( !children_[i] )
        continue;

      if( children_[i]->overlap( minPosition, maxPosition ) )
        children_[i]->getAllNodesInVolumeOnDepth( points, minPosition, maxPosition, maxDepth, higherDepthLeaves );
    }
  }

}


template< typename CoordType, typename ValueType >
inline void spatialaggregate::OcTreeNode< CoordType, ValueType >::getAllNodesInVolumeOnDepth( std::vector< OcTreeNode< CoordType, ValueType >* >& points, const OcTreeKey< CoordType, ValueType >& minPosition, const OcTreeKey< CoordType, ValueType >& maxPosition, int maxDepth, bool higherDepthLeaves ) {

//  if( depth_ == maxDepth ) {
//    if( inRegion( minPosition, maxPosition ) ) {
//      points.push_back( this );
//    }
//  }
//  else if( depth_ < maxDepth ) {
//
//    // for all children
//    // - if regions overlap: call function for the child
//    for( unsigned int i = 0; i < 8; i++ ) {
//      if( children_[i] && children_[i]->overlap( minPosition, maxPosition ) )
//        children_[i]->getAllNodesInVolumeOnDepth( points, minPosition, maxPosition, maxDepth, higherDepthLeaves );
//    }
//
//  }


  if( depth_ > maxDepth )
    return;

  if( type_ == OCTREE_LEAF_NODE ) {

    if( !higherDepthLeaves && depth_ != maxDepth )
      return;

    // check if point in leaf is within region
    if( inRegion( minPosition, maxPosition ) ) {
      points.push_back( this );
      return;
    }

  }
  else {

    if( depth_ == maxDepth ) {
      if( inRegion( minPosition, maxPosition ) )
        points.push_back( this );
      return;
    }

    // for all children
    // - if regions overlap: call function for the child
    for( unsigned int i = 0; i < 8; i++ ) {
      if( !children_[i] )
        continue;

      if( children_[i]->overlap( minPosition, maxPosition ) )
        children_[i]->getAllNodesInVolumeOnDepth( points, minPosition, maxPosition, maxDepth, higherDepthLeaves );
    }
  }

}


template< typename CoordType, typename ValueType >
inline void spatialaggregate::OcTreeNode< CoordType, ValueType >::getNeighbors( std::list< OcTreeNode< CoordType, ValueType >* >& neighbors ) {

  for( int n = 0; n < 27; n++ ) {
    if( neighbors_[n] ) {
      neighbors.push_back( neighbors_[n] );
    }
  }

}


template< typename CoordType, typename ValueType >
inline void spatialaggregate::OcTreeNode< CoordType, ValueType >::getNeighbors( std::vector< OcTreeNode< CoordType, ValueType >* >& neighbors ) {

  neighbors.assign( &neighbors_[0], &neighbors_[26] );

//  neighbors.resize(27);
//  for( int n = 0; n < 27; n++ ) {
//    if( neighbors_[n] ) {
//      neighbors[i++] = neighbors_[n];
//    }
//    neighbors[n] = neighbors_[n];
//  }
//  neighbors.resize(i);

}


template< typename CoordType, typename ValueType >
inline ValueType spatialaggregate::OcTreeNode< CoordType, ValueType >::getValueInVolume( const spatialaggregate::OcTreeKey< CoordType, ValueType >& minPosition, const spatialaggregate::OcTreeKey< CoordType, ValueType >& maxPosition, int maxDepth ) {

  if( type_ == OCTREE_LEAF_NODE ) {

    if( inRegion( minPosition, maxPosition ) )
      return value_;

    return ValueType(0);

  }
  else {

    if( !overlap( minPosition, maxPosition ) )
      return ValueType(0);

    if( containedInRegion( minPosition, maxPosition ) )
      return value_;

    if( depth_ >= maxDepth ) {
      return value_;
    }

    ValueType value = ValueType(0);
    for( unsigned int i = 0; i < 8; i++ ) {
      if(!children_[i])
        continue;
      value += children_[i]->getValueInVolume( minPosition, maxPosition, maxDepth );
    }

    return value;

  }

}




template< typename CoordType, typename ValueType >
inline void spatialaggregate::OcTreeNode< CoordType, ValueType >::applyOperatorInVolume( ValueType& value, void* data, void (*f)( ValueType& v, spatialaggregate::OcTreeNode< CoordType, ValueType >* current, void* data ), const spatialaggregate::OcTreeKey< CoordType, ValueType >& minPosition, const spatialaggregate::OcTreeKey< CoordType, ValueType >& maxPosition, int maxDepth ) {

  if( type_ == OCTREE_LEAF_NODE ) {

    if( inRegion( minPosition, maxPosition ) ) {
      f( value, this, data );
    }

  }
  else {

    if( !overlap( minPosition, maxPosition ) )
      return;

    if( containedInRegion( minPosition, maxPosition ) ) {
      f( value, this, data );
      return;
    }

    if( depth_ >= maxDepth ) {
      // since we check for overlap above, this branching node is only accounted for, when its extent overlaps with the search region
      f( value, this, data );
      return;
    }

    for( unsigned int i = 0; i < 8; i++ ) {
      if(!children_[i])
        continue;

      children_[i]->applyOperatorInVolume( value, data, f, minPosition, maxPosition, maxDepth );

    }

  }

}




template< typename CoordType, typename ValueType >
inline void spatialaggregate::OcTreeNode< CoordType, ValueType >::sweepDown( void* data, void (*f)( spatialaggregate::OcTreeNode< CoordType, ValueType >* current, spatialaggregate::OcTreeNode< CoordType, ValueType >* next, void* data ) ) {

  if( type_ == OCTREE_LEAF_NODE ) {

    f( this, NULL, data );

  }
  else {

    f( this, NULL, data );

    for( unsigned int i = 0; i < 8; i++ ) {
      if( children_[i] )
        children_[i]->sweepDown( data, f );
    }

  }

}

template< typename CoordType, typename ValueType >
inline void spatialaggregate::OcTreeNode< CoordType, ValueType >::initialize( OcTreeNodeType type, const OcTreeKey< CoordType, ValueType >& key, const ValueType& value, int depth, OcTreeNode< CoordType, ValueType >* parent, OcTree< CoordType, ValueType >* tree ) {
  type_ = type;
  value_ = value;
  depth_ = depth;
  parent_ = parent;
  tree_ = tree;
  pos_key_ = key;

  const uint32_t minmask = tree_->minmasks_[depth_];
  const uint32_t maxmask = tree_->maxmasks_[depth_];

  min_key_.x_ = pos_key_.x_ & minmask;
  min_key_.y_ = pos_key_.y_ & minmask;
  min_key_.z_ = pos_key_.z_ & minmask;

  max_key_.x_ = pos_key_.x_ | maxmask;
  max_key_.y_ = pos_key_.y_ | maxmask;
  max_key_.z_ = pos_key_.z_ | maxmask;

  memset( children_, 0, sizeof( children_ ) );
  memset( neighbors_, 0, sizeof( neighbors_ ) );
}


template< typename CoordType, typename ValueType >
inline void spatialaggregate::OcTreeNode< CoordType, ValueType >::getCenterKey( OcTreeKey< CoordType, ValueType >& center_key ) const {

  const uint32_t center_diff =  (max_key_.x_ - min_key_.x_) >> 1;
  center_key.x_ = min_key_.x_ + center_diff;
  center_key.y_ = min_key_.y_ + center_diff;
  center_key.z_ = min_key_.z_ + center_diff;

}



template< typename CoordType, typename ValueType >
inline spatialaggregate::OcTreeNode< CoordType, ValueType >* spatialaggregate::OcTreeNode< CoordType, ValueType >::addPoint( const OcTreeKey< CoordType, ValueType >& position, const ValueType& value, int maxDepth ) {

//  Eigen::Matrix< CoordType, 4, 1 > p = position.getPosition( tree_ );
//
//  std::cout << p << "\n";

  // traverse from root until we found an empty leaf node
  spatialaggregate::OcTreeNode< CoordType, ValueType >* parentNode = parent_;
  unsigned int octant = 0;
  spatialaggregate::OcTreeNode< CoordType, ValueType >* currNode = this;

  while( currNode ) {

    if( currNode->depth_ == maxDepth ) {
      // reached max depth
      currNode->value_ += value;
      return currNode;
    }

    if( currNode->type_ == OCTREE_LEAF_NODE )
      break; // reached a leaf, stop searching

    parentNode = currNode;
    octant = currNode->getOctant( position );
    currNode = currNode->children_[ octant ];

    // add value to integral value of parent node
    parentNode->value_ += value;

  }


  if( currNode == NULL ) {

    // simply add a new leaf node in the parent's octant
    spatialaggregate::OcTreeNode< CoordType, ValueType >* leaf = tree_->allocator_->allocateNode();
    leaf->initialize( OCTREE_LEAF_NODE, position, value, parentNode->depth_ + 1, parentNode, tree_ );

    parentNode->type_ = OCTREE_BRANCHING_NODE;
    parentNode->children_[octant] = leaf;

//    parentNode->finishBranch();
    
    return leaf;

  }
  else {

    if( currNode->type_ == OCTREE_LEAF_NODE ) {

      // branch at parent's octant..
      spatialaggregate::OcTreeNode< CoordType, ValueType >* oldLeaf = currNode;

      if( oldLeaf->pos_key_ == position ) {
        oldLeaf->value_ += value;
        return oldLeaf;
      }

      spatialaggregate::OcTreeNode< CoordType, ValueType >* branch = tree_->allocator_->allocateNode();
      branch->initialize( oldLeaf );
      branch->type_ = OCTREE_BRANCHING_NODE;

      assert( parentNode->children_[octant] == oldLeaf );

      parentNode->children_[octant] = branch;

      // link old leaf to branching node
      unsigned int oldLeafOctant = branch->getOctant( oldLeaf->pos_key_ );
      branch->children_[ oldLeafOctant ] = oldLeaf;
      oldLeaf->initialize( OCTREE_LEAF_NODE, oldLeaf->pos_key_, oldLeaf->value_, branch->depth_ + 1, branch, tree_ );

//      branch->finishBranch();
        
      return branch->addPoint( position, value, maxDepth ); // this could return some older leaf

    }

  }

  assert( false );

  return NULL;

}




template< typename CoordType, typename ValueType >
inline spatialaggregate::OcTreeNode< CoordType, ValueType >* spatialaggregate::OcTreeNode< CoordType, ValueType >::findRepresentative( const OcTreeKey< CoordType, ValueType >& position, int maxDepth ) {
  if( type_ == OCTREE_LEAF_NODE )
    return this;
  else {

    if( depth_ == maxDepth )
      return this;

    spatialaggregate::OcTreeNode< CoordType, ValueType >* n = children_[ getOctant(position) ];
    if( n )
      return n->findRepresentative( position, maxDepth );
    else
      return this;
  }
}



template< typename CoordType, typename ValueType >
inline spatialaggregate::OcTreeNode< CoordType, ValueType >* spatialaggregate::OcTreeNode< CoordType, ValueType >::findClosestNode( const OcTreeKey< CoordType, ValueType >& position, int depth, int& dist ) {

  if( depth_ == depth ) {

    // get distance towards this node
    dist = (pos_key_.x_-position.x_)*(pos_key_.x_-position.x_)+(pos_key_.y_-position.y_)*(pos_key_.y_-position.y_)+(pos_key_.z_-position.z_)*(pos_key_.z_-position.z_);
    return this;

  }
  else {

    spatialaggregate::OcTreeNode< CoordType, ValueType >* closest = NULL;
    int bestDist = std::numeric_limits<int>::max();

    for( int i = 0; i < 8; i++ ) {
      if( children_[i] ) {
        int childDist = 0;
        spatialaggregate::OcTreeNode< CoordType, ValueType >* childClosest = children_[i]->findClosestNode( position, depth, childDist );
        if( childDist < bestDist ) {
          bestDist = childDist;
          closest = childClosest;
        }
      }
    }

    dist = bestDist;
    return closest;

  }

}


template< typename CoordType, typename ValueType >
inline CoordType spatialaggregate::OcTreeNode< CoordType, ValueType >::resolution() {

  return tree_->volumeSizeForDepth( depth_ );

}

template< typename CoordType, typename ValueType >
inline CoordType spatialaggregate::OcTreeNode< CoordType, ValueType >::invResolution() {

  return tree_->inv_resolutions_[ depth_ ];

}

template< typename CoordType, typename ValueType >
inline void spatialaggregate::OcTreeNode< CoordType, ValueType >::finishBranch() {

  const uint32_t mask = tree_->depth_masks_[depth_];

  for( unsigned int i = 0; i < 8; i++ ) {
    if( !children_[i] ) {
      children_[i] = tree_->allocator_->allocateNode();

      OcTreeKey< CoordType, ValueType > key;

      if( i & 0x4 ) // x set?
        key.x_ = pos_key_.x_ | mask;
      else
        key.x_ = pos_key_.x_ & (~mask);

      if( i & 0x2 ) // y set?
        key.y_ = pos_key_.y_ | mask;
      else
        key.y_ = pos_key_.y_ & (~mask);

      if( i & 0x1 ) // z set?
        key.z_ = pos_key_.z_ | mask;
      else
        key.z_ = pos_key_.z_ & (~mask);

      children_[i]->initialize( OCTREE_LEAF_NODE, key, 0, depth_+1, this, tree_ );
    }
  }

}


// uses precomputed look-up tables
template< typename CoordType, typename ValueType >
inline void spatialaggregate::OcTreeNode< CoordType, ValueType >::establishNeighbors() {

  neighbors_[9+3+1] = this;

  if( parent_ ) {

    const uint32_t octant = parent_->getOctant(pos_key_);

    for( unsigned int n = 0; n < 27; n++ ) {

      // set pointers
      spatialaggregate::OcTreeNode< CoordType, ValueType >* parentNeighbor = parent_->neighbors_[tree_->parent_neighbor_[octant][n]];
      if( !parentNeighbor )
        neighbors_[n] = NULL;
      else {
        neighbors_[n] = parentNeighbor->children_[tree_->neighbor_octant_[octant][n]];
      }

    }

  }

  for( unsigned int i = 0; i < 8; i++ ) {
    if( children_[i] )
      children_[i]->establishNeighbors();
  }

}


template< typename CoordType, typename ValueType >
inline spatialaggregate::OcTreeNode< CoordType, ValueType >* spatialaggregate::OcTreeNode< CoordType, ValueType >::getNeighbor( int dx, int dy, int dz ) {

  if( parent_ ) {

    const uint32_t octant = parent_->getOctant(pos_key_);

    const int32_t n = 9*(dx+1) + 3*(dy+1) + (dz+1);

    // set pointers
    spatialaggregate::OcTreeNode< CoordType, ValueType >* parentNeighbor = parent_->neighbors_[tree_->parent_neighbor_[octant][n]];
    if( !parentNeighbor )
      return NULL;
    else {
      return parentNeighbor->children_[tree_->neighbor_octant_[octant][n]];
    }

  }

  return NULL;

}


template< typename CoordType, typename ValueType >
inline bool spatialaggregate::OcTreeNode< CoordType, ValueType >::interpolateTriLinear( double& value, spatialaggregate::OcTreeNode< CoordType, ValueType >* node, const spatialaggregate::OcTreeKey< CoordType, ValueType >& queryKey, double (*f)( spatialaggregate::OcTreeNode< CoordType, ValueType >* n ) ) {

  // determine sign of xyz direction towards query point,
  // check for 8 closest neighbors,
  // if not all neighbors exist, return false

  value = 0.0;

  spatialaggregate::OcTreeKey< CoordType, ValueType > centerKey;
  node->getCenterKey( centerKey );

  int sgnX = 1;
  if( queryKey.x_ < centerKey.x_ )
    sgnX = -1;

  int sgnY = 1;
  if( queryKey.y_ < centerKey.y_ )
    sgnY = -1;

  int sgnZ = 1;
  if( queryKey.z_ < centerKey.z_ )
    sgnZ = -1;

  std::vector< spatialaggregate::OcTreeNode< CoordType, ValueType >* > neighborList( 8 );

  neighborList[0] = node;
  neighborList[1] = node->getNeighbor( 0, 0, sgnZ );
  neighborList[2] = node->getNeighbor( 0, sgnY, 0 );
  neighborList[3] = node->getNeighbor( 0, sgnY, sgnZ );
  neighborList[4] = node->getNeighbor( sgnX, 0, 0 );
  neighborList[5] = node->getNeighbor( sgnX, 0, sgnZ );
  neighborList[6] = node->getNeighbor( sgnX, sgnY, 0 );
  neighborList[7] = node->getNeighbor( sgnX, sgnY, sgnZ );

  bool allNodes = true;
  for( int i = 0; i < 8; i++ )
    if(!neighborList[i])
      allNodes = false;

  if( allNodes ) {

    // tri-linear interpolation
    // iterate through the node list, determine distances in x, y, z from the center keys
    // normalize to resolution and convert to double, weight with (1-diffx)*(1-diffy)*(1-diffz)
    // return weighted average

    double sumWeight = 0.0;

    for( int i = 0; i < 8; i++ ) {

      spatialaggregate::OcTreeKey< CoordType, ValueType > centerKeyN;
      neighborList[i]->getCenterKey( centerKeyN );

      CoordType dx = node->resolution() - tree_->inv_position_normalizer_(0) * (CoordType)abs((int)centerKeyN.x_ - (int)queryKey.x_);
      CoordType dy = node->resolution() - tree_->inv_position_normalizer_(1) * (CoordType)abs((int)centerKeyN.y_ - (int)queryKey.y_);
      CoordType dz = node->resolution() - tree_->inv_position_normalizer_(2) * (CoordType)abs((int)centerKeyN.z_ - (int)queryKey.z_);


      if( dx > 0 && dy > 0 && dz > 0 ) {

        const double weight = dx*dy*dz;
        value += weight * f( neighborList[i] );
        sumWeight += weight;

      }

    }

    if( sumWeight > std::numeric_limits<double>::epsilon() ) {
      value /= sumWeight;
      return true;
    }

  }

  return false;

}


template< typename CoordType, typename ValueType >
inline double spatialaggregate::OcTreeNode< CoordType, ValueType >::getFiniteForwardDifference( int dim, double (*f)( spatialaggregate::OcTreeNode< CoordType, ValueType >* n ) ) {

  // select next neighbor in query dimension
  double vthis = f( this );
  spatialaggregate::OcTreeNode< CoordType, ValueType >* neighborp1 = getNeighbor( dim == 0 ? 1 : 0, dim == 1 ? 1 : 0, dim == 2 ? 1 : 0 );

  if( neighborp1 ) {
//    return invResolution() * (f( neighbors_[ip1] ) - f( this ));
    return (f( neighborp1 ) - vthis);
  }
  else {

//    return 0.0;

    // take average value on positive side in dim
    double v = 0.0;
    double sum = 0.0;
    for( int di = -1; di <= 1; di++ ) {
      for( int dj = -1; dj <= 1; dj++ ) {

        spatialaggregate::OcTreeNode< CoordType, ValueType >* neighbor = NULL;
        const double weight = (1.f-0.5f*fabsf(di)) * (1.f-0.5f*fabsf(dj));

        if( dim == 0 )
          neighbor = getNeighbor( 1, di, dj );
        else if( dim == 1 )
          neighbor = getNeighbor( di, 1, dj );
        else
          neighbor = getNeighbor( di, dj, 1 );

        if( neighbor ) {
          v += weight * f( neighbor );
          sum += weight;
        }
        else {
          v += weight * vthis;
          sum += weight;
        }

      }
    }
    if( sum > 0.0 )
      v /= sum;
    else
      return 0.0;

    return (v - vthis);

  }

  return 0.0;

}


template< typename CoordType, typename ValueType >
inline double spatialaggregate::OcTreeNode< CoordType, ValueType >::getFiniteBackwardDifference( int dim, double (*f)( spatialaggregate::OcTreeNode< CoordType, ValueType >* n ) ) {

  // select next neighbor in query dimension
  double vthis = f( this );
//  int im1 = tree_->neighborhood_m1_map_[dim];
  spatialaggregate::OcTreeNode< CoordType, ValueType >* neighborm1 = getNeighbor( dim == 0 ? -1 : 0, dim == 1 ? -1 : 0, dim == 2 ? -1 : 0 );
  if( neighborm1 ) {
//    return invResolution() * (f( this ) - f( neighbors_[im1] ));
    return (vthis - f( neighborm1 ));
  }
  else {

//    return 0.0;

    // take average value on negative side in dim
    double v = 0.0;
    double sum = 0.0;
    for( int di = -1; di <= 1; di++ ) {
      for( int dj = -1; dj <= 1; dj++ ) {

        spatialaggregate::OcTreeNode< CoordType, ValueType >* neighbor = NULL;
        const double weight = (1.f-0.5f*fabsf(di)) * (1.f-0.5f*fabsf(dj));

        if( dim == 0 )
          neighbor = getNeighbor( -1, di, dj );
        else if( dim == 1 )
          neighbor = getNeighbor( di, -1, dj );
        else
          neighbor = getNeighbor( di, dj, -1 );

        if( neighbor ) {
          v += weight * f( neighbor );
          sum += weight;
        }
        else {
          v += weight * vthis;
          sum += weight;
        }

      }
    }
    if( sum > 0.0 )
      v /= sum;
    else
      return 0.0;

    return (vthis - v);

  }

  return 0.0;

}


template< typename CoordType, typename ValueType >
inline double spatialaggregate::OcTreeNode< CoordType, ValueType >::getFiniteCentralDifference( int dim, double (*f)( spatialaggregate::OcTreeNode< CoordType, ValueType >* n ) ) {

//  double v = f( this );
//
//  double vp1 = 0.0;
//  double sump1 = 0.0;
//
//  double vm1 = 0.0;
//  double summ1 = 0.0;
//
//  for( int di = -1; di <= 1; di++ ) {
//    for( int dj = -1; dj <= 1; dj++ ) {
//
//      spatialaggregate::OcTreeNode< CoordType, ValueType >* neighborp1 = NULL;
//      spatialaggregate::OcTreeNode< CoordType, ValueType >* neighborm1 = NULL;
//
//      const double weight = (1.f-0.5f*fabsf(di)) * (1.f-0.5f*fabsf(dj));
//
//      if( dim == 0 ) {
//        neighborm1 = getNeighbor( -1, di, dj );
//        neighborp1 = getNeighbor( +1, di, dj );
//      }
//      else if( dim == 1 ) {
//        neighborm1 = getNeighbor( di, -1, dj );
//        neighborp1 = getNeighbor( di, +1, dj );
//      }
//      else {
//        neighborm1 = getNeighbor( di, dj, -1 );
//        neighborp1 = getNeighbor( di, dj, +1 );
//      }
//
//      if( neighborm1 ) {
//        vm1 += weight * f( neighborm1 );
//        summ1 += weight;
//      }
//      else {
//        vm1 += weight * v;
//        summ1 += weight;
//      }
//
//      if( neighborp1 ) {
//        vp1 += weight * f( neighborp1 );
//        sump1 += weight;
//      }
//      else {
//        vp1 += weight * v;
//        sump1 += weight;
//      }
//
//    }
//  }
//
//  if( sump1 > 0 )
//    vp1 /= sump1;
//  else
//    vp1 = v;
//
//  if( summ1 > 0 )
//    vm1 /= summ1;
//  else
//    vm1 = v;
//
//  return 0.5 * (vp1-vm1);


  // select next neighbor in query dimension
  spatialaggregate::OcTreeNode< CoordType, ValueType >* neighborm1 = getNeighbor( dim == 0 ? -1 : 0, dim == 1 ? -1 : 0, dim == 2 ? -1 : 0 );
  spatialaggregate::OcTreeNode< CoordType, ValueType >* neighborp1 = getNeighbor( dim == 0 ? 1 : 0, dim == 1 ? 1 : 0, dim == 2 ? 1 : 0 );

//  double vp1 = 0.0;
//  double vm1 = 0.0;
//  if( neighborm1 )
//    vm1 = f( neighborm1 );
//  if( neighborp1 )
//    vp1 = f( neighborp1 );
//
//  return 0.5 * (vp1-vm1);


  if( neighborm1 && neighborp1 ) {
    return 0.5 * (f( neighborp1 ) - f( neighborm1 ));
  }
  else {

    double v = f( this );

    double vp1 = 0.0;

    if( neighborp1 ) {
      vp1 = f( neighborp1 );
    }
    else {

//      vp1 = v;

      double sump1 = 0.0;
      for( int di = -1; di <= 1; di++ ) {
        for( int dj = -1; dj <= 1; dj++ ) {

          spatialaggregate::OcTreeNode< CoordType, ValueType >* neighbor = NULL;

          const double weight = (1.f-0.5f*fabsf(di)) * (1.f-0.5f*fabsf(dj));

          if( dim == 0 )
            neighbor = getNeighbor( 1, di, dj );
          else if( dim == 1 )
            neighbor = getNeighbor( di, 1, dj );
          else
            neighbor = getNeighbor( di, dj, 1 );

          if( neighbor ) {
            vp1 += weight * f( neighbor );
            sump1 += weight;
          }
          else {
            vp1 += weight * v;
            sump1 += weight;
          }

        }
      }

      if( sump1 > 0 )
        vp1 /= sump1;
      else
        vp1 = v;

    }

    double vm1 = 0.0;

    if( neighborm1 ) {
      vm1 = f( neighborm1 );
    }
    else {

//      vm1 = v;

      double summ1 = 0.0;
      for( int di = -1; di <= 1; di++ ) {
        for( int dj = -1; dj <= 1; dj++ ) {

          spatialaggregate::OcTreeNode< CoordType, ValueType >* neighbor = NULL;

          const double weight = (1.f-0.5f*fabsf(di)) * (1.f-0.5f*fabsf(dj));

          if( dim == 0 )
            neighbor = getNeighbor( -1, di, dj );
          else if( dim == 1 )
            neighbor = getNeighbor( di, -1, dj );
          else
            neighbor = getNeighbor( di, dj, -1 );

          if( neighbor ) {
            vm1 += weight * f( neighbor );
            summ1 += weight;
          }
          else {
            vm1 += weight * v;
            summ1 += weight;
          }

        }
      }

      if( summ1 > 0 )
        vm1 /= summ1;
      else
        vm1 = v;

    }

    return 0.5 * ( vp1 - vm1 );
  }

  return 0.0;

}



template< typename CoordType, typename ValueType >
spatialaggregate::OcTree< CoordType, ValueType >::OcTree( const Eigen::Matrix< CoordType, 4, 1 >& center, CoordType minimumVolumeSize, CoordType maxDistance, boost::shared_ptr< spatialaggregate::OcTreeNodeAllocator< CoordType, ValueType > > allocator ) {

  allocator_ = allocator;

  // determine dimensions with 2^k * minimumVolumeSize..
  max_depth_ = ceil( log( 2.0 * maxDistance / minimumVolumeSize ) / log( 2.0 ) );
  const CoordType size = minimumVolumeSize * pow( 2.0, max_depth_ );

  assert( max_depth_ <= MAX_REPRESENTABLE_DEPTH );

  const CoordType minRepresentedVolumeSize = size * pow( 2.0, -MAX_REPRESENTABLE_DEPTH );
//  const CoordType size = minimumVolumeSize * pow( 2.0, max_depth_ );

  min_position_(0) = -0.5*size;
  min_position_(1) = -0.5*size;
  min_position_(2) = -0.5*size;
  min_position_ += center;
  min_position_(3) = 1.0;

  position_normalizer_.setConstant( 1.0 / minRepresentedVolumeSize );
  position_normalizer_(3) = 1.0;

  inv_position_normalizer_.setConstant( minRepresentedVolumeSize );
  inv_position_normalizer_(3) = 1.0;

  root_ = allocator->allocateNode();

  ValueType value;
  root_->initialize( OCTREE_BRANCHING_NODE, OcTreeKey< CoordType, ValueType >( center, this ), value, 0, NULL, this );
//  root_->finishBranch();

  minimum_volume_size_ = minimumVolumeSize;
  inv_minimum_volume_size_ = 1.0 / minimum_volume_size_;

  for( int i = 0; i <= max_depth_; i++ ) {
    resolutions_[i] = minimum_volume_size_ * pow( 2.0, max_depth_ - i );
    inv_resolutions_[i] = 1.0 / resolutions_[i];
    minResolutions_[i] = minimum_volume_size_ * pow( 2.0, (double)(max_depth_ - i) - 0.5 );
    maxResolutions_[i] = minimum_volume_size_ * pow( 2.0, (double)(max_depth_ - i) + 0.5 );
    if( i < max_depth_ )
      depth_masks_[i] = 0x1 << (MAX_REPRESENTABLE_DEPTH-i-1);
    else
      depth_masks_[i] = 0;

    minmasks_[i] = 0xFFFFFFFFLL << (MAX_REPRESENTABLE_DEPTH - i);
    maxmasks_[i] = (0x1LL << (MAX_REPRESENTABLE_DEPTH - i)) - 1LL;

  }

  log_minimum_volume_size_ = log2( minimumVolumeSize );
  log2_inv_ = 1.0 / log( 2.0 );

  neighborhood_p1_map_[0] = 22;
  neighborhood_p1_map_[1] = 16;
  neighborhood_p1_map_[2] = 14;

  neighborhood_m1_map_[0] = 4;
  neighborhood_m1_map_[1] = 10;
  neighborhood_m1_map_[2] = 12;


//  // precompute scale-depth table
//  for( unsigned int i = 0; i < 65536; i++ ) {
//
//    if( i == 0 ) {
//      scale_depth_table_[i] = max_depth_;
//      continue;
//    }
//
//
//    double depth = log2f( (double)i );
//    scale_depth_table_[i] = (double)max_depth_ - std::min( (double)max_depth_, depth );
//  }


  // precompute neighborhood look-up tables ( 8 times 27, two tables of indices for parent neighbor and its child )
  for( unsigned int i = 0; i < 8; i++ ) {

    int32_t lx = !!(i & 4);
    int32_t ly = !!(i & 2);
    int32_t lz = !!(i & 1);

    for( int dx = -1; dx <= 1; dx++ ) {
      for( int dy = -1; dy <= 1; dy++ ) {
        for( int dz = -1; dz <= 1; dz++ ) {

          int32_t nx = lx+dx;
          int32_t ny = ly+dy;
          int32_t nz = lz+dz;

          int32_t px = 0;
          int32_t py = 0;
          int32_t pz = 0;

          if( nx > 1 )
            px = 1, nx = 0;
          if( nx < 0 )
            px = -1, nx = 1;

          if( ny > 1 )
            py = 1, ny = 0;
          if( ny < 0 )
            py = -1, ny = 1;

          if( nz > 1 )
            pz = 1, nz = 0;
          if( nz < 0 )
            pz = -1, nz = 1;

          neighbor_octant_[i][9*(dx+1)+3*(dy+1)+(dz+1)] = (nx << 2) | (ny << 1) | nz;
          parent_neighbor_[i][9*(dx+1)+3*(dy+1)+(dz+1)] = 9*(px+1)+3*(py+1)+(pz+1);

        }
      }
    }
  }

}


template< typename CoordType, typename ValueType >
spatialaggregate::OcTree< CoordType, ValueType >::~OcTree() {

  if( root_ ) {
    allocator_->deallocateNode( root_ );
  }

}




