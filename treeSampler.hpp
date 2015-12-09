#ifndef TREE_SAMPLER_HPP
#define TREE_SAMPLER_HPP

#include <cstdlib>
#include <cstddef>
#include <vector>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <fstream>
#include <random>

#include "stepRegression.hpp"
#include "asymmTree.hpp"

template<class pointType>
class TreeSampler
{

public:

  typedef typename pointType::realScalarType RealScalarType;
  typedef asymmTree<pointType> AsymmTreeType;
  typedef std::vector<pointType> PointsArrayType;
  typedef StepRegression RegressionType;
  typedef Eigen::VectorXd VectorXd;

  TreeSampler()
  {

  }

  void setup(AsymmTreeType* arg_tree_root) {

    // Sorted list of active nodes by ascending volume
    mTreeRoot = arg_tree_root;
    list_active_nodes(mTreeRoot,mActiveNodes);
    std::sort(mActiveNodes.begin(),mActiveNodes.end(),compareVolumes);

    // Number of active nodes
    size_t nnodes = mActiveNodes.size();

    // Interpolation vectors
    VectorXd idx   = VectorXd(nnodes+1);
    VectorXd fcvol = VectorXd(nnodes+1);

    // Initial values
    idx[0]   = 0.;
    fcvol[0] = 0.;

    // Store cumulative volumes and reference indices
    for(size_t ii=0;ii<nnodes;ii++)
    {
      idx[ii+1]   = ii+1;
      fcvol[ii+1] = fcvol[ii] + nodeVolume(mActiveNodes[ii]);
    }

    // Inverse total volume
    RealScalarType icvol = 1. / fcvol[nnodes];

    // Convert to cumulative fractional volumes
    for(size_t ii=0;ii<nnodes;ii++) { fcvol[ii+1] *= icvol; }

    // Fit step regressor
    mRegressionIndex.setup(fcvol,idx);
  }


  template<class RNGType>
  pointType randomPoint(RNGType & rng)
  {
    std::uniform_real_distribution<> distUniReal;
    size_t random_node_index = mRegressionIndex(distUniReal(rng));
    assert(random_node_index >= 0 && random_node_index < mActiveNodes.size());
    return mActiveNodes[random_node_index]->walkRandomPoint(rng);
  }

private:

  // Build a standard vector of all 'active' (leaf) nodes in the tree
  void list_active_nodes(AsymmTreeType* arg_tree_root, std::vector<AsymmTreeType*> &arg_node_list)
  {
    arg_node_list.clear();
    recurse_active_nodes(arg_tree_root,arg_node_list);
  }

  // Recursively add leaf (sub)nodes to a vector of nodes
  void recurse_active_nodes(AsymmTreeType* arg_tree_node, std::vector<AsymmTreeType*> &arg_node_list)
  {
    if(!( arg_tree_node->hasLeftSubTree() || arg_tree_node->hasRightSubTree() ))
    {
      arg_node_list.push_back(arg_tree_node);
      return;
    }

    if(arg_tree_node->hasLeftSubTree())
    {
      recurse_active_nodes(arg_tree_node->leftSubTree(),arg_node_list);
    }

    if(arg_tree_node->hasRightSubTree())
    {
      recurse_active_nodes(arg_tree_node->rightSubTree(),arg_node_list);
    }
  }

  // Return the volume of the given node
  static RealScalarType nodeVolume(AsymmTreeType* arg_node)
  {
    RealScalarType volume = 1.;
    pointType boundMin;
    pointType boundMax;
    arg_node->getBounds(boundMin,boundMax,arg_node->treeIndex());

    assert(boundMin.size() == boundMax.size());
    size_t ndim = boundMin.size();

    for(size_t ii=0; ii<ndim; ii++) { volume *= (boundMax[ii]  - boundMin[ii]); }

    return volume;
  }

  // return true if the volume of the node lhs is greter than that of rhs
  static bool compareVolumes(AsymmTreeType* lhs, AsymmTreeType* rhs)
  {
    return nodeVolume(lhs) > nodeVolume(rhs);
  }

  AsymmTreeType* mTreeRoot;
  std::vector<AsymmTreeType*> mActiveNodes;
  RegressionType mRegressionIndex;

};

#endif // TREE_SAMPLER_HPP
