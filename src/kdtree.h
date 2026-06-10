#ifndef SPHSOL_KDTREE_H_
#define SPHSOL_KDTREE_H_


#include"global.h"
#include"math.h"
#include"kernel.h"
#include"particle.h"


//--------------------------------------------------------------------------------------------------
// Node strcut in k-d tree
//--------------------------------------------------------------------------------------------------
struct Node
{
  // Node(): idx(-1), axis(-1), monopole(0.0), com(Vec3d::Zero()), quadpole(Mat3d::Zero())
  Node(): idx(-1), axis(-1), monopole(0.0), com(Vec3d::Zero())
  {
    vol_min = Vec3d(1.0e99, 1.0e99, 1.0e99);
    vol_max = Vec3d(-1.0e99, -1.0e99, -1.0e99);
    next[0] = next[1] = nullptr;
    
#ifdef LEAF_SIZE
    is_leaf = false;
    leaf_idx = nullptr;
    leaf_num = 0;
#endif
  };

  int idx;        // Index to the original point
  int axis;       // Dimension's axis

#ifdef LEAF_SIZE
  bool is_leaf;
  int* leaf_idx;
  int  leaf_num;
#endif
  
  // use for self-gravity only
  double monopole;        // sum(m_i)
  Vec3d com;              // sum(m_i * pos_i)/sum(m_i)
  // Mat3d quadpole;         // sum(m_i * (pos_i-com) * (pos_i-com).transpose())
  Vec3d vol_min, vol_max; // volume limit

  Node* next[2];  // Pointers to the child nodes
};


//--------------------------------------------------------------------------------------------------
// K-d tree for neighbor searching
//--------------------------------------------------------------------------------------------------
class KdTree
{
public:

  KdTree(const Particle* particles, Kernel& kernel, double hmax, double href, int particles_num, bool calc_grav = false) :
      root_(nullptr), particles_(particles), kernel_(kernel), hmax_(hmax), href_(href), particles_num_(particles_num), calc_grav_(calc_grav)
      {
        href2_ = href_ * href_;
        href3_ = href_ * href2_;
        Build();
      };
  ~KdTree() { Clear(); };

  void Build();
  void Clear();

  void NeighborSearch(int query_idx, std::vector<Neighbor>* neighbors, int* neighbor_num, double kappa = -1.0);

  void CalcSelfGravity(int query_idx, Vec3d* grav);

private:

  Node* root_;      // Root node of the tree

  double hmax_;     // Current max hsml to cut off
  double href_, href2_, href3_;

  Kernel& kernel_;  // SPH kernel to calculate neighbors information

  const Particle* particles_;
  int particles_num_;

  bool calc_grav_; // if use multi-pole moments to calculate self gravity

  Node* BuildRecursive(int* indices, int npoints, int depth);

  void ClearRecursive(Node* node);

  void NeighborSearchRecursive(const Node* node, int query_idx, std::vector<Neighbor>* neighbors, int* neighbor_num, double kappa);

  void CalcSelfGravityRecursive(const Node* node, int query_idx, Vec3d* grav);

};

#endif

