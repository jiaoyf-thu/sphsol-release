#include"kdtree.h"


//--------------------------------------------------------------------------------------------------
// Build a k-d tree
//--------------------------------------------------------------------------------------------------
void KdTree::Build()
{
  Clear();

  // Create a indices list [0 -> particless_num - 1]
  std::vector<int> indices(particles_num_);
  std::iota(std::begin(indices), std::end(indices), 0);

  // Build k-d tree recursively
  root_ = BuildRecursive(indices.data(), particles_num_, 0);
}


//--------------------------------------------------------------------------------------------------
// Clear all pointers in k-d tree
//--------------------------------------------------------------------------------------------------
void KdTree::Clear()
{
  ClearRecursive(root_);

  root_ = nullptr;
}


//--------------------------------------------------------------------------------------------------
// SPH neighbor search with k-d tree
// Default kappa = kernel.kappa, user-defined kappa is also supported
//--------------------------------------------------------------------------------------------------
void KdTree::NeighborSearch(int query_idx, std::vector<Neighbor>* neighbors, int* neighbor_num, double kappa)
{
  if (particles_[query_idx].iflag < 0) // Allow flag = 0
  {
    return;
  }

  *neighbor_num = 0;
  double my_kappa = kappa < 0.0 ? kernel_.get_kappa() : kappa;
  NeighborSearchRecursive(root_, query_idx, neighbors, neighbor_num, my_kappa);
};


//--------------------------------------------------------------------------------------------------
// Calculate self-gravity with k-d tree
//--------------------------------------------------------------------------------------------------
void KdTree::CalcSelfGravity(int query_idx, Vec3d* grav)
{
  if (!calc_grav_ || particles_[query_idx].iflag < 0)
  {
    return;
  }

  CalcSelfGravityRecursive(root_, query_idx, grav);
}


//--------------------------------------------------------------------------------------------------
// Build k-d tree recursively
//--------------------------------------------------------------------------------------------------
Node* KdTree::BuildRecursive(int* indices, int npoints, int depth)
{
  if (npoints <= 0)
    return nullptr;

#ifdef LEAF_SIZE
  if (npoints <= LEAF_SIZE)
  {
    // leaf node
    Node* node = new Node();
    node->is_leaf  = true;
    node->leaf_num = npoints;
    node->leaf_idx = new int[npoints];
    std::copy(indices, indices + npoints, node->leaf_idx);

    if (calc_grav_)
    {
      for (int i = 0; i < npoints; i++)
      {
        if (particles_[indices[i]].iflag < 0)
        {
          continue;
        }
        node->monopole += particles_[indices[i]].mass;
        node->com += particles_[indices[i]].mass * particles_[indices[i]].pos;
        node->vol_min = node->vol_min.cwiseMin(particles_[indices[i]].pos);
        node->vol_max = node->vol_max.cwiseMax(particles_[indices[i]].pos);
      }
      node->com /= node->monopole;
      // for (int i = 0; i < npoints; i++)
      // {
      //   if (particles_[indices[i]].iflag < 0)
      //   {
      //     continue;
      //   }
      //   Vec3d r = particles_[indices[i]].pos - node->com;
      //   node->quadpole += particles_[indices[i]].mass * (3.0 * r * r.transpose() - r.dot(r) * Mat3d::Identity());
      // }
    }

    return node;
  }
#endif

  const int axis = depth % 3;
  const int mid = (npoints - 1) / 2;

  // Search the (mid + 1) value in [indices, indices + npoints]
  // The complexity is O(N)
  std::nth_element(indices, indices + mid, indices + npoints, [&](int lhs, int rhs)
    {
      return particles_[lhs].pos(axis) < particles_[rhs].pos(axis);
    });

  Node* node = new Node();
  node->idx = indices[mid];
  node->axis = axis;
  
  // for the node point itself
  if (calc_grav_ && particles_[node->idx].iflag >= 0)
  {
    node->monopole = particles_[node->idx].mass; // no gravity for iflag < 0
    node->com      = particles_[node->idx].pos;
    node->vol_min = particles_[node->idx].pos;
    node->vol_max = particles_[node->idx].pos;
  }

  // child nodes
  #pragma omp task shared(npoints) firstprivate(indices, mid, depth)
  {
    node->next[0] = BuildRecursive(indices, mid, depth + 1);
  }
  #pragma omp task shared(npoints) firstprivate(indices, mid, depth)
  {
    node->next[1] = BuildRecursive(indices + mid + 1, npoints - mid - 1, depth + 1);
  }
  #pragma omp taskwait
  
  // sumup multi-pole moments & volume limits of child nodes
  if (calc_grav_)
  {
    if (node->next[0])
    {
      double m1 = node->monopole;
      double m2 = node->next[0]->monopole;
      node->monopole = m1 + m2;
      node->vol_min  = node->vol_min.cwiseMin(node->next[0]->vol_min);
      node->vol_max  = node->vol_max.cwiseMax(node->next[0]->vol_max);
      Vec3d r1  = node->com;
      Vec3d r2  = node->next[0]->com;
      node->com      = (m1 * r1 + m2 * r2) / (m1 + m2);
      // r1 -= node->com;
      // r2 -= node->com;
      // node->quadpole = node->quadpole + node->next[0]->quadpole + m1 * (3.0 * r1 * r1.transpose() - r1.dot(r1) * Mat3d::Identity())
      //                 + m2 * (3.0 * r2 * r2.transpose() - r2.dot(r2) * Mat3d::Identity());
    }
    if (node->next[1])
    {
      double m1 = node->monopole;
      double m2 = node->next[1]->monopole;
      node->monopole = m1 + m2;
      node->vol_min  = node->vol_min.cwiseMin(node->next[1]->vol_min);
      node->vol_max  = node->vol_max.cwiseMax(node->next[1]->vol_max);
      Vec3d r1  = node->com;
      Vec3d r2  = node->next[1]->com;
      node->com      = (m1 * r1 + m2 * r2) / (m1 + m2);
      // r1 -= node->com;
      // r2 -= node->com;
      // node->quadpole = node->quadpole + node->next[1]->quadpole + m1 * (3.0 * r1 * r1.transpose() - r1.dot(r1) * Mat3d::Identity())
      //                 + m2 * (3.0 * r2 * r2.transpose() - r2.dot(r2) * Mat3d::Identity());
    }
  }

  return node;
}


//--------------------------------------------------------------------------------------------------
// Clear k-d tree recursively
//--------------------------------------------------------------------------------------------------
void KdTree::ClearRecursive(Node* node)
{
  if (node == nullptr)
    return;

  if (node->next[0])
    ClearRecursive(node->next[0]);

  if (node->next[1])
    ClearRecursive(node->next[1]);

#ifdef LEAF_SIZE
  if (node->is_leaf)
    delete[] node->leaf_idx;
#endif

  delete node;
}


//--------------------------------------------------------------------------------------------------
// SPH neighbor search recursively
//--------------------------------------------------------------------------------------------------
void KdTree::NeighborSearchRecursive(const Node* node, int query_idx, std::vector<Neighbor>* neighbors, int* neighbor_num, double kappa)
{
  if (node == nullptr)
  {
    return;
  }

#ifdef LEAF_SIZE
  if (node->is_leaf)
  {
    // for leaf node
    for (int i = 0; i < node->leaf_num; i++)
    {
      int train_idx = node->leaf_idx[i];
      const Particle& query = particles_[query_idx];
      const Particle& train = particles_[train_idx];

      if (query_idx != train_idx)
      {
        Vec3d xij = query.pos - train.pos;
        double dist2 = xij.dot(xij);
        double hij = 0.5 * (query.hsml + train.hsml);
        if (dist2 < pow2(hij * kappa))
        {
          double dist = sqrt(dist2);
          double kernel = kernel_.KernelValue(dist, hij);
          Vec3d gradKernel = kernel_.KernelGradient(dist, hij, xij);
          Vec3d vij = query.vel - train.vel;

          *neighbor_num += 1;
          if (*neighbor_num <= neighbors->size())
          {
            (*neighbors)[*neighbor_num - 1].ResetNeighbor(train_idx, dist, hij, xij, vij, kernel, gradKernel);
          }
          else
          {
            neighbors->push_back(Neighbor(train_idx, dist, hij, xij, vij, kernel, gradKernel));
          }
        }
      }
    }
    return;
  }
#endif

  const Particle& query = particles_[query_idx];
  const Particle& train = particles_[node->idx];

  if (query_idx != node->idx)
  {
    Vec3d xij = query.pos - train.pos;

    double dist2 = xij.dot(xij);
    double hij = 0.5 * (query.hsml + train.hsml);
    if (dist2 < pow2(hij * kappa))
    {
      double dist = sqrt(dist2);
      double kernel = kernel_.KernelValue(dist, hij);
      Vec3d gradKernel = kernel_.KernelGradient(dist, hij, xij);
      Vec3d vij = query.vel - train.vel;

      *neighbor_num += 1;
      if (*neighbor_num <= neighbors->size())
      {
        (*neighbors)[*neighbor_num - 1].ResetNeighbor(node->idx, dist, hij, xij, vij, kernel, gradKernel);
      }
      else
      {
        neighbors->push_back(Neighbor(node->idx, dist, hij, xij, vij, kernel, gradKernel));
      }
    }
  }

  int axis = node->axis;
  int dir = query.pos(axis) < train.pos(axis) ? 0 : 1;
  NeighborSearchRecursive(node->next[dir], query_idx, neighbors, neighbor_num, kappa);

  // Calculate the distance to hyperplane that is defined by the node point and the assigned dimension
  // Cut off if larger than (0.5 * (h + hmax) * kappa), which is the upper limit for a potential neighbor inside (kappa * hij)
  // Note that a tree-based neighbor search effiency depends on the cut-off, very much
  double diff = fabs(query.pos(axis) - train.pos(axis));
  if (diff < 0.5 * (hmax_ + query.hsml) * kappa)
  {
    NeighborSearchRecursive(node->next[!dir], query_idx, neighbors, neighbor_num, kappa);
  }
}


//--------------------------------------------------------------------------------------------------
// Calculate self-gravity recursively
//--------------------------------------------------------------------------------------------------
void KdTree::CalcSelfGravityRecursive(const Node* node, int query_idx, Vec3d* grav)
{
  if (node == nullptr)
  {
    return;
  }

#ifdef LEAF_SIZE
  if (node->is_leaf)
  {
    // for leaf 
    const Particle& query = particles_[query_idx];
    Vec3d xij = query.pos - node->com;
    double dist2 = xij.dot(xij);
    double box2  = pow2((node->vol_max - node->vol_min).maxCoeff());

    // open angle criterion
    if (box2 < dist2 * kOpenAngle2)
    {
      double dist = sqrt(dist2);
      Vec3d ir = xij / dist;
      // // monopole term
      // *grav -= node->monopole / dist2 * ir;
      // // quadpole term
      // *grav += (node->quadpole * ir - 2.5 * ir.dot(node->quadpole * ir) * ir) / pow2(dist2);
      *grav -= node->monopole * kernel_.KernelSoften(dist, query.hsml) * ir;
      return;
    }

    // cannot approximate
    for (int i = 0; i < node->leaf_num; i++)
    {
      int train_idx = node->leaf_idx[i];
      const Particle& train = particles_[train_idx];

      if (query_idx != train_idx && train.iflag >= 0)
      {
        Vec3d xij = query.pos - train.pos;
        double dist2 = xij.dot(xij);
        double dist = sqrt(dist2);
        // *grav -= train.mass / (dist * dist2) * xij;
        *grav -= train.mass * kernel_.KernelSoften(dist, query.hsml) * xij / dist;
      }
    }
    return;
  }
#endif

  const Particle& query = particles_[query_idx];
  const Particle& train = particles_[node->idx];

  // exclude the particle itself
  if (query_idx != node->idx)
  {
    // open angle criterion
    Vec3d xij = query.pos - node->com;
    double dist2 = xij.dot(xij);
    double box2  = pow2((node->vol_max - node->vol_min).maxCoeff()); // 0 for single particle node is allowed here
    
    if (box2 < dist2 * kOpenAngle2)
    {
      double dist = sqrt(dist2);
      Vec3d  ir   = xij / dist;
      // // monopole term
      // *grav -= node->monopole / dist2 * ir;
      // // quadpole term
      // *grav += (node->quadpole * ir - 2.5 * ir.dot(node->quadpole * ir) * ir) / pow2(dist2) ;
      *grav -= node->monopole * kernel_.KernelSoften(dist, query.hsml) * ir;
      return;
    }
    else
    {
      // add the node point's gravity if cannot approximate
      xij = query.pos - train.pos;
      dist2 = xij.dot(xij);
      double dist = sqrt(dist2);
      // *grav -= train.mass / (dist * dist2) * xij;
      *grav -= train.mass * kernel_.KernelSoften(dist, query.hsml) * xij / dist;
    }
  }

  // if cannot approximate
  if (node->next[0])
  {
    CalcSelfGravityRecursive(node->next[0], query_idx, grav);
  }
  if (node->next[1])
  {
    CalcSelfGravityRecursive(node->next[1], query_idx, grav);
  }
}