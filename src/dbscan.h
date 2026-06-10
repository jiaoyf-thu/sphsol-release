#ifndef SPHSOL_DBSCAN_H_
#define SPHSOL_DBSCAN_H_

#include"global.h"
#include"particle.h"


//--------------------------------------------------------------------------------------------------
// DBSCAN algorithm
//--------------------------------------------------------------------------------------------------
class DBSCAN
{
public:

  DBSCAN(const Particle* particles, const std::vector<Neighbor>* neighbors, int particles_num,
      std::vector<int>& cluster_list, double hsml, int min_cluster_size = 4) :
      particles_(particles), neighbors_(neighbors), particles_num_(particles_num),
      cluster_list_(cluster_list), hsml_ref_(hsml), min_cluster_size_(min_cluster_size)
  {
    // Initialize the cluster list by 0
    cluster_list_.resize(particles_num_, 0);
    cluster_history_.resize(particles_num_, 0);
  };
  ~DBSCAN() {};

  void Run();

private:

  double hsml_ref_;

  int min_cluster_size_; // Default: 4
  std::vector<int>& cluster_list_; // Store all particles' cluster id (0: INIT, -1: NOISE, positive int: clusterID)
  std::vector<int> cluster_history_; // Store if a particle has been a NOISE (0: not so far, -1: yes it is)

  const Particle* particles_;
  const std::vector<Neighbor>* neighbors_;
  int particles_num_;

  // Functions
  int ExpandCluster(int p, int clusterID);
  std::vector<int> CalculateCluster(int p, bool undamaged_only=false);
};


#endif
