#include "dbscan.h"


//--------------------------------------------------------------------------------------------------
// DBSCAN main loop
//--------------------------------------------------------------------------------------------------
void DBSCAN::Run()
{
  // All damaged particles are initialized as NOISE
  for(int p = 0; p < particles_num_; p++)
  {
    cluster_list_[p] = (particles_[p].damage > 1.0 - kEps) ? -1 : 0;
    cluster_history_[p] = cluster_list_[p];
  }
  
  // Cluster algorithm
  int cluster_id = 1;
  for(int p = 0; p < particles_num_; p++)
  {
    if (cluster_list_[p] == 0)
    {
      if (ExpandCluster(p, cluster_id) != 0)
      {
        cluster_id += 1;
      }
    }
  }

  // Record cluster size & id
  struct ClusterOrder
  {
    ClusterOrder(int size_ = 0, int id_ = 0) : size(size_), id(id_) {};
    int size;
    int id;
  };

  // Sort clusters by size descending order except NOISE
  std::vector<ClusterOrder> cluster_order(cluster_id);
  for (int id = 0; id < cluster_id; id++)
  {
    cluster_order[id].id = id;
  }
  cluster_order[0].id = -1;
  cluster_order[0].size = -1;
  for (auto cluster : cluster_list_)
  {
    if (cluster == -1)
    {
      // cluster_order[0].size++;
      continue;
    }
    cluster_order[cluster].size++;
  }
  std::sort(cluster_order.begin() + 1, cluster_order.end(), [](ClusterOrder& A, ClusterOrder& B){ return (A.size > B.size); });
  // Now reset the cluster_size = new index
  for (int id = 1; id < cluster_id; id++)
  {
    cluster_order[id].size = id;
  }
  // Change the cluster_id in the list
  for (int p = 0; p < cluster_list_.size(); p++)
  {
    int value{cluster_list_[p]};
    cluster_list_[p] = std::find_if(cluster_order.begin(), cluster_order.end(), [value](ClusterOrder& A){ return (A.id == value); })->size;
  }
}


//--------------------------------------------------------------------------------------------------
// DBSCAN expand cluster from CORE point
//--------------------------------------------------------------------------------------------------
int DBSCAN::ExpandCluster(int p, int cluster_id)
{    
  std::vector<int> cluster_seeds = CalculateCluster(p,true);

  // Too small CORE is marked as NOISE
  if (cluster_seeds.size() < min_cluster_size_)
  {
    cluster_list_[p] = -1;
    cluster_history_[p] = -1;
    return 0;
  }
  else
  {
    // Assign clusterID to all neighbors
    for(auto iter_seeds = cluster_seeds.begin(); iter_seeds != cluster_seeds.end(); ++iter_seeds)
    {
      cluster_list_.at(*iter_seeds) = cluster_id;
    }
    // Delete itself
    cluster_seeds.erase(cluster_seeds.begin());
    
    // Here clusterseed's size is increasing by adding neighbors' neighbors' ...
    for(int i = 0, n = cluster_seeds.size(); i < n; ++i)
    {
      if (cluster_history_[cluster_seeds[i]] == -1)
      {
        continue;
      }
      std::vector<int> cluster_neighors = CalculateCluster(cluster_seeds[i]);
      if (cluster_neighors.size() < min_cluster_size_)
      {
        cluster_history_[cluster_seeds[i]] = -1;
        continue;
      }

      for (auto iter_neighors = cluster_neighors.begin(); iter_neighors != cluster_neighors.end(); ++iter_neighors)
      {
        // Noise point cannot expand the cluster
        if (cluster_history_.at(*iter_neighors) == -1)
        { 
          cluster_list_.at(*iter_neighors) = cluster_id;
          continue;
        }
        // Add current neighbors into clusterSeeds resursively
        if (cluster_list_.at(*iter_neighors) == 0)
        {
          cluster_seeds.push_back(*iter_neighors);
          n = cluster_seeds.size();
          cluster_list_.at(*iter_neighors) = cluster_id;
        }
      }

    }
    return 1;
  }
}


//--------------------------------------------------------------------------------------------------
// Calculate neighbors of specific point including itself
//--------------------------------------------------------------------------------------------------
std::vector<int> DBSCAN::CalculateCluster(int p, bool undamaged_only)
{
  std::vector<int> cluster_index;
  cluster_index.reserve(10);
  cluster_index.push_back(p); // The first neighbor is the point itself

  for(auto neighbor : neighbors_[p])
  {
    // if (neighbor.dist < 1.01 * neighbor.hij)
    if (neighbor.dist < 1.01 * hsml_ref_)
    {
      if (undamaged_only && particles_[neighbor.id].damage < 1.0)
      {
        cluster_index.push_back(neighbor.id);
      }
      else if (!undamaged_only)
      {
        cluster_index.push_back(neighbor.id);
      }
    }
  }
  return cluster_index;
}

