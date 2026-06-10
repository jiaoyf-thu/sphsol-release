#ifndef SPHSOL_GRID_H_
#define SPHSOL_GRID_H_

#include"global.h"
#include"math.h"
#include"kernel.h"
#include"particle.h"


//--------------------------------------------------------------------------------------------------
// Grid for spatial partition
//--------------------------------------------------------------------------------------------------
class Grid
{
public:

  Grid(): particles_(nullptr), cells_(nullptr), is_initialised_(false), cell_num_(0) {};
  ~Grid() { Clear(); };

  void InitGrid(double cell_size, Vec3d vol_min, Vec3d vol_max, Particle* particles, int particles_num, double hmax, Kernel& kernel);
  void Clear(bool free_memeory = true);

  void Insert(int pid);

  void NeighborSearch(int pid, std::vector<Neighbor>* neighbors, int* neighbor_num);
  void getCellNeighbor(int pid, std::vector<int>* cell_neighbor);

  void setHmax(double hmax) { hmax_ = hmax; };
  bool isInitialized() { return is_initialised_; };

private:


  int getCellIndex(Vec3i grid_id);
  int getCellIndex(int ix, int iy, int iz);
  Vec3i getGridIndex(const Vec3d& pos);

  std::vector<int>*  cells_; // Each cell contains its inside particles' index

  Particle* particles_;
  int particles_num_;

  Kernel kernel_;  // SPH kernel to calculate neighbors information

  Vec3d vol_min_; // Where the grid starts
  Vec3i grid_num_; // Grid numbers in XYZ directions
  int cell_num_; // Cell numbers = grid_num_(0) * grid_num_(1) * grid_num_(2)
  double cell_size_; // Length of each cell
  double hmax_;

  bool is_initialised_;
};  


#endif