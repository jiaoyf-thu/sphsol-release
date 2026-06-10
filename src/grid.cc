#include"grid.h"


//------------------------------------------------------------------------------
// Initialize grid
//------------------------------------------------------------------------------
void Grid::InitGrid(double cell_size, Vec3d vol_min, Vec3d vol_max, Particle* particles, int particles_num, double hmax, Kernel& kernel)
{
  // Clear previous grid
  Clear();

  // Initialize the grid
  particles_ = particles;
  particles_num_ = particles_num;
  hmax_ = hmax;
  kernel_ = kernel;
  vol_min_ = vol_min;
  cell_size_ = cell_size;

  // Compute number of grid on each axis
  int resX = static_cast<int>(ceil((vol_max(0) - vol_min(0)) / cell_size)) + 1;
  int resY = static_cast<int>(ceil((vol_max(1) - vol_min(1)) / cell_size)) + 1;
  int resZ = static_cast<int>(ceil((vol_max(2) - vol_min(2)) / cell_size)) + 1;
  grid_num_ << resX, resY, resZ;
  cell_num_ = resX * resY * resZ;

  // Allocate memory
  cells_ = new std::vector<int>[cell_num_];
  for(int c = 0; c < cell_num_; c++)
  {
    cells_[c].reserve(10);
  }

  is_initialised_ = true;
}


//------------------------------------------------------------------------------
// Clear all grid infomation
//------------------------------------------------------------------------------
void Grid::Clear(bool free_memeory)
{
  if (free_memeory)
  {
    for(int c = 0; c < cell_num_; c++)
    {
      cells_[c].clear();
      std::vector<int>().swap(cells_[c]);
    }
  }
  else
  {
    for(int c = 0; c < cell_num_; c++)
    {
      cells_[c].clear();
    }
  }
}


//------------------------------------------------------------------------------
// Insert particle into grid
//------------------------------------------------------------------------------
void Grid::Insert(int pid)
{
  // Get Index
  Vec3i grid_idx = getGridIndex(particles_[pid].pos);

  int cell_idx = getCellIndex(grid_idx);
  cells_[cell_idx].push_back(pid);

  // // Make sure we are not out of the grid
  // if (grid_idx(0) >= 0 && grid_idx(0) < grid_num_(0) && grid_idx(1) >= 0 && grid_idx(1) < grid_num_(1) 
  //     && grid_idx(2) >= 0 && grid_idx(2) < grid_num_(2) && cell_idx >= 0 && cell_idx < cell_num_)
  // {
  //   cells_[cell_idx].push_back(pid);
  // }
}


//------------------------------------------------------------------------------
// Get all cells around a particle inside kappa*(h+hmax)/2
//------------------------------------------------------------------------------
void Grid::getCellNeighbor(int pid, std::vector<int>* cell_neighbor)
{
  // Get min/max indices
  double radii = 0.5 * (particles_[pid].hsml + hmax_) * kernel_.get_kappa();
  Vec3d pos_min = particles_[pid].pos - radii * Vec3d::Ones();
  Vec3d pos_max = particles_[pid].pos + radii * Vec3d::Ones();
  Vec3i grid_min = getGridIndex(pos_min);
  Vec3i grid_max = getGridIndex(pos_max);

  // // Clamp indices
  // for (int i = 0; i < 3; i++)
  // {
  //   grid_min(i) = max(grid_min(i), 0);
  //   grid_max(i) = min(grid_max(i), grid_num_(i) - 1);
  // }

  // Traverse cells
  cell_neighbor->clear();
  for (int ix = grid_min(0); ix <= grid_max(0); ix++)
  {
    for (int iy = grid_min(1); iy <= grid_max(1); iy++)
    {
      for (int iz = grid_min(2); iz <= grid_max(2); iz++)
      {
        cell_neighbor->push_back(getCellIndex(ix, iy, iz));
      }
    }
  }

}


//------------------------------------------------------------------------------
// Search neighor particles using the grid method
//------------------------------------------------------------------------------
void Grid::NeighborSearch(int pid, std::vector<Neighbor>* neighbors, int* neighbor_num)
{
  double kappa = kernel_.get_kappa();

  // Get cell neighbor
  std::vector<int> cell_neighbor;
  cell_neighbor.reserve(27);
  getCellNeighbor(pid, &cell_neighbor);

  // Search neighbors from cells
  const Particle& query = particles_[pid];
  *neighbor_num = 0;
  for (int c : cell_neighbor)
  {
    for (int nid : cells_[c])
    {
      const Particle& train = particles_[nid];

      if (pid != nid)
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
            (*neighbors)[*neighbor_num - 1].ResetNeighbor(nid, dist, hij, xij, vij, kernel, gradKernel);
          }
          else
          {
            neighbors->push_back(Neighbor(nid, dist, hij, xij, vij, kernel, gradKernel));
          }
        }
      }
    }
  }

}


//------------------------------------------------------------------------------
// Calculate cell index from grid index
//------------------------------------------------------------------------------
int Grid::getCellIndex(int ix, int iy, int iz)
{
  return ix + (iy * grid_num_(0)) + (iz * grid_num_(1) * grid_num_(0));
}
int Grid::getCellIndex(Vec3i grid_id)
{
  return grid_id(0) + (grid_id(1) * grid_num_(0)) + (grid_id(2) * grid_num_(1) * grid_num_(0));
}


//------------------------------------------------------------------------------
// Calculate grid index from position
//------------------------------------------------------------------------------
Vec3i Grid::getGridIndex(const Vec3d& pos)
{
  Vec3i grid_id;

  // Clamp indices
  for (int i = 0; i < 3; i++)
  {
    grid_id(i) = static_cast<int>(floor((pos(i) - vol_min_(i)) / cell_size_));
    grid_id(i) = middleOfThree(grid_id(i), 0, grid_num_(i) - 1);
  }

  return grid_id;
}
