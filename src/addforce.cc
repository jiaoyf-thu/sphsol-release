#include"solver.h"


// //--------------------------------------------------------------------------------------------------
// // Costomed function to add extra forces: do nothing
// //--------------------------------------------------------------------------------------------------
// void Solver::AddForce() {}


// //--------------------------------------------------------------------------------------------------
// // Costomed function to add extra forces: constant gravity
// //--------------------------------------------------------------------------------------------------
// void Solver::AddForce()
// {
//   Vec3d grav_const = Vec3d(0.0, 0.0, -9.81);
  
//   #pragma omp parallel for schedule(guided)
//   for (int p = 0; p < particles_num_; p++)
//   {
//     if (particles_[p].iflag < 0) { continue; }

//     particles_[p].grav += grav_const;
//   }
// }


// //--------------------------------------------------------------------------------------------------
// // Costomed function to add extra forces: spherical gravity
// //--------------------------------------------------------------------------------------------------
// void Solver::AddForce()
// {
//   struct SphereGravity
//   {
//     Vec3d center;
//     double radius;
//     double rho;
//     double coef;
//   };

//   static const SphereGravity sphere_gravity = {
//     Vec3d(0.0, 0.0, 0.0),
//     6.9911E7,
//     1.326E3,
//     -4.0 / 3.0 * kPI * kG * 1.326E3   // coef = -4/3 π G ρ
//   };

//   #pragma omp parallel for schedule(guided)
//   for (int p = 0; p < particles_num_; p++)
//   {
//     if (particles_[p].iflag < 0) { continue; }

//     Vec3d grav_sphere = sphere_gravity.coef * (particles_[p].pos-sphere_gravity.center);
//     double dr = (particles_[p].pos-sphere_gravity.center).norm();
//     grav_sphere *= (dr > sphere_gravity.radius) ? pow3(sphere_gravity.radius/dr) : 1.0;
//     particles_[p].grav += grav_sphere;
//   }
// }


//--------------------------------------------------------------------------------------------------
// Costomed function to add extra forces: tidal forces
//--------------------------------------------------------------------------------------------------
// Table for interpolation
struct TimePos {
  std::vector<double> t;
  std::vector<Vec3d> pos;
};

// Interpolation using BinarySearch
Vec3d interpolatePos(double t, const TimePos& table)
{
  int idx = BinarySearch(table.t, t);

  // Clamp idx to valid interval
  if (idx >= static_cast<int>(table.t.size() - 1))
  {
    idx = table.t.size() - 2;
  }
  if (idx < 0)
  {
    idx = 0;
  }

  double u = (t - table.t[idx]) / (table.t[idx + 1] - table.t[idx]);
  return (1.0 - u) * table.pos[idx] + u * table.pos[idx + 1];
}

// Load table from CSV into struct-of-vectors
TimePos loadTableFromCSV(const std::string& filename)
{
  TimePos table;
  std::ifstream fin(filename);

  std::string line;
  bool header = true;
  while (std::getline(fin, line))
  {
    if (header) { header = false; continue; } // skip header

    std::stringstream ss(line);
    std::string cell;
    double t, x, y, z;

    std::getline(ss, cell, ','); t = std::stod(cell);
    std::getline(ss, cell, ','); x = std::stod(cell);
    std::getline(ss, cell, ','); y = std::stod(cell);
    std::getline(ss, cell, ','); z = std::stod(cell);

    table.t.push_back(t);
    table.pos.push_back(Vec3d(x, y, z));
  }

  return table;
}

// First-order tidal force calculation
void Solver::AddForce()
{
  // Load once and keep in function-local static
  static const TimePos table = loadTableFromCSV("../Input/ADD-Trajectory.csv");

  static const double kGmass = kG * 5.6834E26; // kg, Saturn
  static const double kRadii2 = pow2(6.0268E7); // m, Saturn
  // static const double kGmass = kG * 1.02409E26; // kg, Neptune
  // static const double kGmass = kG * 1.8982E27; // kg, Jupiter

  Vec3d com = interpolatePos(current_time_, table);
  // Vec3d com_old = interpolatePos(current_time_ - dt_, table);

  #pragma omp parallel for schedule(guided)
  for (int p = 0; p < particles_num_; p++)
  {
    // check if particle impact Saturn
    double dist2 = (particles_[p].pos + com).squaredNorm();
    particles_[p].iflag = (dist2 <= kRadii2) ? -1 : particles_[p].iflag;

    if (particles_[p].iflag < 0)
    {
      // freeze impacted particles in saturn centered frame
      particles_[p].pos = -com;
      continue;
    }

    double r_norm = com.norm();
    Vec3d  r_hat  = (r_norm > kEps) ? com / r_norm : Vec3d(0.0, 0.0, 0.0);
    Vec3d  tide   = (kGmass / pow3(r_norm)) * ((3.0*r_hat*r_hat.transpose() - Mat3d::Identity()) * particles_[p].pos);
    particles_[p].grav += tide;
  }

}



// //--------------------------------------------------------------------------------------------------
// // Costomed function to add extra forces: tidal forces exact expression
// //--------------------------------------------------------------------------------------------------
// // Table for interpolation
// struct TimePos {
//   std::vector<double> t;
//   std::vector<Vec3d> pos;
// };

// // Interpolation using BinarySearch
// Vec3d interpolatePos(double t, const TimePos& table)
// {
//   int idx = BinarySearch(table.t, t);

//   // Clamp idx to valid interval
//   if (idx >= static_cast<int>(table.t.size() - 1))
//   {
//     idx = table.t.size() - 2;
//   }
//   if (idx < 0)
//   {
//     idx = 0;
//   }

//   double u = (t - table.t[idx]) / (table.t[idx + 1] - table.t[idx]);
//   return (1.0 - u) * table.pos[idx] + u * table.pos[idx + 1];
// }

// // Load table from CSV into struct-of-vectors
// TimePos loadTableFromCSV(const std::string& filename)
// {
//   TimePos table;
//   std::ifstream fin(filename);

//   std::string line;
//   bool header = true;
//   while (std::getline(fin, line))
//   {
//     if (header) { header = false; continue; } // skip header

//     std::stringstream ss(line);
//     std::string cell;
//     double t, x, y, z;

//     std::getline(ss, cell, ','); t = std::stod(cell);
//     std::getline(ss, cell, ','); x = std::stod(cell);
//     std::getline(ss, cell, ','); y = std::stod(cell);
//     std::getline(ss, cell, ','); z = std::stod(cell);

//     table.t.push_back(t);
//     table.pos.push_back(Vec3d(x, y, z));
//   }

//   return table;
// }

// // Exact tidal force calculation
// void Solver::AddForce()
// {
//   // Load once and keep in function-local static
//   static const TimePos table = loadTableFromCSV("../Input/ADD-Trajectory.csv");
//   // static const double kGmass = kG * 1.8982E27; // kg, Jupiter
//   static const double kGmass = kG * 5.6834E26; // kg, Saturn
//   // static const double kGmass = kG * 1.02409E26; // kg, Neptune

//   Vec3d com = interpolatePos(current_time_, table);
//   double r_norm = com.norm();

//   #pragma omp parallel for schedule(guided)
//   for (int p = 0; p < particles_num_; p++)
//   {
//     if (particles_[p].iflag < 0) { continue; }
    
//     Vec3d  r_vec  = particles_[p].pos + com;
//     double r_vec_norm = r_vec.norm();
    
//     Vec3d  tide   = kGmass*(r_vec/pow3(r_vec_norm) - com/pow3(r_norm));
//     particles_[p].grav += tide;
//   }

// }

