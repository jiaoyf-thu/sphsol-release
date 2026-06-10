#ifndef SPHSOL_PARTICLE_H_
#define SPHSOL_PARTICLE_H_

#include"global.h"
#include"material.h"


//--------------------------------------------------------------------------------------------------
// SPH particle data
//--------------------------------------------------------------------------------------------------
struct Particle
{
  Vec3d pos;                      // position for physics equations
  Vec3d dxdt;

  Vec3d vel;                      // velocity
  Vec3d dvdt;

  Mat3d dev_stress;               // deviatoric stress tensor
  Mat3d ddev_stress_dt;
  Mat3d stress;                   // stress tensor
  Mat3d tension;                  // tension tensor, for calculating artificial stress

  Mat3d strain_rate;
  Mat3d rotate_rate;

  Mat3d correction;               // correction tensor

  int iflag;                      // particle state (1: normal, 0: outside, -1: freeze)

  int part_id;

  double mass;                    // mass

  double rho0;                    // initial density
  double rho;                     // density
  double drhodt;

  double pressure;                // pressure
  double pressure_peak;           // peak pressure

  Vec3d grav;                     //  gravity term: g_i = sum(G*m_j/r/r) * unit_vec(r_j-r_i)
  
  double energy;                  // internal energy
  double energy_cold;             // cold internal energy
  double energy_thermal_peak;     // peak <thermal internal energy = energy - energy_cold>
  double dedt;

  double hsml;                    // smooth length
  double c_sound;                 // sound speed

  double damage;                  // damage
  double damage_root3;            // damage^(1/3)
  double ddamage_root3_dt;
  double plastic_strain;

  double alpha;                   // distention for porous state
  double alpha_inv;               // inverse of alpha
  double dalpha_dt;
  double dpdt;
  double dalpha_drho;

  double Hi;                      // artifitial heat

  struct Flaw
  {
    int ntot;
    double strain_min, strain_max;
    double damage_root3_max;
    double m0;
  }flaws;                         // weibull flaws

  struct Frag
  {
    double area, dAdt;
    double Lm;                    // peak fragment size
  }frag;                          // fragment size estimate

};


//--------------------------------------------------------------------------------------------------
// Neighbor
//--------------------------------------------------------------------------------------------------
struct Neighbor
{
  Neighbor() {};
  Neighbor(int id, double dist, double hij, Vec3d xij, Vec3d vij, double kernel = 0.0, Vec3d gradKernel = Vec3d::Zero(), double IIij = 0.0, Mat3d zeta = Mat3d::Zero()) :
      id(id), dist(dist), hij(hij), kernel_value(kernel), xij(xij), vij(vij), kernel_gradient(gradKernel), IIij(IIij), zeta(zeta), iflag(1) {};
  
  // Reset neighbor data
  void ResetNeighbor(int id_, double dist_, double hij_, Vec3d xij_, Vec3d vij_, double kernel = 0.0, Vec3d gradKernel = Vec3d::Zero(), double IIij_ = 0.0, Mat3d zeta_ = Mat3d::Zero())
  {
    id              = id_;
    dist            = dist_;
    hij             = hij_;
    IIij            = IIij_;
    zeta            = zeta_;
    kernel_value    = kernel;
    xij             = xij_;
    vij             = vij_;
    kernel_gradient = gradKernel;
  };

  int     id;               // neighbor particle id
  int     iflag;            // 1: normal neighbor, -1: ghost neighbor
  double  dist;             // distance
  double  hij;              // 0.5 * (hi + hj)
  double  IIij;             // aitificial viscosity
  Mat3d   zeta;             // artificial stress
  double  kernel_value;     // kernel value

  Vec3d   xij;              // (xi - xj)
  Vec3d   vij;              // (vi - vj)
  Vec3d   kernel_gradient;  // kernel gardient
};


//--------------------------------------------------------------------------------------------------
// Part
//--------------------------------------------------------------------------------------------------
struct Part
{
  int mat_id; // material id

  std::vector<int> particle_id; // all linked particles' id
};


#endif
