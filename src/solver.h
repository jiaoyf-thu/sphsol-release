#ifndef SPHSOL_SOLVER_H_
#define SPHSOL_SOLVER_H_

#include"kernel.h"
#include"particle.h"
#include"math.h"
#include"output.h"
#include"physics.h"
#include"kdtree.h"
#include"dbscan.h"
#include"grid.h"
#include"boundary.h"


enum class IntegratorEnum
{
  LEAP_FROG,
  PREDICT_CORRECT
};
enum class NeighborSearchEnum
{
  ALL_PAIR,
  GRID,
  KDTREE,
  HYBRID
};
enum class PredictCorrectEnum
{
  PREDICT,
  CORRECT
};
enum class ArtificialViscosityEnum
{
  NONE,
  STANDARD
};
enum class ArtificialHeatEnum
{
  NONE,
  STANDARD
};
enum class ArtificialStressEnum
{
  NONE,
  STANDARD
};
enum class XsphEnum
{
  NONE,
  STANDARD
};
enum class DensityEnum
{
  DENSITY_SUMMATION,
  CONTINUOUS_DENSITY
};
enum class CalcGravityEnum
{
  NONE,
  CONST,
  SPHERE,
  SELF_GRAVITY
};
enum class RunMode
{
  RELAX,
  SHOCK,
  LATERUN
};


//--------------------------------------------------------------------------------------------------
// SPH solver
//--------------------------------------------------------------------------------------------------
class Solver
{
public:
  Solver() : particles_(nullptr), neighbors_(nullptr), parts_(nullptr) {};
  ~Solver() { delete[] particles_; delete[] neighbors_; delete[] parts_; };

  void Run();
  void Initialize();

  Output out_;
  std::vector<int> input_flag_;

private:

  // Particles data
  Particle* particles_;
  std::vector<Neighbor>* neighbors_;
  Part* parts_;

  int particles_num_;  // Total particles number
  int parts_num_;

  // Current State
  int    current_step_;
  double current_time_;
  double current_hmax_;
  void   CheckStep(timeval* start, double total_time);

  double current_internal_energy_;
  double current_kinetic_energy_;

  double dt_;

  PredictCorrectEnum current_pc_step_; // only for predict-correct integrator

  // Settings
  std::vector<RunMode> run_modes_;
  RunMode current_run_mode_;
  
  IntegratorEnum integrator_;

  double total_time_;
  double dt_min_;
  double dt_max_;

  double dt_ref_;
  double dt_stp_;

  double late_time_;

  double hsml_ref_;
  double hsml_min_;
  double hsml_max_;

  Vec3d  volume_min_; // Particle is active only if inside the box
  Vec3d  volume_max_;

  Kernel kernel_;

  ArtificialViscosityEnum  artificial_viscosity_;

  ArtificialHeatEnum  artificial_heat_;

  ArtificialStressEnum  artificial_stress_;
  double avg_kernel_; // for artificial stress only

  NeighborSearchEnum neighbor_search_;

  XsphEnum xsph_;

  DensityEnum density_update_;

  Grid grid_;

  // // Gravity const
  // Vec3d gravity_;
  CalcGravityEnum calc_gravity_;
  // struct SphereGravity
  // {
  //   Vec3d center;
  //   double radius;
  //   double rho;
  //   // double escape_vel;
  // }sphere_gravity_;

  bool additional_force_; // whether to add additional forces

  // use theoretical pressure estimation before relaxation
  bool use_sphere_pressure_before_relax_;

  // Symmetric plane
  SymmetricPlane symmetric_plane_;

  // Rotation
  bool is_rotating_;
  struct Rotation
  {
    int part_id; // default is 0 as the target
    Vec3d pos0; // CoM position
    Vec3d omg0; // body's angular velocity
  }rotation_;

  double damp_time_;

  // Function in main loop
  void LoadParticles();
  void LoadSettings();
  void LoadMatDatabase();

  void UpdateTimeStep();
  void UpdateHsml();

  void BuildNeighbors(bool calc_grav = false);
  void UpdateNeighbors();

  void UpdatePressure();
  void UpdateArtificialViscosity();
  void UpdateArtificialHeat();
  void UpdateArtificialStress();
  void UpdateGravity();

  void Correction();

  void ComputeStressDerivative();
  void ComputeDamageDerivative();
  void ComputeDistentionDerivative();
  void ComputeDrhoDt();
  void ComputeDvDt();
  void ComputeDeDt();
  void ComputeDxDt();

  void AddForce();

  void Relax();
  void ShockRun();
  void LateRun();

  void WeibullDistribution();
  void UpdateCurrentEnergy();
};


#endif
