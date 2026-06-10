#ifndef SPHSOL_BOUNDARY_H_
#define SPHSOL_BOUNDARY_H_

#include"global.h"
#include"math.h"
#include"particle.h"


//--------------------------------------------------------------------------------------------------
// Symmetric plane (Ax + By +Cz + D = 0, where A^2 + B^2 + C^2 > 0)
//--------------------------------------------------------------------------------------------------
// Note that all real particles must be (Ax + By +Cz + D > 0) to make sure positive distance values
//--------------------------------------------------------------------------------------------------
class SymmetricPlane
{
public:
  SymmetricPlane(): is_initialised_(false) {};
  ~SymmetricPlane(){};

  void Init(double* coef);
  bool isInitialized() { return is_initialised_; };
  
  double GetDistance(const Particle& real_particle);
  void GetGhostParticle(const Particle& real_particle, Particle& ghost_particle);
  void GetMirrorParticle(Particle& real_particle);

private:
  double a_, b_, c_, d_;
  Vec3d unit_n;
  double sum_square_; // A^2 + B^2 + C^2
  double sqrt_sum_square_; // sqrt(A^2 + B^2 + C^2)
  bool is_initialised_;
};


#endif