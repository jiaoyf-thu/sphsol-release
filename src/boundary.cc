#include"boundary.h"


//--------------------------------------------------------------------------------------------------
// Initialize the symmetric plane (Ax + By +Cz + D = 0)
//--------------------------------------------------------------------------------------------------
void SymmetricPlane::Init(double* coef)
{
  a_ = coef[0];
  b_ = coef[1];
  c_ = coef[2];
  d_ = coef[3];

  sum_square_ = pow2(a_) + pow2(b_) + pow2(c_);
  sqrt_sum_square_ = sqrt(sum_square_);

  unit_n << a_, b_, c_;
  unit_n.normalize();

  is_initialised_ = true;
}


//--------------------------------------------------------------------------------------------------
// Calculate the particle's distance to the symmetric plane (Ax + By +Cz + D = 0)
//--------------------------------------------------------------------------------------------------
// Note that all real particles must be (Ax + By +Cz + D > 0) to make sure positive distance values
// While particles at (Ax + By +Cz + D < 0) return negetive distance values and errors
//--------------------------------------------------------------------------------------------------
double SymmetricPlane::GetDistance(const Particle& particle)
{
  double dist = (a_ * particle.pos(0) + b_ * particle.pos(1) + c_ * particle.pos(2) + d_) / sqrt_sum_square_;
  // return max(dist, -dist);
  return dist;
}


//--------------------------------------------------------------------------------------------------
// Generate a mirroring ghost particle
//--------------------------------------------------------------------------------------------------
// Note that the ghost particle is not stored here
// So that we don't need to change other ghost particle states, except its hsml, pos & vel
//--------------------------------------------------------------------------------------------------
void SymmetricPlane::GetGhostParticle(const Particle& real_particle, Particle& ghost_particle)
{
  // // Copy all values
  // ghost_particle = real_particle;
  ghost_particle.hsml = real_particle.hsml;

  // Mirror ghost particle's pos & vel
  double coef = 2.0 * (a_ * real_particle.pos(0) + b_ * real_particle.pos(1) + c_ * real_particle.pos(2) + d_) / sum_square_;
  ghost_particle.pos = real_particle.pos - Vec3d(a_, b_, c_) * coef;
  ghost_particle.vel = real_particle.vel - 2.0 * (unit_n.dot(real_particle.vel)) * unit_n;

  // // Change ghost particle's stress
  // ghost_particle.dev_stress = 2.0 * (Mat3d::Identity().array() * real_particle.dev_stress.array()).matrix() - real_particle.dev_stress;
  // ghost_particle.stress = 2.0 * (Mat3d::Identity().array() * real_particle.stress.array()).matrix() - real_particle.stress;

}


//--------------------------------------------------------------------------------------------------
// Generate a mirroring real particle
//--------------------------------------------------------------------------------------------------
void SymmetricPlane::GetMirrorParticle(Particle& real_particle)
{
  // Mirror ghost particle's pos & vel
  double coef = 2.0 * (a_ * real_particle.pos(0) + b_ * real_particle.pos(1) + c_ * real_particle.pos(2) + d_) / sum_square_;
  real_particle.pos = real_particle.pos - Vec3d(a_, b_, c_) * coef;
  real_particle.vel = real_particle.vel - 2.0 * (unit_n.dot(real_particle.vel)) * unit_n;

  // Change ghost particle's stress
  real_particle.dev_stress = 2.0 * (Mat3d::Identity().array() * real_particle.dev_stress.array()).matrix() - real_particle.dev_stress;
  real_particle.stress = 2.0 * (Mat3d::Identity().array() * real_particle.stress.array()).matrix() - real_particle.stress;
}
