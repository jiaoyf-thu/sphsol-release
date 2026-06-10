#ifndef SPHSOL_KERNEL_H_
#define SPHSOL_KERNEL_H_

#include"global.h"
#include"math.h"


enum class KernelEnum
{
  CUBIC_SPLINE
};


//--------------------------------------------------------------------------------------------------
// Cubic spline parameters
//--------------------------------------------------------------------------------------------------
struct CubicSplinePara
{
  CubicSplinePara() 
  {
    kappa = 2.0;
    k1    = 1.0 / (4.0 * kPI);
    k2    = 3.0 / (4.0 * kPI);
  };

  double kappa;
  double k1;
  double k2;
};


//--------------------------------------------------------------------------------------------------
// SPH kernel
//--------------------------------------------------------------------------------------------------
class Kernel
{
public:
  Kernel() { kernel_ = KernelEnum::CUBIC_SPLINE; };
  ~Kernel() {};

  double get_kappa();

  double KernelValue(double dist, double h);
  Vec3d  KernelGradient(double dist, double h, const Vec3d& xij);

  // gravity softening kernel
  double KernelSoften(double dist, double hij);

  // Kernel selection
  KernelEnum kernel_;

private:

  CubicSplinePara  cubic_spline_;

};


#endif
