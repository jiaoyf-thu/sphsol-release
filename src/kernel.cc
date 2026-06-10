#include"kernel.h"


//--------------------------------------------------------------------------------------------------
// SPH kernel range (kappa * hsml)
//--------------------------------------------------------------------------------------------------
double Kernel::get_kappa()
{
  switch (kernel_)
  {
  case(KernelEnum::CUBIC_SPLINE):
  {
    return(cubic_spline_.kappa);
    break;
  }
  default:
  {
    return 0.0;
    break;
  }
  }
}


//--------------------------------------------------------------------------------------------------
// Caluculate kernel value
//--------------------------------------------------------------------------------------------------
double Kernel::KernelValue(double dist, double hij)
{
  double wij = 0.0;
  double R = dist / hij;

  switch (kernel_)
  {
  case(KernelEnum::CUBIC_SPLINE):
  {
    if (R >= 0.0 && R < 1.0)
    {
      wij = 4.0 - 6.0 * pow2(R) + 3.0 * pow3(R);
    }
    else if (R >= 1.0 && R <= 2.0)
    {
      wij = pow3(2.0 - R);
    }
    wij *= cubic_spline_.k1 / pow3(hij);

    break;
  }
  default:
    break;
  }

  return wij;
}


//--------------------------------------------------------------------------------------------------
// Calculate kernel gradient
//--------------------------------------------------------------------------------------------------
Vec3d Kernel::KernelGradient(double dist, double hij, const Vec3d& xij)
{
  Vec3d dwdx = Vec3d::Zero();
  double R = dist / hij;
  double dwdr = 0.0;

  switch (kernel_)
  {
  case(KernelEnum::CUBIC_SPLINE):
  {
    if (R >= 0.0 && R < 1.0)
    {
      dwdr = -4.0 * R + 3.0 * pow2(R);
    }
    else if (R >= 1.0 && R < 2.0)
    {
      dwdr = -pow2(2.0 - R);
    }
    dwdr *= cubic_spline_.k2;
    dwdx = xij * dwdr / (dist * pow4(hij));
    break;
  }
  default:
    break;
  }

  return dwdx;
}


//--------------------------------------------------------------------------------------------------
// Caluculate gravity softening kernel
//--------------------------------------------------------------------------------------------------
double Kernel::KernelSoften(double dist, double hij)
{
  double phi_p = 0.0;
  double R = dist / hij;
  double inv_h2 = 1.0 / (hij * hij);

  switch (kernel_)
  {
  case(KernelEnum::CUBIC_SPLINE):
  {
    if (R < 1.0)
      {
        phi_p = (4.0/3.0 * R - 1.2 * R*R*R + 0.5 * R*R*R*R) * inv_h2;
      }
      else if (R < 2.0)
      {
        phi_p = (8.0/3.0 * R - 3.0 * R*R + 1.2 * R*R*R - 1.0/6.0 * R*R*R*R - 1.0/(15.0 * R*R)) * inv_h2;
      }
      else
      {
        phi_p = 1.0 / (dist * dist);
      }
    break;
  }
  default:
    break;
  }

  return phi_p;
}
