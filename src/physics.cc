#include"physics.h"


//--------------------------------------------------------------------------------------------------
// Equation of state
//--------------------------------------------------------------------------------------------------
void EquatonOfState(Particle* particle, const Material& mat)
{
  double alpha = (mat.porosity == PorosityEnum::P_ALPHA) ? particle->alpha : 1.0;
  double alpha0 = (mat.porosity == PorosityEnum::P_ALPHA) ? mat.p_alpha.alpha0 : 1.0;

  double rho = particle->rho * alpha;
  double rho0 = particle->rho0 * alpha0;
  double energy = particle->energy;
  MY_ASSERT(std::isfinite(rho));
  MY_ASSERT(std::isfinite(rho0));
  MY_ASSERT(std::isfinite(energy));

  double pressure    = 0.0;
  double c_pressure2 = 0.0;
  double pmin        = (mat.yield == YieldEnum::LUND) ? -mat.lund.Yd0 * max(0.0, 1.0 - particle->plastic_strain) / mat.lund.mud : 0.0;

  double eta = rho / rho0;
  double mu = eta - 1.0;

  switch (mat.eos)
  {
  case(EosEnum::TILLOTSON):
  {
    const TillotsonPara& para = mat.tillotson;

    double nu = 1.0 / eta - 1.0;
    double omg = energy / (para.E0 * eta * eta) + 1.0;

    // There is a density limit for cold expanded state
    if (eta < 0.8 && energy < para.Eiv)
    {
      pressure = 0.0;
      // particle->damage = 1.0; // Directly set particle failure if over (cold or hot) expanded
    }
    else
    {
      // Common Tillotson equation of state
      if (eta >= 1.0 || energy <= para.Eiv)
      {
        pressure = (para.a + para.b / omg) * rho * energy + mu * (para.A + mu * para.B);
        // pressure = pressure < 0.0 ? pressure * (1.0 - particle->damage) : pressure;
        pressure = max(pressure, pmin);

        c_pressure2 = (1.0 + para.a + para.b / omg) * pressure / rho + (2.0 * energy - pressure / rho) * 
            para.b * (omg - 1.0) / pow2(omg) + (para.A + para.B * (pow2(eta) - 1.0)) / rho;

        MY_ASSERT(std::isfinite(pressure));
        MY_ASSERT(std::isfinite(c_pressure2));
      }
      else if (energy >= para.Ecv)
      {
        pressure = para.a * rho * energy + (para.b * rho * energy / omg + para.A * mu * exp(-para.beta * nu)) *
            exp(-para.alpha * nu * nu);
        pressure = max(pressure, pmin);
        
        c_pressure2 = (1.0 + para.a + para.b / omg * exp(-para.alpha * nu * nu)) * pressure / rho + 
            exp(-para.alpha * nu * nu) * (((2 * energy - pressure / rho) / (para.E0 * rho) + 
            2.0 * para.alpha * nu * omg / rho0) * para.b * rho * energy / (pow2(omg) * pow2(eta)) +
            (1.0 + mu / pow2(eta) * (para.beta + 2.0 * para.alpha * nu - eta)) * para.A / rho0 * exp(-para.beta * nu));

        // particle->damage = 1.0; // Directly set particle failure if over cold or hot expanded

        MY_ASSERT(std::isfinite(pressure));
        MY_ASSERT(std::isfinite(c_pressure2));
      }
      else if (energy > para.Eiv && energy < para.Ecv)
      {
        double pc = (para.a + para.b / omg) * rho * energy + mu * (para.A + mu * para.B);
        double pe = para.a * rho * energy + (para.b * rho * energy / omg + para.A * mu * exp(-para.beta * nu)) * 
            exp(-para.alpha * nu * nu);

        pressure = (pe * (energy - para.Eiv) + pc * (para.Ecv - energy)) / (para.Ecv - para.Eiv);
        pressure = max(pressure, pmin);

        double cc2 = (1.0 + para.a + para.b / omg) * pressure / rho + (2.0 * energy - pressure / rho) * para.b * 
            (omg - 1.0) / pow2(omg) + (para.A + para.B * (pow2(eta) - 1.0)) / rho;
        double ce2 = (1.0 + para.a + para.b / omg * exp(-para.alpha * nu * nu)) * pressure / rho + 
            exp(-para.alpha * nu * nu) * (((2 * energy - pressure / rho) / (para.E0 * rho) + 
            2.0 * para.alpha * nu * omg / rho0) * para.b * rho * energy / (pow2(omg) * pow2(eta)) +
            (1.0 + mu / pow2(eta) * (para.beta + 2.0 * para.alpha * nu - eta)) * para.A / rho0 * exp(-para.beta * nu));
        c_pressure2 = (ce2 * (energy - para.Eiv) + cc2 * (para.Ecv - energy)) / (para.Ecv - para.Eiv);

        // particle->damage = 1.0; // Directly set particle failure if over cold or hot expanded

        MY_ASSERT(std::isfinite(pressure));
        MY_ASSERT(std::isfinite(c_pressure2));
      }
    }
    break;
  }
  case(EosEnum::SIMPLIFIED_TILLOTSON):
  {
    pressure = mat.sim_tillotson.c * rho * energy + mat.sim_tillotson.A * mu;
    // pressure = mat.sim_tillotson.A * mu;
    pressure = max(pressure, pmin);

    // c_pressure2 = mat.sim_tillotson.c * (energy + pressure / rho) + mat.sim_tillotson.A / rho0;
    c_pressure2 = mat.sim_tillotson.A / rho;
    break;
  }
  case(EosEnum::GRUNEISEN):
  {
    break;
  }
  case(EosEnum::ELASTIC):
  {
    pressure = mat.bulk_modulus * mu;
    pressure = max(pressure, pmin);
    c_pressure2 = mat.bulk_modulus / rho;
    break;
  }
  default:
    break;
  }

  // Record the symbol of dPdt to dalpha_dt temporarily
  // particle->dalpha_dt = (pressure / alpha > particle->pressure) ? 1.0 : -1.0;

  // Update the pressure & sound speed
  // double min_pressure = mat.yield == YieldEnum::LUND ? -mat.lund.Yi0 : -mat.von_mises.Y0;
  // particle->pressure = max(pressure / alpha, min_pressure);
  particle->pressure = pressure / alpha;
  // Set the longitudinal sound speed no less than elastic sound speed
  c_pressure2 = max(c_pressure2, mat.bulk_modulus / rho0);
  particle->c_sound = sqrt(c_pressure2);

}



//--------------------------------------------------------------------------------------------------
// Equation of state
//--------------------------------------------------------------------------------------------------
double EquatonOfStateCold(const Particle* particle, const Material& mat)
{
  double alpha = (mat.porosity == PorosityEnum::P_ALPHA) ? particle->alpha : 1.0;
  double alpha0 = (mat.porosity == PorosityEnum::P_ALPHA) ? mat.p_alpha.alpha0 : 1.0;

  double rho = particle->rho * alpha;
  double rho0 = particle->rho0 * alpha0;
  double energy = particle->energy_cold;
  MY_ASSERT(std::isfinite(rho));
  MY_ASSERT(std::isfinite(rho0));
  MY_ASSERT(std::isfinite(energy));

  double pressure    = 0.0;
  double c_pressure2 = 0.0;

  double eta = rho / rho0;
  double mu = eta - 1.0;

  switch (mat.eos)
  {
  case(EosEnum::TILLOTSON):
  {
    const TillotsonPara& para = mat.tillotson;

    double nu = 1.0 / eta - 1.0;
    double omg = energy / (para.E0 * eta * eta) + 1.0;

    // There is a density limit for cold expanded state
    if (eta < 0.8 && energy < para.Eiv)
    {
      pressure = 0.0;
    }
    else
    {
      // Common Tillotson equation of state
      if (eta >= 1.0 || energy <= para.Eiv)
      {
        pressure = (para.a + para.b / omg) * rho * energy + mu * (para.A + mu * para.B);

        MY_ASSERT(std::isfinite(pressure));
      }
      else if (energy >= para.Ecv)
      {
        pressure = para.a * rho * energy + (para.b * rho * energy / omg + para.A * mu * exp(-para.beta * nu)) *
            exp(-para.alpha * nu * nu);
        
        MY_ASSERT(std::isfinite(pressure));
      }
      else if (energy > para.Eiv && energy < para.Ecv)
      {
        double pc = (para.a + para.b / omg) * rho * energy + mu * (para.A + mu * para.B);
        double pe = para.a * rho * energy + (para.b * rho * energy / omg + para.A * mu * exp(-para.beta * nu)) * 
            exp(-para.alpha * nu * nu);

        pressure = (pe * (energy - para.Eiv) + pc * (para.Ecv - energy)) / (para.Ecv - para.Eiv);

        MY_ASSERT(std::isfinite(pressure));
      }
    }
    break;
  }
  case(EosEnum::SIMPLIFIED_TILLOTSON):
  {
    pressure = mat.sim_tillotson.c * rho * energy + mat.sim_tillotson.A * mu;
    // pressure = mat.sim_tillotson.A * mu;
    
    break;
  }
  case(EosEnum::GRUNEISEN):
  {
    break;
  }
  case(EosEnum::ELASTIC):
  {
    pressure = mat.bulk_modulus * mu;
    break;
  }
  default:
    break;
  }

  pressure = pressure < 0.0 ? pressure * (1.0 - particle->damage) : pressure;

  return pressure / alpha;
}


//--------------------------------------------------------------------------------------------------
// Strength model
//--------------------------------------------------------------------------------------------------
double StrengthModel(const Material& mat, double damage, double J2, double pressure, double plastic_strain)
{
  double f = 1.0;
  switch (mat.yield)
  {
  case(YieldEnum::VON_MISES):
  {
    double yield_strength = mat.von_mises.Y0 * (1.0 - damage);
    // Yield conditon: Y / sqrt(3 * J2) > 1
    f = J2 > kEps ? min(1.0, yield_strength / sqrt(3.0 * J2)) : 1.0;
    break;
  }
  case(YieldEnum::LUND):
  {
    double weaken_factor = max(0.0, 1.0 - plastic_strain);
    double yield_strength_intact = mat.lund.Yi0 * weaken_factor + mat.lund.mui * pressure / (1.0 + mat.lund.mui * pressure / (mat.lund.Ym - mat.lund.Yi0));
    // Limit Yd <= Yi
    double yield_strength_damage = min(mat.lund.Yd0 * weaken_factor + mat.lund.mud * pressure, yield_strength_intact);
    // Limit Y >= 0.0
    double yield_strength = max(yield_strength_intact * (1.0 - damage) + yield_strength_damage * damage, 0.0);
    // Yield conditon: Y / sqrt(J2) > 1
    f = J2 > kEps ? min(1.0, yield_strength / sqrt(J2)) : 1.0;
    break;
  }
  default:
    break;
  }
  return f;
}


//--------------------------------------------------------------------------------------------------
// Shear damage
//--------------------------------------------------------------------------------------------------
double ShearDamage(double plastic_strain, double pressure, const Material& mat)
{
  if (mat.damage == DamageEnum::GRADY_KIPP)
  {
    double failure_strain = 0.01;
    if (pressure<mat.grady_kipp.p_bd)
    {
      failure_strain = max(0.01, 0.01+0.04*pressure/mat.grady_kipp.p_bd);
    }
    else if (pressure > mat.grady_kipp.p_bp)
    {
      failure_strain = 0.1 + 0.5*(pressure-mat.grady_kipp.p_bp)/mat.grady_kipp.p_bp;
    }
    else
    {
      failure_strain = 0.05 + 0.05*(pressure-mat.grady_kipp.p_bd)/(mat.grady_kipp.p_bp-mat.grady_kipp.p_bd);
    }
    return min(1.0, plastic_strain/failure_strain);
  }
  else
  {
    return 0.0;
  }
}
