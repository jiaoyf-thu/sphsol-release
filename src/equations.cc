#include"solver.h"


//--------------------------------------------------------------------------------------------------
// Initialize the weibull flaws distribution for grady-kipp damage model
//--------------------------------------------------------------------------------------------------
void Solver::WeibullDistribution()
{
  for (int partID = 0; partID < parts_num_; partID++)
  {
    const Part& part = parts_[partID];

    if (mat_database[part.mat_id].damage == DamageEnum::GRADY_KIPP)
    {
      int particles_num = part.particle_id.size();

      // Total volume
      double volume = 0.0;
      for (int p = 0; p < particles_num; p++)
      {
        volume += particles_[part.particle_id[p]].mass / particles_[part.particle_id[p]].rho0;
      }
      // Use (Vs = V / alpha0)
      if (mat_database[part.mat_id].porosity == PorosityEnum::P_ALPHA)
      {
        volume = volume / mat_database[part.mat_id].p_alpha.alpha0;
      }

      double strain_activate = 0.0; // activated strain
      double m_inverse = 1.0 / mat_database[part.mat_id].grady_kipp.m;
      double weibull_coef = 1.0 / (pow(mat_database[part.mat_id].grady_kipp.k, m_inverse) * pow(volume, m_inverse));

      srand(MY_SEED);
      for (int p = 0; p < particles_num; p++)
      {
        int particle_id = part.particle_id[p];
        Particle& particle = particles_[particle_id];

        // double rand_flaw_num = static_cast<double>(rand()) / RAND_MAX;
        // particle.flaws.ntot = BinarySerarch(poisson_list, rand_flaw_num) + 1;
        particle.flaws.ntot = max(PoissonDistribution(log((double)particles_num)), 1); // Poisson distribution with lambda = log(N)

        double flaw_first_activate, flaw_last_activate;
        double rand_num = static_cast<double>(rand()) / RAND_MAX;

        flaw_first_activate = -static_cast<double>(particles_num) * log(1.0 - rand_num);
        flaw_last_activate  = static_cast<double>(particles_num) * log(1.0 + static_cast<double>(particles_num) * rand_num);

        particle.flaws.strain_min = weibull_coef * pow(min(flaw_first_activate, flaw_last_activate), m_inverse);
        particle.flaws.strain_max = weibull_coef * pow(max(flaw_first_activate, flaw_last_activate), m_inverse);

        // Ensure that m0 >= 1
        particle.flaws.strain_max = min(particle.flaws.strain_max, particle.flaws.strain_min * particle.flaws.ntot);
      }

      // Here: m0 = log(ntot) / log(strainMAX / strainMIN)
      for (int p = 0; p < part.particle_id.size(); p++)
      {
        Particle& particle = particles_[part.particle_id[p]];

        if (particle.flaws.ntot == 1)
        {
          particle.flaws.m0 = 1.0;
        }
        else
        {
          particle.flaws.m0 = log(particle.flaws.ntot) / log(particle.flaws.strain_max / particle.flaws.strain_min);
        }

        MY_ASSERT(particle.flaws.ntot > 0);
        MY_ASSERT(particle.flaws.m0 >= 1.0);
        MY_ASSERT(particle.flaws.strain_min > 0.0 && particle.flaws.strain_max >= particle.flaws.strain_min);
      }
    }

  }
}


//--------------------------------------------------------------------------------------------------
// Update pressure with porosity and eos
//--------------------------------------------------------------------------------------------------
void Solver::UpdatePressure()
{
  switch (current_run_mode_)
  {

    // use EOS
    case (RunMode::SHOCK):
    {
      #pragma omp parallel for schedule(guided)
      for (int p = 0; p < particles_num_; p++)
      {
        if (particles_[p].iflag < 0)
        {
          // particles_[p].pressure = 0.0;
          continue;
        }
        Material& mat = mat_database[parts_[particles_[p].part_id].mat_id];

        EquatonOfState(&particles_[p], mat);
      }
      break;
    }

    // suppose elastic
    case (RunMode::RELAX):
    case (RunMode::LATERUN):
    {
      #pragma omp parallel for schedule(guided)
      for (int p = 0; p < particles_num_; p++)
      {
        if (particles_[p].iflag < 0)
        {
          continue;
        }
        Material& mat = mat_database[parts_[particles_[p].part_id].mat_id];
        
        double alpha = (mat.porosity == PorosityEnum::P_ALPHA) ? particles_[p].alpha : 1.0;
        double alpha0 = (mat.porosity == PorosityEnum::P_ALPHA) ? mat.p_alpha.alpha0 : 1.0;
        double mu = (particles_[p].rho * alpha) / (particles_[p].rho0 * alpha0) - 1.0;
        double pressure = mat.bulk_modulus * particles_[p].alpha_inv * mu;

        particles_[p].pressure = pressure < 0.0 ? pressure * (1.0 - particles_[p].damage) : pressure;
        particles_[p].c_sound = sqrt(mat.bulk_modulus * particles_[p].alpha_inv / particles_[p].rho);
      }
      break;
    }
    
    default:
      break;
  }

}


//--------------------------------------------------------------------------------------------------
// Update artificial viscosity
//--------------------------------------------------------------------------------------------------
void Solver::UpdateArtificialViscosity()
{
  switch (artificial_viscosity_)
  {
  case(ArtificialViscosityEnum::STANDARD):
  {
    double alpha = 1.0;
    double beta  = 2.0;

    #pragma omp parallel for schedule(guided)
    for (int p = 0; p < particles_num_; p++)
    {
      if (particles_[p].iflag < 0)
      {
        continue;
      }
      for (int n = 0; n < neighbors_[p].size(); n++)
      {
        Neighbor* neighbor = &neighbors_[p][n];

        double temp = neighbor->vij.dot(neighbor->xij);
        if (temp < 0.0)
        {
          double phi     = 0.1 * neighbor->hij;
          double uij     = neighbor->hij * temp / (pow2(neighbor->dist) + pow2(phi));
          double cij     = 0.5 * (particles_[p].c_sound + particles_[neighbor->id].c_sound);
          double rhoij   = 0.5 * (particles_[p].rho + particles_[neighbor->id].rho);

          neighbor->IIij = (-alpha * cij * uij + beta * pow2(uij)) / rhoij;

          MY_ASSERT(std::isfinite(neighbor->IIij));
        }
        else
        {
          neighbor->IIij = 0.0;
        }
      }
    }
    break;
  }
  default:
    break;
  }
}


//--------------------------------------------------------------------------------------------------
// Update artificial heat
//--------------------------------------------------------------------------------------------------
void Solver::UpdateArtificialHeat()
{
  switch (artificial_heat_)
  {
  case(ArtificialHeatEnum::STANDARD):
  {
    break;
  }
  default:
    break;
  }
}


//--------------------------------------------------------------------------------------------------
// Update artificial stress
//--------------------------------------------------------------------------------------------------
void Solver::UpdateArtificialStress()
{
  switch (artificial_stress_)
  {
  case(ArtificialStressEnum::STANDARD):
  {
    #pragma omp parallel for schedule(guided)
    for (int p = 0; p < particles_num_; p++)
    {
      if (particles_[p].iflag < 0)
      {
        continue;
      }
      for (int n = 0; n < neighbors_[p].size(); n++)
      {
        neighbors_[p][n].zeta = (particles_[p].tension + particles_[neighbors_[p][n].id].tension) * pow4(neighbors_[p][n].kernel_value / avg_kernel_);
      }
    }
    break;
  }
  default:
    break;
  }
}


//--------------------------------------------------------------------------------------------------
// Compute dsdt only
//--------------------------------------------------------------------------------------------------
void Solver::ComputeStressDerivative()
{
  #pragma omp parallel for schedule(guided)
  for (int p = 0; p < particles_num_; p++)
  {
    particles_[p].ddev_stress_dt = Mat3d::Zero();
    if (particles_[p].iflag < 0)
    {
      continue;
    }
    Material& mat = mat_database[parts_[particles_[p].part_id].mat_id];

    // Strain rate & rotate rate
    particles_[p].strain_rate = Mat3d::Zero();
    particles_[p].rotate_rate = Mat3d::Zero();
    for (int n = 0; n < neighbors_[p].size(); n++)
    {
      // Mat3d temp = -neighbors_[p][n].vij * neighbors_[p][n].kernel_gradient.transpose();
      // Mat3d temp = -neighbors_[p][n].vij * (0.5 * (particles_[p].correction + particles_[neighbors_[p][n].id].correction) * neighbors_[p][n].kernel_gradient).transpose();
      Mat3d temp = -neighbors_[p][n].vij * ((particles_[p].correction) * neighbors_[p][n].kernel_gradient).transpose();
      particles_[p].strain_rate += particles_[neighbors_[p][n].id].mass / particles_[neighbors_[p][n].id].rho * (temp + temp.transpose());
      particles_[p].rotate_rate += particles_[neighbors_[p][n].id].mass / particles_[neighbors_[p][n].id].rho * (temp - temp.transpose());
    }
    particles_[p].strain_rate *= 0.5;
    particles_[p].rotate_rate *= 0.5;

    // Elastic-perfactly plastic
    particles_[p].ddev_stress_dt = particles_[p].dev_stress * particles_[p].rotate_rate.transpose() + particles_[p].rotate_rate * particles_[p].dev_stress + 
        2.0 * mat.shear_modulus * particles_[p].alpha_inv * (particles_[p].strain_rate - particles_[p].strain_rate.trace() / 3.0 * Mat3d::Identity());
  }
}


//--------------------------------------------------------------------------------------------------
// Compute dD^(1/3)dt only
//--------------------------------------------------------------------------------------------------
void Solver::ComputeDamageDerivative()
{
  #pragma omp parallel for schedule(guided)
  for (int p = 0; p < particles_num_; p++)
  {
    particles_[p].ddamage_root3_dt = 0.0;
    particles_[p].frag.dAdt = 0.0;
    if (particles_[p].iflag < 0)
    {
      continue;
    }
    Material& mat = mat_database[parts_[particles_[p].part_id].mat_id];

    // Grady kipp damage
    if (mat.damage == DamageEnum::GRADY_KIPP && particles_[p].damage < 1.0 - kEps)
    {
      double max_principle_stress = particles_[p].stress.eigenvalues().real().maxCoeff();
      if (max_principle_stress > 0.0)
      {
        double strain = max_principle_stress / ((1.0 - particles_[p].damage) * mat.elastic_modulus * particles_[p].alpha_inv);

        if (strain >= particles_[p].flaws.strain_min)
        {
          double nact = 0.0;
          if (strain >= particles_[p].flaws.strain_max)
          {
            nact = static_cast<double>(particles_[p].flaws.ntot);
          }
          else
          {
            nact = pow(strain / particles_[p].flaws.strain_min, particles_[p].flaws.m0);
          }
          // double c_growth  = mat.grady_kipp.cg_ce * sqrt((mat.bulk_modulus + 4.0 / 3.0 * mat.shear_modulus * 
          //     (1.0 - particles_[p].damage)) * particles_[p].alpha_inv / particles_[p].rho);
          double c_growth  = mat.grady_kipp.cg_ce * sqrt((mat.bulk_modulus + 4.0 / 3.0 * mat.shear_modulus) * particles_[p].alpha_inv / particles_[p].rho0);

          particles_[p].flaws.damage_root3_max = pow(nact / static_cast<double>(particles_[p].flaws.ntot), 1.0 / 9.0);
          particles_[p].ddamage_root3_dt = nact * c_growth / (kernel_.get_kappa() * particles_[p].hsml);

          // record fracture area rate
          double m_ = mat.grady_kipp.m;
          double alpha = 8.0 * kPI * pow3(c_growth) * mat.grady_kipp.k / ((m_ + 1.0) * (m_ + 2.0) * (m_ + 3.0));
          particles_[p].frag.dAdt = (m_ + 2.0) * (m_ + 3.0) / (2.0*c_growth) * pow(alpha, 2.0/3.0) * pow(strain, m_*2.0/3.0); // * D^(1/3)

          MY_ASSERT(particles_[p].ddamage_root3_dt >= 0.0);
        }
      }
    }
  }
}


//--------------------------------------------------------------------------------------------------
// Compute distention derivative
//--------------------------------------------------------------------------------------------------
void Solver::ComputeDistentionDerivative()
{
  #pragma omp parallel for schedule(guided)
  for (int p = 0; p < particles_num_; p++)
  {
    if (particles_[p].iflag < 0)
    {
      continue;
    }
    Material& mat = mat_database[parts_[particles_[p].part_id].mat_id];
    switch (mat.porosity)
    {
    case(PorosityEnum::P_ALPHA):
    {
      particles_[p].dalpha_dt = 0.0;
      particles_[p].dpdt = 0.0;
      if (mat.eos == EosEnum::TILLOTSON || mat.eos == EosEnum::SIMPLIFIED_TILLOTSON)
      {
        // Partial derivative of eos
        double pdp_drho, pdp_de;
        pdp_drho = mat.sim_tillotson.c * particles_[p].energy + mat.sim_tillotson.A / (particles_[p].rho0 * mat.p_alpha.alpha0);
        pdp_de   = mat.sim_tillotson.c * particles_[p].rho * particles_[p].alpha;

        // Crush curve
        double dalpha_dp = 0.0;
        double p_crush = mat.p_alpha.ps - (mat.p_alpha.ps-mat.p_alpha.pe) * sqrt((particles_[p].alpha-1.0)/(mat.p_alpha.alpha0-1.0));
        // Plastic loading compaction: here dalpha_dt is only used to judge if dPdt > 0
        // if (particles_[p].dalpha_dt > 0.0 && particles_[p].pressure > mat.p_alpha.pe && particles_[p].pressure < mat.p_alpha.ps)
        if (particles_[p].pressure >= p_crush && particles_[p].pressure <= mat.p_alpha.ps)
        {
          dalpha_dp = 2.0 * (mat.p_alpha.alpha0 - 1.0) * (particles_[p].pressure - mat.p_alpha.ps) / pow2(mat.p_alpha.ps - mat.p_alpha.pe);
        }
        else if (particles_[p].pressure > mat.p_alpha.ps)
        {
          // dalpha_dp = (1.0-mat.p_alpha.alpha0) / (particles_[p].pressure - mat.p_alpha.ps);
          dalpha_dp = (1.0-mat.p_alpha.alpha0) / (mat.p_alpha.ps - mat.p_alpha.pe);
        }
        MY_ASSERT(dalpha_dp < kEps);

        // particles_[p].dalpha_dt = dalpha_dp * particles_[p].drhodt * (pdp_de * particles_[p].pressure / pow2(particles_[p].rho) +
        //    particles_[p].alpha * pdp_drho) / (particles_[p].alpha + dalpha_dp * (particles_[p].pressure - particles_[p].rho * pdp_drho));
        particles_[p].dpdt = (particles_[p].dedt * pdp_de + particles_[p].drhodt * particles_[p].alpha * pdp_drho) /
            (particles_[p].alpha + dalpha_dp * (particles_[p].pressure - particles_[p].rho * pdp_drho));
        particles_[p].dalpha_dt = dalpha_dp * particles_[p].dpdt;

        // // for updating dsdt
        // particles_[p].dalpha_drho = dalpha_dp * (particles_[p].pressure / pow2(particles_[p].rho) * pdp_de + 
        //     particles_[p].alpha * pdp_drho) / (particles_[p].alpha + dalpha_dp * (particles_[p].pressure - particles_[p].rho * pdp_drho));

        MY_ASSERT(std::isfinite(particles_[p].dalpha_dt));
        particles_[p].dalpha_dt = min(particles_[p].dalpha_dt, 0.0);
      }
      else
      {
        out_.PrintLogger("Error: p-alpha model is only implemented to tillotson eos!\n");
      }
      
      break;
    }
    default:
      break;
    }
  }
}


//--------------------------------------------------------------------------------------------------
// Correction
//--------------------------------------------------------------------------------------------------
void Solver::Correction()
{
  // update correction tensor
  #pragma omp parallel for schedule(guided)
  for (int p = 0; p < particles_num_; p++)
  {
    if (particles_[p].iflag < 0 || neighbors_[p].size() == 0)
    {
      particles_[p].correction = Mat3d::Identity();
      continue;
    }
    particles_[p].correction = Mat3d::Zero();
    double completeness = 0.0;
    for (int n = 0; n < neighbors_[p].size(); n++)
    {
      particles_[p].correction += particles_[neighbors_[p][n].id].mass / particles_[neighbors_[p][n].id].rho *
          neighbors_[p][n].kernel_gradient * (-neighbors_[p][n].xij).transpose();
      completeness += particles_[neighbors_[p][n].id].mass / particles_[neighbors_[p][n].id].rho * neighbors_[p][n].kernel_value;
    }

    // // compute the eigenvalues of correction tensor
    // double lmd1 = particles_[p].correction.eigenvalues().real().maxCoeff();
    // double lmd3 = particles_[p].correction.eigenvalues().real().minCoeff();
    // double ratio = lmd1 / lmd3;
    
    // Blend with identity matrix based on completeness
    completeness = middleOfThree((completeness-0.6)*2.5, 0.0, 1.0);
    // if (ratio > 1.0e4)
    // {
    //   if (ratio < 1.0e8)
    //     completeness *= (1.0e8 - ratio) * 1.0e-4 / 9999.0;
    //   else
    //     completeness = 0.0;
    // }

    // Compute the SVD of the matrix
    particles_[p].correction = pseudoInverse(particles_[p].correction);
    particles_[p].correction = completeness * particles_[p].correction + (1.0 - completeness) * Mat3d::Identity();
  }

  // // correct the kernel gradient
  // #pragma omp parallel for schedule(guided)
  // for (int p = 0; p < particles_num_; p++)
  // {
  //   if (particles_[p].iflag < 0 || neighbors_[p].size() == 0)
  //   {
  //     continue;
  //   }
  //   for (int n = 0; n < neighbors_[p].size(); n++)
  //   {
  //     neighbors_[p][n].kernel_gradient = particles_[p].correction * neighbors_[p][n].kernel_gradient;
  //   }
  // }
  
}


//--------------------------------------------------------------------------------------------------
// Compute density derivative
//--------------------------------------------------------------------------------------------------
void Solver::ComputeDrhoDt()
{
  switch (density_update_)
  {
  case(DensityEnum::CONTINUOUS_DENSITY):
  {
    // density derivative
    #pragma omp parallel for schedule(guided)
    for (int p = 0; p < particles_num_; p++)
    {
      particles_[p].drhodt = 0.0;
      if (particles_[p].iflag < 0)
      {
        continue;
      }
      for (int n = 0; n < neighbors_[p].size(); n++)
      {
        // particles_[p].drhodt += particles_[neighbors_[p][n].id].mass * neighbors_[p][n].vij.dot(neighbors_[p][n].kernel_gradient);
      //   particles_[p].drhodt += particles_[neighbors_[p][n].id].mass / particles_[neighbors_[p][n].id].rho *
      //       neighbors_[p][n].vij.dot(0.5 * (particles_[p].correction + particles_[neighbors_[p][n].id].correction) * neighbors_[p][n].kernel_gradient);
        particles_[p].drhodt += particles_[neighbors_[p][n].id].mass / particles_[neighbors_[p][n].id].rho *
            neighbors_[p][n].vij.dot((particles_[p].correction) * neighbors_[p][n].kernel_gradient);
      }
      particles_[p].drhodt *= particles_[p].rho;

      MY_ASSERT(std::isfinite(particles_[p].drhodt));
    }
    break;
  }

  case(DensityEnum::DENSITY_SUMMATION):
  {
    // create rho_tmp
    double* rho_tmp = new double[particles_num_];

    #pragma omp parallel for schedule(guided)
    for (int p = 0; p < particles_num_; p++)
    {
      if (particles_[p].iflag < 0)
      {
        // particles_[p].rho = particles_[p].mass * kernel_.KernelValue(0.0, particles_[p].hsml);
        continue;
      }
      // A regularized version of density summation
      double rho_sum = particles_[p].mass * kernel_.KernelValue(0.0, particles_[p].hsml);
      double regular = particles_[p].mass / particles_[p].rho * kernel_.KernelValue(0.0, particles_[p].hsml);
      for (int n = 0; n < neighbors_[p].size(); n++)
      {
        rho_sum += particles_[neighbors_[p][n].id].mass * neighbors_[p][n].kernel_value;
        regular += particles_[neighbors_[p][n].id].mass / particles_[neighbors_[p][n].id].rho * neighbors_[p][n].kernel_value;
      }
      rho_tmp[p] = rho_sum / regular; // this can modify neighbor's when looping
    }
        #pragma omp parallel for schedule(guided)
    for (int p = 0; p < particles_num_; p++)
    {
      if (particles_[p].iflag < 0)
      {
        continue;
      }
      particles_[p].rho = rho_tmp[p];
    }
    delete[] rho_tmp;
    break;
  }

  default:
    break;
  }
}


//--------------------------------------------------------------------------------------------------
// Compute velocity derivative
//--------------------------------------------------------------------------------------------------
void Solver::ComputeDvDt()
{
  switch (calc_gravity_)
  {
  case (CalcGravityEnum::SELF_GRAVITY):
  {
  #pragma omp parallel for schedule(guided)
  for (int p = 0; p < particles_num_; p++)
  {
    particles_[p].dvdt = Vec3d::Zero();
    if (particles_[p].iflag < 0)
    {
      continue;
    };
    // normal case
    particles_[p].dvdt += particles_[p].grav;
    for (int n = 0; n < neighbors_[p].size(); n++)
    {
      if (neighbors_[p][n].iflag == 1)
      {
        // // standard version
        // particles_[p].dvdt += particles_[neighbors_[p][n].id].mass * ((particles_[p].stress / pow2(particles_[p].rho) +
        //     particles_[neighbors_[p][n].id].stress / pow2(particles_[neighbors_[p][n].id].rho) - neighbors_[p][n].IIij *
        //     Mat3d::Identity() + neighbors_[p][n].zeta) * (0.5 * (particles_[p].correction + particles_[neighbors_[p][n].id].correction) * neighbors_[p][n].kernel_gradient));

        // Geometric density average version
        // particles_[p].dvdt += particles_[neighbors_[p][n].id].mass * (((particles_[p].stress + particles_[neighbors_[p][n].id].stress) /
        //     (particles_[neighbors_[p][n].id].rho * particles_[p].rho) - neighbors_[p][n].IIij * Mat3d::Identity() + 
        //     neighbors_[p][n].zeta) * (0.5 * (particles_[p].correction + particles_[neighbors_[p][n].id].correction) * neighbors_[p][n].kernel_gradient));
        particles_[p].dvdt += particles_[neighbors_[p][n].id].mass * (((particles_[p].stress + particles_[neighbors_[p][n].id].stress) /
            (particles_[neighbors_[p][n].id].rho * particles_[p].rho) - neighbors_[p][n].IIij * Mat3d::Identity() + 
            neighbors_[p][n].zeta) * ((particles_[p].correction) * neighbors_[p][n].kernel_gradient));
      }
      else
      {
        // Modify ghost neighbor stress
        Mat3d ghost_stress = 2.0 * (Mat3d::Identity().array() * particles_[p].stress.array()).matrix() - particles_[p].stress;
      //   particles_[p].dvdt += particles_[neighbors_[p][n].id].mass * ((particles_[p].stress / pow2(particles_[p].rho) +
      //       ghost_stress / pow2(particles_[neighbors_[p][n].id].rho) - neighbors_[p][n].IIij *
      //       Mat3d::Identity() + neighbors_[p][n].zeta) * (0.5 * (particles_[p].correction + particles_[neighbors_[p][n].id].correction) * neighbors_[p][n].kernel_gradient));
        particles_[p].dvdt += particles_[neighbors_[p][n].id].mass * ((particles_[p].stress / pow2(particles_[p].rho) +
            ghost_stress / pow2(particles_[neighbors_[p][n].id].rho) - neighbors_[p][n].IIij *
            Mat3d::Identity() + neighbors_[p][n].zeta) * ((particles_[p].correction) * neighbors_[p][n].kernel_gradient));
      }
    }
    MY_ASSERT(std::isfinite(particles_[p].dvdt(0)));
    MY_ASSERT(std::isfinite(particles_[p].dvdt(1)));
    MY_ASSERT(std::isfinite(particles_[p].dvdt(2)));
  }
    break;
  }

  default: // no gravity
  {
  #pragma omp parallel for schedule(guided)
  for (int p = 0; p < particles_num_; p++)
  {
    particles_[p].dvdt = Vec3d::Zero();
    if (particles_[p].iflag < 0)
    {
      continue;
    }
    for (int n = 0; n < neighbors_[p].size(); n++)
    {
      if (neighbors_[p][n].iflag == 1)
      {
        particles_[p].dvdt += particles_[neighbors_[p][n].id].mass * ((particles_[p].stress / pow2(particles_[p].rho) +
            particles_[neighbors_[p][n].id].stress / pow2(particles_[neighbors_[p][n].id].rho) - neighbors_[p][n].IIij *
            Mat3d::Identity() + neighbors_[p][n].zeta) * neighbors_[p][n].kernel_gradient);
      }
      else
      {
        // Modify ghost neighbor stress
        Mat3d ghost_stress = 2.0 * (Mat3d::Identity().array() * particles_[p].stress.array()).matrix() - particles_[p].stress;
        // particles_[p].dvdt += particles_[neighbors_[p][n].id].mass * ((particles_[p].stress / pow2(particles_[p].rho) +
        //     ghost_stress / pow2(particles_[neighbors_[p][n].id].rho) - neighbors_[p][n].IIij *
        //     Mat3d::Identity() + neighbors_[p][n].zeta) * (0.5 * (particles_[p].correction + particles_[neighbors_[p][n].id].correction) * neighbors_[p][n].kernel_gradient));
        particles_[p].dvdt += particles_[neighbors_[p][n].id].mass * ((particles_[p].stress / pow2(particles_[p].rho) +
            ghost_stress / pow2(particles_[neighbors_[p][n].id].rho) - neighbors_[p][n].IIij *
            Mat3d::Identity() + neighbors_[p][n].zeta) * ((particles_[p].correction) * neighbors_[p][n].kernel_gradient));
      }
    }
    MY_ASSERT(std::isfinite(particles_[p].dvdt(0)));
    MY_ASSERT(std::isfinite(particles_[p].dvdt(1)));
    MY_ASSERT(std::isfinite(particles_[p].dvdt(2)));
  }
    break;
  }
  }
}


//--------------------------------------------------------------------------------------------------
// Compute energy derivative
//--------------------------------------------------------------------------------------------------
void Solver::ComputeDeDt()
{
  #pragma omp parallel for schedule(guided)
  for (int p = 0; p < particles_num_; p++)
  {
    particles_[p].dedt = 0.0;
    if (particles_[p].iflag < 0)
    {
      continue;
    }
    // // standard version
    // for (int n = 0; n < neighbors_[p].size(); n++)
    // {
    //   particles_[p].dedt += particles_[neighbors_[p][n].id].mass * (particles_[p].pressure / pow2(particles_[p].rho) +
    //       particles_[neighbors_[p][n].id].pressure / pow2(particles_[neighbors_[p][n].id].rho) + neighbors_[p][n].IIij) *
    //       neighbors_[p][n].vij.dot((0.5 * (particles_[p].correction + particles_[neighbors_[p][n].id].correction) * neighbors_[p][n].kernel_gradient));
    // }
    // particles_[p].dedt *= 0.5;

    // Geometric density average version
    for (int n = 0; n < neighbors_[p].size(); n++)
    {
      // particles_[p].dedt += particles_[neighbors_[p][n].id].mass * (particles_[p].pressure / (particles_[p].rho * particles_[neighbors_[p][n].id].rho) +
      //     neighbors_[p][n].IIij) *
      //     neighbors_[p][n].vij.dot((0.5 * (particles_[p].correction + particles_[neighbors_[p][n].id].correction) * neighbors_[p][n].kernel_gradient));
        particles_[p].dedt += particles_[neighbors_[p][n].id].mass * (particles_[p].pressure / (particles_[p].rho * particles_[neighbors_[p][n].id].rho) +
          neighbors_[p][n].IIij) * neighbors_[p][n].vij.dot(((particles_[p].correction) * neighbors_[p][n].kernel_gradient));
      }
    
    particles_[p].dedt += (particles_[p].dev_stress.array() * particles_[p].strain_rate.array()).sum() / particles_[p].rho + 
        particles_[p].Hi;

    MY_ASSERT(std::isfinite(particles_[p].dedt));
  }
}


//--------------------------------------------------------------------------------------------------
// Compute position derivative
//--------------------------------------------------------------------------------------------------
void Solver::ComputeDxDt()
{
  switch (xsph_)
  {
  case(XsphEnum::STANDARD):
  {
    #pragma omp parallel for schedule(guided)
    for (int p = 0; p < particles_num_; p++)
    {
      particles_[p].dxdt = particles_[p].vel;
      if (particles_[p].iflag < 0)
      {
        continue;
      }
      for (int n = 0; n < neighbors_[p].size(); n++)
      {
        particles_[p].dxdt -= particles_[neighbors_[p][n].id].mass / (particles_[p].rho + particles_[neighbors_[p][n].id].rho) *
          neighbors_[p][n].vij * neighbors_[p][n].kernel_value;
      }

      MY_ASSERT(std::isfinite(particles_[p].dxdt(0)));
      MY_ASSERT(std::isfinite(particles_[p].dxdt(1)));
      MY_ASSERT(std::isfinite(particles_[p].dxdt(2)));
    }
    break;
  }
  case(XsphEnum::NONE):
  {
    #pragma omp parallel for schedule(guided)
    for (int p = 0; p < particles_num_; p++)
    {
      particles_[p].dxdt = particles_[p].vel;
    }
    break;
  }
  default:
    break;
  }

}


//--------------------------------------------------------------------------------------------------
// Update smooth length to make sure proper neighbors size
//--------------------------------------------------------------------------------------------------
void Solver::UpdateHsml()
{
  current_hmax_ = hsml_min_;

  #pragma omp parallel for schedule(guided) reduction(max:current_hmax_)
  for (int p = 0; p < particles_num_; p++)
  {
    if (particles_[p].iflag < 0)
    {
      // Particle with iflag < 0 does not change hsml as well as current hmax
      continue;
    }

    // Update smooth length differently depending on the density update
    switch (density_update_)
    {
    case(DensityEnum::CONTINUOUS_DENSITY):
    {
      particles_[p].hsml -= dt_ * particles_[p].hsml * particles_[p].drhodt / (3.0 * particles_[p].rho);
      break;
    }
    case(DensityEnum::DENSITY_SUMMATION):
    {
      particles_[p].hsml = kEta * pow(particles_[p].mass / particles_[p].rho, 1.0 / 3.0);
      break;
    }
    default:
      break;
    }
    
    // Limit smooth length
    particles_[p].hsml = max(hsml_min_, particles_[p].hsml);
    particles_[p].hsml = min(hsml_max_, particles_[p].hsml);

    // Current max hsml for neighbor search
    current_hmax_ = max(current_hmax_, particles_[p].hsml);

    MY_ASSERT(std::isfinite(particles_[p].hsml));
    MY_ASSERT(particles_[p].hsml > 0.0);
  }
}


//--------------------------------------------------------------------------------------------------
// Update current internal / kinetic energy
//--------------------------------------------------------------------------------------------------
void Solver::UpdateCurrentEnergy()
{
  current_internal_energy_ = 0.0;
  current_kinetic_energy_  = 0.0;
  for (int p = 0; p < particles_num_; p++)
  {
    current_internal_energy_ += particles_[p].mass * particles_[p].energy;
    current_kinetic_energy_  += particles_[p].mass * particles_[p].vel.dot(particles_[p].vel);
  }
  current_kinetic_energy_ *= 0.5;
}

