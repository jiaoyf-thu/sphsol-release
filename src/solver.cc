#include"solver.h"


//--------------------------------------------------------------------------------------------------
// Run the simulation
//--------------------------------------------------------------------------------------------------
void Solver::Run()
{
  timeval start;
  gettimeofday(&start, NULL);

  for (const auto& mode : run_modes_)
  {
    current_run_mode_ = mode;
    switch (current_run_mode_)
    {
      case (RunMode::RELAX):
      {
        Relax();
        break;
      }
      case (RunMode::SHOCK):
      {
        ShockRun();
        break;
      }
      case (RunMode::LATERUN):
      {
        LateRun();
        break;
      }

      default:
        break;
    }
  }

  out_.PrintLogger("  Calculation Time = %f s\n\n", GetTimeInterval(&start));
};



//--------------------------------------------------------------------------------------------------
// fast intergrate to a stable state
//--------------------------------------------------------------------------------------------------
void Solver::Relax()
{
  if ((calc_gravity_ == CalcGravityEnum::NONE) && (is_rotating_ == false))
  {
    out_.PrintLogger(" SKIP RELAX.\n\n");
    return;
  }

  out_.PrintLogger(" RELAX TARGET BODY:\n\n");

  int temp_pnum = particles_num_; // for the whole body!

  // particles_num_ = parts_[0].particle_id.size(); // for the target only
  // particles_num_ -= parts_[parts_num_-1].particle_id.size();

  BuildNeighbors(calc_gravity_ == CalcGravityEnum::SELF_GRAVITY);

  // main loop
  switch (integrator_)
  {
    // leapfrog integration
    case (IntegratorEnum::LEAP_FROG):
    {
      timeval start;
      gettimeofday(&start, NULL);

      while (current_time_ < damp_time_)
      {
        #pragma omp parallel for schedule(guided)
        for (int p = 0; p < particles_num_; p++)
        {
          particles_[p].vel += 0.5*dt_*particles_[p].dvdt;
          particles_[p].pos += dt_*particles_[p].vel;
        }

        // UpdateNeighbors();

        // if (dt_stp_ >= dt_ref_)
        // {
        //   dt_stp_ = 0.0;
        //   UpdateGravity(); // after UpdateNeighbors
        // }

        // General case for shock simulation
        if (dt_stp_ >= dt_ref_)
        {
          dt_stp_ = 0.0;
          BuildNeighbors(calc_gravity_ == CalcGravityEnum::SELF_GRAVITY); // Update neighbors & gravity
        }
        else
        {
          UpdateNeighbors(); // Update neighbors only
        }

        // artificial terms
        UpdateArtificialViscosity();
        // UpdateArtificialHeat();
        UpdateArtificialStress();

        // Time derivative for integration
        Correction();
        ComputeDrhoDt();
        ComputeDeDt();
        ComputeStressDerivative();

        // // This is a previous version for elastic material only
        // #pragma omp parallel for schedule(guided) reduction(max:current_hmax_)
        // for (int p = 0; p < particles_num_; p++)
        // {
        //   double rho_frac = 1.0 + dt_ * particles_[p].drhodt / particles_[p].rho;
        //   particles_[p].rho *= rho_frac;
        //   particles_[p].mass *= rho_frac;  // also modify the particle mass here
        //   particles_[p].energy += dt_ * particles_[p].dedt;

        //   // pressure update: assuming elastic
        //   Material& mat = mat_database[parts_[particles_[p].part_id].mat_id];
        //   double alpha = (mat.porosity == PorosityEnum::P_ALPHA) ? particles_[p].alpha : 1.0;
        //   double alpha0 = (mat.porosity == PorosityEnum::P_ALPHA) ? mat.p_alpha.alpha0 : 1.0;
        //   double mu = (particles_[p].rho * alpha) / (particles_[p].rho0 * alpha0) - 1.0;
        //   double pressure = mat.bulk_modulus * particles_[p].alpha_inv * mu;
        //   particles_[p].pressure = pressure < 0.0 ? pressure * (1.0 - particles_[p].damage) : pressure;
        //   particles_[p].c_sound = sqrt(mat.bulk_modulus * particles_[p].alpha_inv / particles_[p].rho);

        //   // stress update
        //   Mat3d dev_stress_temp = particles_[p].dev_stress + particles_[p].ddev_stress_dt * dt_;
        //   double J2 = 0.5 * (dev_stress_temp.array() * dev_stress_temp.array()).sum();
        //   double f = StrengthModel(mat, particles_[p].damage, J2, particles_[p].pressure);
        //   particles_[p].dev_stress = f * dev_stress_temp;
        //   particles_[p].stress = particles_[p].dev_stress - particles_[p].pressure * Mat3d::Identity();
        // }

        // This is the same as the shock run, except damage, distention
        #pragma omp parallel for schedule(guided) reduction(max:current_hmax_)
        for (int p = 0; p < particles_num_; p++)
        {
          Material& mat = mat_database[parts_[particles_[p].part_id].mat_id];

          // Density evolution limit
          // double rho_limit = particles_[p].mass * kernel_.KernelValue(0.0, hsml_ref_);
          double rho_limit = 0.8 * particles_[p].rho0;
          particles_[p].rho = particles_[p].rho + dt_ * particles_[p].drhodt;
          if (particles_[p].rho < rho_limit)
          {
            particles_[p].rho = rho_limit;
            // particles_[p].hsml = hsml_ref_;
          }
          else
          {
            // hsml update
            particles_[p].hsml *= 1.0 - dt_ * particles_[p].drhodt / (3.0 * particles_[p].rho);
            particles_[p].hsml = middleOfThree(particles_[p].hsml, hsml_min_, hsml_max_);
          }
          current_hmax_ = max(current_hmax_, particles_[p].hsml);

          // energy update
          particles_[p].energy += dt_ * particles_[p].dedt;
          // peak thermal energy update only for shock phase
          particles_[p].energy_cold += dt_ * EquatonOfStateCold(&particles_[p], mat) / pow2(particles_[p].rho) * particles_[p].drhodt;
          particles_[p].energy_thermal_peak = max(particles_[p].energy_thermal_peak, particles_[p].energy-particles_[p].energy_cold);

          // eos
          EquatonOfState(&particles_[p], mat);
          particles_[p].pressure_peak = max(particles_[p].pressure_peak, particles_[p].pressure);

          // stress update
          Mat3d dev_stress_temp = particles_[p].dev_stress + particles_[p].ddev_stress_dt * dt_;
          double J2 = 0.5 * (dev_stress_temp.array() * dev_stress_temp.array()).sum();
          double f = StrengthModel(mat, particles_[p].damage, J2, particles_[p].pressure);
          particles_[p].dev_stress = f * dev_stress_temp;
          particles_[p].stress = particles_[p].dev_stress - particles_[p].pressure * Mat3d::Identity();

          // tension stress update
          if (artificial_stress_ == ArtificialStressEnum::STANDARD)
          {
            // eigrnvalue and eigenvector of stress
            Eigen::SelfAdjointEigenSolver<Mat3d> es(particles_[p].stress);
            Vec3d eigval = es.eigenvalues();
            Mat3d eigvec = es.eigenvectors();
            for (int i = 0; i < 3; i++)
            {
              eigval(i) = max(eigval(i), 0.0);
            }
            particles_[p].tension = -0.2 / pow2(particles_[p].rho) * eigvec * eigval.asDiagonal() * eigvec.transpose();
          }
        }

        ComputeDvDt();
        UpdateTimeStep();


        // damping term
        double damping = (1.0 - 0.99 * current_time_ / damp_time_) / dt_;

        if (is_rotating_)
        {
          #pragma omp parallel for schedule(guided)
          for (int p = 0; p < particles_num_; p++)
          {
            Vec3d r_rel = particles_[p].pos - rotation_.pos0;
            Vec3d vel_rigid = rotation_.omg0.cross(r_rel);
            
            particles_[p].dvdt -= damping * (particles_[p].vel - vel_rigid); // damping term
            particles_[p].vel += 0.5*dt_*particles_[p].dvdt;
          }
        }
        else
        {
          #pragma omp parallel for schedule(guided)
          for (int p = 0; p < particles_num_; p++)
          {
            // apply damping
            particles_[p].dvdt -= damping * particles_[p].vel;
            particles_[p].vel += 0.5*dt_*particles_[p].dvdt;
          }
        }

        current_step_++;
        current_time_ += dt_;
        dt_stp_ += dt_;
        
        CheckStep(&start, damp_time_);
      }
      break;
    }

    case (IntegratorEnum::PREDICT_CORRECT):
    {
      out_.PrintLogger("  - Predict-Correct Integrator Not Implemented Yet!\n\n");
      break;
    }

    default:
      break;
  }
  
  // reset
  particles_num_ = temp_pnum;
  current_step_ = 0;
  current_time_ = 0.0;
  dt_stp_ = 0.0;
  dt_ = dt_min_;
};


//--------------------------------------------------------------------------------------------------
// SPH run for shock simulation
//--------------------------------------------------------------------------------------------------
void Solver::ShockRun()
{
  out_.PrintLogger(" RUN SHOCK STAGE:\n\n");

  BuildNeighbors(calc_gravity_ == CalcGravityEnum::SELF_GRAVITY); // Update neighbors & gravity

  switch (integrator_)
  {
    case (IntegratorEnum::LEAP_FROG):
    {
      timeval start;
      gettimeofday(&start, NULL);

      while (current_time_ < total_time_ * (1 + kEps))
      {
        // leapfrog integration
        #pragma omp parallel for schedule(guided)
        for (int p = 0; p < particles_num_; p++)
        {
          particles_[p].vel += 0.5*dt_*particles_[p].dvdt;
          particles_[p].pos += dt_*particles_[p].vel;
        }
        
        // General case for shock simulation
        if (dt_stp_ >= dt_ref_)
        // // Tidal case, we use a much less frequent neighbor/gravity update
        // if (dt_stp_ >= dt_ref_*100.0)
        {
          dt_stp_ = 0.0;
          BuildNeighbors(calc_gravity_ == CalcGravityEnum::SELF_GRAVITY); // Update neighbors & gravity

          // Add force if needed
          if (additional_force_) { AddForce(); }
        }
        else
        {
          UpdateNeighbors(); // Update neighbors only
        }

        // artificial terms
        UpdateArtificialViscosity();
        UpdateArtificialHeat();
        UpdateArtificialStress();

        // Time derivative for integration
        Correction();
        ComputeDrhoDt();
        ComputeDeDt();
        ComputeStressDerivative();
        ComputeDamageDerivative();
        ComputeDistentionDerivative();

        #pragma omp parallel for schedule(guided) reduction(max:current_hmax_)
        for (int p = 0; p < particles_num_; p++)
        {
          Material& mat = mat_database[parts_[particles_[p].part_id].mat_id];

          // Density evolution limit
          switch (density_update_)
          {
          case(DensityEnum::CONTINUOUS_DENSITY):
          {
            double rho_limit = 0.8 * particles_[p].rho0;
            particles_[p].rho = particles_[p].rho + dt_ * particles_[p].drhodt;
            if (particles_[p].rho < rho_limit)
            {
              particles_[p].rho = rho_limit;
              // particles_[p].hsml = hsml_ref_;
            }
            else
            {
              // hsml update
              particles_[p].hsml *= 1.0 - dt_ * particles_[p].drhodt / (3.0 * particles_[p].rho);
              particles_[p].hsml = middleOfThree(particles_[p].hsml, hsml_min_, hsml_max_);
            }
            break;
          }
          case(DensityEnum::DENSITY_SUMMATION):
          {
            particles_[p].hsml = kEta * pow(particles_[p].mass / particles_[p].rho, 1.0/3.0);
            particles_[p].hsml = middleOfThree(particles_[p].hsml, hsml_min_, hsml_max_);
            break;
          }
          default:
            break;
          }
          current_hmax_ = max(current_hmax_, particles_[p].hsml);

          // energy update
          particles_[p].energy += dt_ * particles_[p].dedt;
          // peak thermal energy update only for shock phase
          particles_[p].energy_cold += dt_ * EquatonOfStateCold(&particles_[p], mat) / pow2(particles_[p].rho) * particles_[p].drhodt;
          particles_[p].energy_thermal_peak = max(particles_[p].energy_thermal_peak, particles_[p].energy-particles_[p].energy_cold);

          // eos
          EquatonOfState(&particles_[p], mat);
          particles_[p].pressure_peak = max(particles_[p].pressure_peak, particles_[p].pressure);

          // Distention update
          if (mat.porosity == PorosityEnum::P_ALPHA)
          {
            particles_[p].alpha = max(particles_[p].alpha + dt_ * particles_[p].dalpha_dt, 1.0);
            particles_[p].alpha_inv = 1.0 / particles_[p].alpha;
            // particles_[p].ddev_stress_dt *= 1.0 + particles_[p].dalpha_drho * particles_[p].rho / particles_[p].alpha;
          }

          // stress update
          Mat3d dev_stress_temp = particles_[p].dev_stress + particles_[p].ddev_stress_dt * dt_;
          double J2 = 0.5 * (dev_stress_temp.array() * dev_stress_temp.array()).sum();
          double f = StrengthModel(mat, particles_[p].damage, J2, particles_[p].pressure, particles_[p].plastic_strain);
          particles_[p].dev_stress = f * dev_stress_temp;
          particles_[p].stress = particles_[p].dev_stress - particles_[p].pressure * Mat3d::Identity();

          // plastic strain accumulate
          particles_[p].plastic_strain += (1.0-f) * sqrt(J2/3.0) / mat.shear_modulus;

          // Damage accumulate
          if (mat.damage == DamageEnum::GRADY_KIPP)
          {
            if (particles_[p].damage < 1.0-kEps)
            {
              particles_[p].damage_root3 = middleOfThree(particles_[p].damage_root3 + particles_[p].ddamage_root3_dt * dt_,
                    particles_[p].flaws.damage_root3_max, particles_[p].damage_root3);
              particles_[p].damage = min(1.0, pow3(particles_[p].damage_root3));

              // record damage process to estimate fragment size
              if (particles_[p].damage > kEps && particles_[p].damage < 1.0-kEps)
              {
                // particles_[p].frag.tf += dt_;
                particles_[p].frag.dAdt *= particles_[p].damage_root3;
                particles_[p].frag.area += particles_[p].frag.dAdt * dt_;
              }
            }
            particles_[p].frag.Lm = 3.0 * (mat.grady_kipp.m+3.0) / (mat.grady_kipp.m+2.0) / particles_[p].frag.area;
          }

          // tension stress update
          if (artificial_stress_ == ArtificialStressEnum::STANDARD)
          {
            // eigrnvalue and eigenvector of stress
            Eigen::SelfAdjointEigenSolver<Mat3d> es(particles_[p].stress);
            Vec3d eigval = es.eigenvalues();
            Mat3d eigvec = es.eigenvectors();
            for (int i = 0; i < 3; i++)
            {
              eigval(i) = max(eigval(i), 0.0);
            }
            particles_[p].tension = -0.2 / pow2(particles_[p].rho) * eigvec * eigval.asDiagonal() * eigvec.transpose();
          }
        }

        ComputeDvDt();
        UpdateTimeStep();

        #pragma omp parallel for schedule(guided)
        for (int p = 0; p < particles_num_; p++)
        {
          particles_[p].vel += 0.5*dt_*particles_[p].dvdt;
        }

        current_time_ += dt_;
        current_step_ += 1;
        dt_stp_ += dt_;

        CheckStep(&start, total_time_);
      }
      break;
    }

    case (IntegratorEnum::PREDICT_CORRECT):
    {
      out_.PrintLogger("  - Predict-Correct Integrator Not Implemented Yet!\n\n");
      break;
    }
  
    default:
      break;
  }

  // reset
  current_step_ = 0;
  dt_stp_ = 0.0;
  dt_ = dt_min_;
}


//--------------------------------------------------------------------------------------------------
// Late run
//--------------------------------------------------------------------------------------------------
void Solver::LateRun()
{
  out_.PrintLogger(" RUN LATE STAGE:\n\n");

  BuildNeighbors(calc_gravity_ == CalcGravityEnum::SELF_GRAVITY); // Update neighbors & gravity

  switch (integrator_)
  {
    case (IntegratorEnum::LEAP_FROG):
    {
      timeval start;
      gettimeofday(&start, NULL);

      while (current_time_ < late_time_ * (1 + kEps))
      {
        // leapfrog integration
        #pragma omp parallel for schedule(guided)
        for (int p = 0; p < particles_num_; p++)
        {
          particles_[p].vel += 0.5*dt_*particles_[p].dvdt;
          particles_[p].pos += dt_*particles_[p].vel;
        }
        
        if (dt_stp_ >= dt_ref_ * 100.0)
        {
          dt_stp_ = 0.0;
          BuildNeighbors(calc_gravity_ == CalcGravityEnum::SELF_GRAVITY); // Update neighbors & gravity every 10 steps

          // Add force if needed
          if (additional_force_) { AddForce(); }
        }
        else
        {
          UpdateNeighbors(); // Update neighbors only
        }

        // artificial terms
        UpdateArtificialViscosity();

        // Time derivative for integration
        Correction();
        ComputeDrhoDt();
        ComputeDeDt();
        ComputeStressDerivative();

        #pragma omp parallel for schedule(guided) reduction(max:current_hmax_)
        for (int p = 0; p < particles_num_; p++)
        {
          // Density evolution limit
          double rho_limit = 0.8 * particles_[p].rho0;
          particles_[p].rho = particles_[p].rho + dt_ * particles_[p].drhodt;
          if (particles_[p].rho < rho_limit)
          {
            particles_[p].rho = rho_limit;
            // particles_[p].hsml = hsml_ref_;
          }
          else
          {
            // hsml update
            particles_[p].hsml *= 1.0 - dt_ * particles_[p].drhodt / (3.0 * particles_[p].rho);
            particles_[p].hsml = middleOfThree(particles_[p].hsml, hsml_min_, hsml_max_);
          }
          current_hmax_ = max(current_hmax_, particles_[p].hsml);

          // energy update
          particles_[p].energy += dt_ * particles_[p].dedt;

          // // eos
          // Material& mat = mat_database[parts_[particles_[p].part_id].mat_id];
          // EquatonOfState(&particles_[p], mat);
          // pressure update: assuming elastic
          Material& mat = mat_database[parts_[particles_[p].part_id].mat_id];
          double alpha = (mat.porosity == PorosityEnum::P_ALPHA) ? particles_[p].alpha : 1.0;
          double alpha0 = (mat.porosity == PorosityEnum::P_ALPHA) ? mat.p_alpha.alpha0 : 1.0;
          double mu = (particles_[p].rho * alpha) / (particles_[p].rho0 * alpha0) - 1.0;
          double pressure = mat.bulk_modulus * particles_[p].alpha_inv * mu;
          double pmin = (mat.yield == YieldEnum::LUND) ? -mat.lund.Yd0 * max(0.0, 1.0 - particles_[p].plastic_strain) / mat.lund.mud : 0.0;

          // particles_[p].pressure = pressure < 0.0 ? pressure * (1.0 - particles_[p].damage) : pressure;
          particles_[p].pressure = max(pressure, pmin);
          particles_[p].c_sound  = sqrt(mat.bulk_modulus * particles_[p].alpha_inv / particles_[p].rho);

          // stress update
          Mat3d dev_stress_temp = particles_[p].dev_stress + particles_[p].ddev_stress_dt * dt_;
          double J2 = 0.5 * (dev_stress_temp.array() * dev_stress_temp.array()).sum();
          double f = StrengthModel(mat, particles_[p].damage, J2, particles_[p].pressure, particles_[p].plastic_strain);
          particles_[p].dev_stress = f * dev_stress_temp;
          particles_[p].stress = particles_[p].dev_stress - particles_[p].pressure * Mat3d::Identity();

          // assume no strain increase in the late stage
        }

        ComputeDvDt();
        UpdateTimeStep();

        #pragma omp parallel for schedule(guided)
        for (int p = 0; p < particles_num_; p++)
        {
          particles_[p].vel += 0.5*dt_*particles_[p].dvdt;
        }

        current_time_ += dt_;
        current_step_ += 1;
        dt_stp_ += dt_;

        CheckStep(&start, late_time_);
      }
      break;
    }

    case (IntegratorEnum::PREDICT_CORRECT):
    {
      out_.PrintLogger("  - Predict-Correct Integrator Not Implemented Yet!\n\n");
      break;
    }
  
    default:
      break;
  }
}


//--------------------------------------------------------------------------------------------------
// Check the timestep for output
//--------------------------------------------------------------------------------------------------
void Solver::CheckStep(timeval* start, double total_time)
{
  // Export the time and energy for reference only
  if (current_step_ == 1 || current_step_ % 1000 == 0)
  {
    Vec3i real_time = GetCurrentTime();
    out_.PrintLogger("  - Current Time Step = %d, dt = %e, Real Time = %02d:%02d:%02d\n", current_step_, dt_, 
        real_time(0), real_time(1), real_time(2));

    UpdateCurrentEnergy();
    out_.PrintLogger("  - Total Energy = %e (Internal: %e, Kinetic: %e)\n\n", current_internal_energy_ + current_kinetic_energy_,
        current_internal_energy_, current_kinetic_energy_);
  }

  // export particles data only for shock and late run
  if (current_run_mode_ != RunMode::RELAX)
  {
    if (current_time_ >= out_.output_time_[out_.output_step_])
    {
    out_.PrintParticles(particles_, particles_num_);
    out_.PrintLogger("  - Write File: Simulation Time = %e, Calculation Time = %f s\n\n", current_time_, GetTimeInterval(start, true));
    }  
    // Estimated calculation time
    if (current_step_ == 1)
    {
      double tmp = GetTimeInterval(start, true); // Refresh time record because the first step costs more time
    }
    if (current_step_ == 2)
    {
      out_.PrintLogger("  Estimated Calculation Time = %f s\n\n", GetTimeInterval(start) * total_time / dt_);
    }
  }
};


//--------------------------------------------------------------------------------------------------
// Initialize the solver
//--------------------------------------------------------------------------------------------------
void Solver::Initialize()
{
  out_.PrintMyLogo();
  out_.PrintLogger(" INITIALIZE SPH SOLVER :\n\n");

  // Initialize material database
  LoadMatDatabase();
  out_.PrintLogger("  - Material Database Initializied\n");

  // Initialize SPH settings
  LoadSettings();
  out_.PrintLogger("  - SPH Settings Initializied\n");

  // Initialize particles data
  LoadParticles();
  out_.PrintLogger("  - Particles Initializied\n");

  out_.PrintTime();
  
  // Set weibull flaws when using grady-kipp model
  WeibullDistribution();

  // initialize grids
  if (neighbor_search_ == NeighborSearchEnum::GRID || neighbor_search_ == NeighborSearchEnum::HYBRID)
  {
    grid_.InitGrid(hsml_ref_ * 2.0, volume_min_, volume_max_, particles_, particles_num_, current_hmax_, kernel_);
  }

  // set average kernel for artificial stress
  avg_kernel_ = artificial_stress_ == ArtificialStressEnum::STANDARD ? kernel_.KernelValue(hsml_ref_/kEta, hsml_ref_) : 0.0;

  current_step_ = 0;
  current_time_ = 0.0;
  dt_ = dt_min_;

  out_.PrintLogger("\n");
}


//--------------------------------------------------------------------------------------------------
// Build all neighbors
//--------------------------------------------------------------------------------------------------
void Solver::BuildNeighbors(bool calc_grav)
{
  switch (neighbor_search_)
  {
  // All pair neighbor search: O(N^2), not recommended
  case(NeighborSearchEnum::ALL_PAIR):
  {
    #pragma omp parallel for schedule(guided)
    for (int p = 0; p < particles_num_; p++)
    {
      if (particles_[p].iflag < 0)
      {
        // particle with iflag < 0 has 0 neighbors
        std::vector<Neighbor>().swap(neighbors_[p]);
        MY_ASSERT(neighbors_[p].size() == 0);
        continue;
      }
      if (neighbors_[p].size() == 0)
      {
        neighbors_[p].reserve(60);
      }

      const Particle& particle = particles_[p];
      int neighbor_num = 0;
      for (int n = 0; n < particles_num_; n++)
      {
        const Particle& neighborParticle = particles_[n];

        if (n == p) // particles with flag == -1 should also be considered as neighbors
        {
          continue;
        }

        Vec3d xij  = particle.pos - neighborParticle.pos;
        double hij = 0.5 * (particle.hsml + neighborParticle.hsml);

        // Judge the distance
        double dist2 = xij.dot(xij);
        if (dist2 < pow2(hij * kernel_.get_kappa()))
        {
          double dist       = sqrt(dist2);
          double kernel     = kernel_.KernelValue(dist, hij);
          Vec3d gradKernel  = kernel_.KernelGradient(dist, hij, xij);
          Vec3d vij         = particle.vel - neighborParticle.vel;

          neighbor_num += 1;
          if (neighbor_num <= neighbors_[p].size())
          {
            neighbors_[p][neighbor_num - 1].ResetNeighbor(n, dist, hij, xij, vij, kernel, gradKernel);
          }
          else
          {
            neighbors_[p].push_back(Neighbor(n, dist, hij, xij, vij, kernel, gradKernel));
          }
        }
      }
      neighbors_[p].resize(neighbor_num);
    }
    if (calc_grav) { UpdateGravity(); }
    break;
  }

  // K-d tree neighbor search O(NlogN)
  case(NeighborSearchEnum::KDTREE):
  {
    // Build new kd tree
    KdTree tree(particles_, kernel_, current_hmax_, hsml_ref_, particles_num_, calc_grav);

    #pragma omp parallel for schedule(guided)
    for (int p = 0; p < particles_num_; p++)
    {
      if (particles_[p].iflag < 0)
      {
        // particle with iflag < 0 has 0 neighbors
        std::vector<Neighbor>().swap(neighbors_[p]);
        MY_ASSERT(neighbors_[p].size() == 0);
        continue;
      }
      if (neighbors_[p].size() == 0)
      {
        neighbors_[p].reserve(60);
      }

      // Neighbor search with k-d tree
      int neighbor_num = 0;
      tree.NeighborSearch(p, &neighbors_[p], &neighbor_num);
      neighbors_[p].resize(neighbor_num);

      // Calculate gravity with k-d tree
      if (calc_grav)
      {
        particles_[p].grav = Vec3d::Zero();
        tree.CalcSelfGravity(p, &particles_[p].grav);
        particles_[p].grav *= kG;
      }
    }

    // Clear kd tree
    tree.Clear();

    break;
  }

  // Grid neighbor search O(N): more efficient for shock-only case
  case(NeighborSearchEnum::GRID):
  {
    // Rebuild the grid
    grid_.Clear(false);
    grid_.setHmax(current_hmax_);
    for (int p = 0; p < particles_num_; p++)
    {
      grid_.Insert(p);
    }

    // Neighbor search with grid
    #pragma omp parallel for schedule(guided)
    for (int p = 0; p < particles_num_; p++)
    {
      if (particles_[p].iflag < 0)
      {
        // particle with iflag < 0 has 0 neighbors
        std::vector<Neighbor>().swap(neighbors_[p]);
        MY_ASSERT(neighbors_[p].size() == 0);
        continue;
      }
      if (neighbors_[p].size() == 0)
      {
        neighbors_[p].reserve(60);
      }

      // Neighbor search with k-d tree
      int neighbor_num = 0;
      grid_.NeighborSearch(p, &neighbors_[p], &neighbor_num);
      neighbors_[p].resize(neighbor_num);
    }
    if (calc_grav) { UpdateGravity(); }
    break;
  }

  // Hybrid neighbor search
  case (NeighborSearchEnum::HYBRID):
  {
    if (calc_grav) // use k-d tree for both neighbor search and gravity calculation
    {
      // Build new kd tree
      KdTree tree(particles_, kernel_, current_hmax_, hsml_ref_, particles_num_, calc_grav);

      #pragma omp parallel for schedule(guided)
      for (int p = 0; p < particles_num_; p++)
      {
        if (particles_[p].iflag < 0)
        {
          // particle with iflag < 0 has 0 neighbors
          std::vector<Neighbor>().swap(neighbors_[p]);
          MY_ASSERT(neighbors_[p].size() == 0);
          continue;
        }
        if (neighbors_[p].size() == 0)
        {
          neighbors_[p].reserve(60);
        }

        // Neighbor search with k-d tree
        int neighbor_num = 0;
        tree.NeighborSearch(p, &neighbors_[p], &neighbor_num);
        neighbors_[p].resize(neighbor_num);

        // Calculate gravity with k-d tree
        particles_[p].grav = Vec3d::Zero();
        tree.CalcSelfGravity(p, &particles_[p].grav);
        particles_[p].grav *= kG;
      }

      // Clear kd tree
      tree.Clear();
    }
    else // use grid for neighbor search only
    {
      // Rebuild the grid
      grid_.Clear(false);
      grid_.setHmax(current_hmax_);
      for (int p = 0; p < particles_num_; p++)
      {
        grid_.Insert(p);
      }

      // Neighbor search with grid
      #pragma omp parallel for schedule(guided)
      for (int p = 0; p < particles_num_; p++)
      {
        if (particles_[p].iflag < 0)
        {
          // particle with iflag < 0 has 0 neighbors
          std::vector<Neighbor>().swap(neighbors_[p]);
          MY_ASSERT(neighbors_[p].size() == 0);
          continue;
        }
        if (neighbors_[p].size() == 0)
        {
          neighbors_[p].reserve(60);
        }

        // Neighbor search with k-d tree
        int neighbor_num = 0;
        grid_.NeighborSearch(p, &neighbors_[p], &neighbor_num);
        neighbors_[p].resize(neighbor_num);
      }
    }
    break;
  }

  default:
    break;
  }

  // Add symmetric ghost neighbors
  if (symmetric_plane_.isInitialized())
  {
    #pragma omp parallel for schedule(guided)
    for (int p = 0; p < particles_num_; p++)
    {
      double dist_to_plane = symmetric_plane_.GetDistance(particles_[p]);
      // Note that only (Ax + By + Cz +D > 0) is mirrored
      if (dist_to_plane > 0.0 && dist_to_plane < 0.5 * (current_hmax_ + particles_[p].hsml))
      {
        int real_neighbor_size = neighbors_[p].size();
        neighbors_[p].reserve(real_neighbor_size + 50);

        Particle ghost_particle; // Only temporarily used
        for (int n = 0; n < real_neighbor_size; n++)
        {
          symmetric_plane_.GetGhostParticle(particles_[neighbors_[p][n].id], ghost_particle);

          Vec3d xij  = particles_[p].pos - ghost_particle.pos;
          double hij = 0.5 * (particles_[p].hsml + ghost_particle.hsml);
          // Judge the distance
          double dist2 = xij.dot(xij);
          if (dist2 < pow2(hij * kernel_.get_kappa()))
          {
            double dist       = sqrt(dist2);
            double kernel     = kernel_.KernelValue(dist, hij);
            Vec3d gradKernel  = kernel_.KernelGradient(dist, hij, xij);
            Vec3d vij         = particles_[p].vel - ghost_particle.vel;
            neighbors_[p].push_back(Neighbor(neighbors_[p][n].id, dist, hij, xij, vij, kernel, gradKernel));
            neighbors_[p].back().iflag = -1; // Set ghost flag
          }
        }

        // Note that the particle itself also has a ghost particle, which is of great importance
        symmetric_plane_.GetGhostParticle(particles_[p], ghost_particle);
        Vec3d xij  = particles_[p].pos - ghost_particle.pos;
        double hij = 0.5 * (particles_[p].hsml + ghost_particle.hsml);
        double dist2 = xij.dot(xij);
        if (dist2 < pow2(hij * kernel_.get_kappa()))
        {
          double dist       = sqrt(dist2);
          double kernel     = kernel_.KernelValue(dist, hij);
          Vec3d gradKernel  = kernel_.KernelGradient(dist, hij, xij);
          Vec3d vij         = particles_[p].vel - ghost_particle.vel;
          neighbors_[p].push_back(Neighbor(p, dist, hij, xij, vij, kernel, gradKernel));
          neighbors_[p].back().iflag = -1; // Set ghost flag
        }
      }

    }
  }
  
}


//--------------------------------------------------------------------------------------------------
// Update the neighbors' info only without re-search
//--------------------------------------------------------------------------------------------------
void Solver::UpdateNeighbors()
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
      neighbors_[p][n].xij = particles_[p].pos - particles_[neighbors_[p][n].id].pos;
      neighbors_[p][n].dist = neighbors_[p][n].xij.norm();
      neighbors_[p][n].hij = 0.5 * (particles_[p].hsml + particles_[neighbors_[p][n].id].hsml);
      neighbors_[p][n].kernel_value = kernel_.KernelValue(neighbors_[p][n].dist, neighbors_[p][n].hij);
      neighbors_[p][n].kernel_gradient = kernel_.KernelGradient(neighbors_[p][n].dist, neighbors_[p][n].hij, neighbors_[p][n].xij);
      neighbors_[p][n].vij = particles_[p].vel - particles_[neighbors_[p][n].id].vel;
    }
    // ghost neighbors not considered here?
  }
};


//--------------------------------------------------------------------------------------------------
// Update time step
//--------------------------------------------------------------------------------------------------
void Solver::UpdateTimeStep()
{
  dt_ = dt_max_;
  double coef1 = 0.2;
  double coef2 = 0.7;

  switch (current_run_mode_)
  {
    case (RunMode::SHOCK):
    {
      #pragma omp parallel for schedule(guided) reduction(min:dt_)
      for (int p = 0; p < particles_num_; p++)
      {
        if (particles_[p].iflag <= 0)
        {
          continue;
        }

        // dt = min{h / c}
        // General case for shock simulation
        double dt1 = coef1 * particles_[p].hsml / (particles_[p].c_sound + particles_[p].vel.norm());
        
        // // Tidal case, we consider sound speed only 
        // double dt1 = coef1 * particles_[p].hsml / particles_[p].c_sound;

        // dt = min{c / dvdt}
        double dvdt = sqrt(particles_[p].vel.dot(particles_[p].vel));
        double dt2 = dvdt > kEps ? (coef2 * particles_[p].c_sound / dvdt) : dt_max_;

        // dt = min{(rho + 0.1 * rho0) / drhodt}
        double dt3 = fabs(particles_[p].drhodt) > kEps ? (coef2 * (particles_[p].rho + 0.1 * particles_[p].rho0) / fabs(particles_[p].drhodt)) : dt_max_;

        // dt = min{(e + 0.1 * c^2) / dedt}
        double dt4 = fabs(particles_[p].dedt) > kEps ? (coef2 * (fabs(particles_[p].energy) + 0.1 * pow2(particles_[p].c_sound)) / fabs(particles_[p].dedt)) : dt_max_;

        // dt = min{(D^(1/3) + 0.1) / dD^(1/3)dt}
        double dt5 = particles_[p].ddamage_root3_dt > kEps ? (coef2 * (particles_[p].damage_root3 + 0.1) / particles_[p].ddamage_root3_dt) : dt_max_;

        // dt = min{alpha / dalpha_dt}
        double dt6 = dt_max_, dt7 = dt_max_;
        Material& mat = mat_database[parts_[particles_[p].part_id].mat_id];
        if ((mat.porosity == PorosityEnum::P_ALPHA) && (particles_[p].alpha > 1.0+kEps))
        {
          dt6 = particles_[p].dalpha_dt < -kEps ? (-coef2 * particles_[p].alpha / particles_[p].dalpha_dt) : dt_max_;
          dt7 = fabs(particles_[p].dpdt) > kEps ? (0.1 * mat.p_alpha.ps / fabs(particles_[p].dpdt)) : dt_max_;
        }

        dt_ = max(dt_min_, min(dt_, dt1, dt3, dt4, dt5, dt6, dt7));
      }
      break;
    }

    case (RunMode::RELAX):
    case (RunMode::LATERUN):
    {
      #pragma omp parallel for schedule(guided) reduction(min:dt_)
      for (int p = 0; p < particles_num_; p++)
      {
        if (particles_[p].iflag <= 0)
        {
          continue;
        }

        // dt = min{h / c}
        double dt1 = coef2 * particles_[p].hsml / (particles_[p].c_sound + particles_[p].vel.norm());

        // dt = min{(rho + 0.1 * rho0) / drhodt}
        double dt3 = fabs(particles_[p].drhodt) > kEps ? (coef2 * (particles_[p].rho + 0.1 * particles_[p].rho0) / fabs(particles_[p].drhodt)) : dt_max_;

        // dt = min{(e + 0.1 * c^2) / dedt}
        double dt4 = fabs(particles_[p].dedt) > kEps ? (coef2 * (fabs(particles_[p].energy) + 0.1 * pow2(particles_[p].c_sound)) / fabs(particles_[p].dedt)) : dt_max_;

        dt_ = max(dt_min_, min(dt_, dt1, dt4));
      }
      break;
    }

    default:
      break;

  }

  MY_ASSERT(std::isfinite(dt_));
}


//--------------------------------------------------------------------------------------------------
// update gravity
//--------------------------------------------------------------------------------------------------
void Solver::UpdateGravity()
{
  switch (calc_gravity_)
  {
  case (CalcGravityEnum::SELF_GRAVITY):
  {
    // // All pair gravity calculate: O(N^2)
    // #pragma omp parallel for schedule(guided)
    // for (int p = 0; p < particles_num_; p++)
    // {
    //   if (particles_[p].iflag < 0)
    //   {
    //     continue;
    //   };
    //   particles_[p].grav = Vec3d::Zero();
    //   for (int n = 0; n < particles_num_; n++)
    //   {
    //     if (n == p || particles_[n].iflag < 0)
    //     {
    //       continue;
    //     };

    //     Vec3d xij  = particles_[p].pos - particles_[n].pos;
    //     double hsml = 0.5 * (particles_[p].hsml + particles_[n].hsml);
    //     double dist = max(xij.norm(), hsml); // add soften gravity
    //     particles_[p].grav -= kG*particles_[n].mass/pow3(dist) * xij;
    //   }
    // }

    // kdTree gravity calculate: O(NlogN)
    KdTree tree(particles_, kernel_, current_hmax_, hsml_ref_, particles_num_, true);
    #pragma omp parallel for schedule(guided)
    for (int p = 0; p < particles_num_; p++)
    {
      if (particles_[p].iflag < 0)
      {
        continue;
      }
      particles_[p].grav = Vec3d::Zero();
      tree.CalcSelfGravity(p, &particles_[p].grav);
      particles_[p].grav *= kG;
    }
    tree.Clear();

    break;
  };
  
  default:
    break;
  }
};
