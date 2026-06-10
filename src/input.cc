#include"solver.h"


//--------------------------------------------------------------------------------------------------
// Load settings from file
//--------------------------------------------------------------------------------------------------
void Solver::LoadSettings()
{
  FILE* fp = fopen(kPathSetting.c_str(), "r");

  fscanf(fp, "%*[^\n]\n");
  fscanf(fp, "%*[^\n]\n");
  fscanf(fp, "%*[^\n]\n");
  fscanf(fp, "%*[^\n]\n");
  fscanf(fp, "%*[^\n]\n");
  fscanf(fp, "%*[^\n]\n");

  // Load solver settings
  char value[128];
  fscanf(fp, "%*30c%s\n", &value);
  if (std::strcmp(value, "LEAP_FROG") == 0) { integrator_ = IntegratorEnum::LEAP_FROG; };
  if (std::strcmp(value, "PREDICT_CORRECT") == 0) { integrator_ = IntegratorEnum::PREDICT_CORRECT; };

  fscanf(fp, "%*30c%s\n", &value);
  if (std::strcmp(value, "STANDARD") == 0) { artificial_viscosity_ = ArtificialViscosityEnum::STANDARD; };
  if (std::strcmp(value, "NONE") == 0) { artificial_viscosity_ = ArtificialViscosityEnum::NONE; };

  fscanf(fp, "%*30c%s\n", &value);
  if (std::strcmp(value, "STANDARD") == 0) { artificial_heat_ = ArtificialHeatEnum::STANDARD; };
  if (std::strcmp(value, "NONE") == 0) { artificial_heat_ = ArtificialHeatEnum::NONE; };

  fscanf(fp, "%*30c%s\n", &value);
  if (std::strcmp(value, "STANDARD") == 0) { artificial_stress_ = ArtificialStressEnum::STANDARD; };
  if (std::strcmp(value, "NONE") == 0) { artificial_stress_ = ArtificialStressEnum::NONE; };

  fscanf(fp, "%*30c%s\n", &value);
  if (std::strcmp(value, "STANDARD") == 0) { xsph_ = XsphEnum::STANDARD; };
  if (std::strcmp(value, "NONE") == 0) { xsph_ = XsphEnum::NONE; };

  fscanf(fp, "%*30c%s\n", &value);
  if (std::strcmp(value, "CONTINUOUS_DENSITY") == 0) { density_update_ = DensityEnum::CONTINUOUS_DENSITY; };
  if (std::strcmp(value, "DENSITY_SUMMATION") == 0) { density_update_ = DensityEnum::DENSITY_SUMMATION; };

  fscanf(fp, "%*30c%s\n", &value);
  if (std::strcmp(value, "ALL_PAIR") == 0) { neighbor_search_ = NeighborSearchEnum::ALL_PAIR; };
  if (std::strcmp(value, "GRID") == 0) { neighbor_search_ = NeighborSearchEnum::GRID; };
  if (std::strcmp(value, "KDTREE") == 0) { neighbor_search_ = NeighborSearchEnum::KDTREE; };
  if (std::strcmp(value, "HYBRID") == 0) { neighbor_search_ = NeighborSearchEnum::HYBRID; };

  // gravity
  fscanf(fp, "%*30c%s\n", &value);

  // gravity_ << 0.0, 0.0, 0.0;
  // if (std::strcmp(value, "CONST") == 0) 
  // {
  //   calc_gravity_ = CalcGravityEnum::CONST;
  //   fscanf(fp, "%*30c%lf,%lf,%lf\n", &gravity_(0), &gravity_(1), &gravity_(2));
  //   fscanf(fp, "%*[^\n]\n");
  // };
  // if (std::strcmp(value, "SPHERE") == 0) 
  // {
  //   calc_gravity_ = CalcGravityEnum::SPHERE;
  //   fscanf(fp, "%*[^\n]\n");
  //   fscanf(fp, "%*30c%lf,%lf,%lf,%lf,%lf\n", &sphere_gravity_.center(0), &sphere_gravity_.center(1), &sphere_gravity_.center(2), &sphere_gravity_.radius, &sphere_gravity_.rho);
  //   // sphere_gravity_.escape_vel = sqrt(2.0*kG * 4.0/3.0*kPI*pow3(sphere_gravity_.radius) * sphere_gravity_.rho / sphere_gravity_.radius);
  // };
  if (std::strcmp(value, "SELF_GRAVITY") == 0) 
  {
    // we will add the const and sphere gravity to the self-gravity
    // keep them as zero if not used !!!
    calc_gravity_ = CalcGravityEnum::SELF_GRAVITY;
    // fscanf(fp, "%*30c%lf,%lf,%lf\n", &gravity_(0), &gravity_(1), &gravity_(2));
    // fscanf(fp, "%*30c%lf,%lf,%lf,%lf,%lf\n", &sphere_gravity_.center(0), &sphere_gravity_.center(1), &sphere_gravity_.center(2), &sphere_gravity_.radius, &sphere_gravity_.rho);
  };
  if (std::strcmp(value, "NONE") == 0) 
  {
    calc_gravity_ = CalcGravityEnum::NONE;
    // fscanf(fp, "%*[^\n]\n");
    // fscanf(fp, "%*[^\n]\n");
  };

  // Add force
  int tmp;
  fscanf(fp, "%*30c%d\n", &tmp);
  additional_force_ = (tmp == 1 ? true : false);

  // rotation
  fscanf(fp, "%*30c%s\n", &value);
  if (std::strcmp(value, "TARGET") == 0)
  {
    is_rotating_ = true;
    rotation_.part_id = 0; // defalut the target
    fscanf(fp, "%*30c%lf,%lf,%lf,%lf,%lf,%lf\n", &rotation_.pos0(0), &rotation_.pos0(1), &rotation_.pos0(2), 
              &rotation_.omg0(0), &rotation_.omg0(1), &rotation_.omg0(2));
  };
  if (std::strcmp(value, "NONE") == 0)
  {
    is_rotating_ = false;
    fscanf(fp, "%*[^\n]\n");
  };

  fscanf(fp, "%*30c%lf\n", &damp_time_);

  // symmetric plane
  fscanf(fp, "%*30c%s\n", &value);
  if (std::strcmp(value, "PLANE") == 0) 
  {
    double coef[4];
    fscanf(fp, "%*30c%lf,%lf,%lf,%lf\n", &coef[1], &coef[1], &coef[2], &coef[3]);
    symmetric_plane_.Init(coef);
  }
  else
  {
    fscanf(fp, "%*[^\n]\n");
  };

  // relax
  fscanf(fp, "%*30c%d\n", &tmp);
  if (tmp == 1 || tmp == 2) // 1 for zero pressure; 2 for given pressure
  {
    run_modes_.push_back(RunMode::RELAX);
    use_sphere_pressure_before_relax_ = (tmp == 1 ? true : false);
  }
  
  // run_modes_.push_back(RunMode::SHOCK); // default mode

  // shock stage
  fscanf(fp, "%*30c%d\n", &tmp);
  if (tmp == 1)
  {
    run_modes_.push_back(RunMode::SHOCK);
  }

  // late stage
  fscanf(fp, "%*30c%d\n", &tmp);
  if (tmp == 1)
  {
    run_modes_.push_back(RunMode::LATERUN);
  }

  fscanf(fp, "%*[^\n]\n");
  fscanf(fp, "%*[^\n]\n");
  fscanf(fp, "%*[^\n]\n");

  // Load limits
  fscanf(fp, "%*30c%d\n", &particles_num_);
  particles_ = new Particle[particles_num_];
  neighbors_ = new std::vector<Neighbor>[particles_num_];

  fscanf(fp, "%*30c%d\n", &parts_num_);
  parts_ = new Part[parts_num_];

  for (int part = 0; part < parts_num_; part++)
  {
    switch (part)
    {
    case(0):
    {
      fscanf(fp, "%*30c%d", &parts_[part].mat_id);
      break;
    }
    default:
    {
      fscanf(fp, ",%d", &parts_[part].mat_id);
      break;
    }
    }
  }
  fscanf(fp, "\n");

  // run time
  double dt_ref;
  if (tmp == 0) // tmp represents the late_stage
  {
    fscanf(fp, "%*30c%lf\n", &total_time_);
  }
  else
  {
    fscanf(fp, "%*30c%lf,%lf\n", &total_time_, &late_time_);
  }
  fscanf(fp, "%*30c%lf\n", &dt_ref);
  fscanf(fp, "%*30c%lf,%lf\n", &dt_min_, &dt_max_);
  dt_max_ *= dt_ref;
  dt_min_ *= dt_ref;
  dt_ref_  = dt_ref;
  dt_stp_  = 0.0;

  // hsml
  fscanf(fp, "%*30c%lf\n", &hsml_ref_);
  fscanf(fp, "%*30c%lf,%lf\n", &hsml_min_, &hsml_max_);
  hsml_max_ *= hsml_ref_;
  hsml_min_ *= hsml_ref_;

  fscanf(fp, "%*30c%lf,%lf,%lf\n", &volume_min_(0), &volume_min_(1), &volume_min_(2));
  fscanf(fp, "%*30c%lf,%lf,%lf\n", &volume_max_(0), &volume_max_(1), &volume_max_(2));

  fscanf(fp, "%*[^\n]\n");
  fscanf(fp, "%*[^\n]\n");
  fscanf(fp, "%*[^\n]\n");

  // Load input settings
  input_flag_.resize(13);
  for (int i = 0; i < 13; i++)
  {
    fscanf(fp, "%*30c%d\n", &input_flag_[i]);
  }

  fscanf(fp, "%*[^\n]\n");
  fscanf(fp, "%*[^\n]\n");
  fscanf(fp, "%*[^\n]\n");

  // Load output settings
  int output_step1, output_step2;
  if (tmp == 0)
  {
    fscanf(fp, "%*30c%d\n", &output_step1);
  }
  else
  {
    fscanf(fp, "%*30c%d,%d\n", &output_step1, &output_step2);
  }

  std::vector<int> output_flag(19);
  for (int i = 0; i < 19; i++)
  {
    fscanf(fp, "%*30c%d\n", &output_flag[i]);
  }

  out_.Initialize(total_time_, output_step1, output_flag);

  if (tmp == 1)
  {
    std::vector<double> add_output_time;
    add_output_time.resize(output_step2);
    for (int i = 1; i <= output_step2; i++)
    {
      add_output_time[i-1] = total_time_ + (late_time_-total_time_) * static_cast<double>(i) / static_cast<double>(output_step2);
    }
    out_.AddOutputTime(add_output_time);
  }

  fclose(fp);
}


//--------------------------------------------------------------------------------------------------
// Load particles from file
//--------------------------------------------------------------------------------------------------
void Solver::LoadParticles()
{
  FILE* fp = fopen(kPathInput.c_str(), "r");

  char title[128];

  fscanf(fp, "%s\n", title);

  int id = 0;

  double oneThird  = 1.0 / 3.0;
  double fourThird = 4.0 / 3.0;

  Particle pp; // temp particle  

  while (!feof(fp))
  {
    fscanf(fp, "%d,%lf,%lf,%lf,%lf,%lf,%lf,%d,%lf,%lf", &id, &pp.pos(0), &pp.pos(1), &pp.pos(2),
      &pp.vel(0), &pp.vel(1), &pp.vel(2), &pp.part_id, &pp.mass, &pp.rho);

    // Check input_flag
    if (input_flag_[6] == 1)
    {
      fscanf(fp, ",%lf", &pp.rho0);
    }
    if (input_flag_[7] == 1)
    {
      fscanf(fp, ",%d", &pp.iflag);
    }
    if (input_flag_[8] == 1)
    {
      fscanf(fp, ",%lf", &pp.damage);
    }
    if (input_flag_[9] == 1)
    {
      fscanf(fp, ",%lf", &pp.pressure);
    }
    if (input_flag_[10] == 1)
    {
      fscanf(fp, ",%lf", &pp.alpha);
    }
    if (input_flag_[11] == 1)
    {
      fscanf(fp, ",%lf", &pp.energy);
    }
    if (input_flag_[12] == 1)
    {
      fscanf(fp, ",%lf", &pp.hsml);
    }
    fscanf(fp, "\n");

    // Input data
    particles_[id].pos            = pp.pos;
    particles_[id].vel            = pp.vel;
    particles_[id].part_id        = pp.part_id;
    particles_[id].mass           = pp.mass;
    particles_[id].rho            = pp.rho;
    particles_[id].rho0           = input_flag_[6]  == 1 ? pp.rho0 : pp.rho;
    particles_[id].iflag          = input_flag_[7]  == 1 ? pp.iflag : 1;
    particles_[id].damage         = input_flag_[8]  == 1 ? pp.damage : 0.0;
    particles_[id].damage_root3   = input_flag_[8]  == 1 ? pow(pp.damage,oneThird) : 0.0;
    particles_[id].pressure       = input_flag_[9]  == 1 ? pp.pressure : 0.0;
    // input_flag_[10] is for porosity alpha only
    particles_[id].energy         = input_flag_[11] == 1 ? pp.energy : 0.0;
    particles_[id].hsml           = input_flag_[12] == 1 ? pp.hsml : kEta * pow(particles_[id].mass / particles_[id].rho, oneThird);

    // Auto initialized data
    Material& mat                    = mat_database[parts_[particles_[id].part_id].mat_id];

    parts_[particles_[id].part_id].particle_id.push_back(id);

    // particles_[id].rho0              = particles_[id].rho;
    particles_[id].dxdt              = Vec3d::Zero();
    particles_[id].dev_stress        = Mat3d::Zero();
    particles_[id].ddev_stress_dt    = Mat3d::Zero();
    particles_[id].stress            = - particles_[id].pressure * Mat3d::Identity();
    particles_[id].tension           = Mat3d::Zero();
    particles_[id].strain_rate       = Mat3d::Zero();
    particles_[id].rotate_rate       = Mat3d::Zero();
    particles_[id].correction        = Mat3d::Identity();
    particles_[id].grav              = Vec3d::Zero();
    particles_[id].pressure_peak     = 0.0;
    particles_[id].energy_cold       = 0.0;
    particles_[id].energy_thermal_peak = 0.0;
    particles_[id].c_sound           = sqrt((mat.bulk_modulus) / particles_[id].rho);
    particles_[id].drhodt            = 0.0;
    particles_[id].dedt              = 0.0;
    particles_[id].ddamage_root3_dt  = 0.0;
    particles_[id].plastic_strain    = 0.0;
    particles_[id].Hi                = 0.0;
    particles_[id].flaws.ntot        = 0;
    particles_[id].flaws.strain_min  = 0.0;
    particles_[id].flaws.strain_max  = 0.0;
    particles_[id].flaws.m0          = 0.0;
    particles_[id].frag.Lm           = 0.0;
    particles_[id].frag.area         = 0.0;
    particles_[id].frag.dAdt         = 0.0;
    particles_[id].flaws.damage_root3_max = 0.0;

    switch (mat.porosity)
    {
    case(PorosityEnum::P_ALPHA):
    {
      particles_[id].alpha     = input_flag_[10] == 1 ? pp.alpha : mat.p_alpha.rho_solid / particles_[id].rho;
      // particles_[id].alpha     = mat.p_alpha.rho_solid / particles_[id].rho;
      particles_[id].alpha_inv = 1.0 / particles_[id].alpha;
      particles_[id].rho0      = mat.p_alpha.rho_solid / mat.p_alpha.alpha0; // in case we are re-starting from a previous run
      break;
    }
    default:
    {
      particles_[id].alpha     = 1.0;
      particles_[id].alpha_inv = 1.0;
      break;
    }
    }
    particles_[id].dalpha_dt   = 0.0;
    particles_[id].dpdt        = 0.0;
    particles_[id].dalpha_drho = 0.0;

    current_hmax_ = max(current_hmax_, particles_[id].hsml);
  }

  fclose(fp);
}


// Material database
std::vector<Material> mat_database;


//--------------------------------------------------------------------------------------------------
// Initialize material database
//--------------------------------------------------------------------------------------------------
void Solver::LoadMatDatabase()
{
  FILE* fp = fopen(kPathMaterials.c_str(), "r");

  fscanf(fp, "%*[^\n]\n");
  fscanf(fp, "%*[^\n]\n");
  fscanf(fp, "%*[^\n]\n");
  fscanf(fp, "%*[^\n]\n");
  fscanf(fp, "%*[^\n]\n");
  fscanf(fp, "%*[^\n]\n");

  // Load materials number
  int materials_num;
  fscanf(fp, "%*30c%d\n", &materials_num);
  mat_database.resize(materials_num);

  char title[128];
  char value[128];
  int mat_id = -1;
  // Load all materials
  for (int i = 0; i < materials_num; i++)  
  {
    fscanf(fp, "%*[^\n]\n");
    fscanf(fp, "%*[^\n]\n");
    fscanf(fp, "%*[^\n]\n");

    while (!feof(fp))
    {
      fscanf(fp, "%s\n", title);
      if (std::strcmp(title, "@END") == 0) { break; };

      if (std::strcmp(title, "@MATERIAL") == 0)
      {
        // General parameters
        fscanf(fp, "%*30c%d\n", &mat_id);
        fscanf(fp, "%*30c%s\n", &value);
        mat_database[mat_id].mat_name = value;

        fscanf(fp, "%*30c%s\n", &value);
        if (std::strcmp(value, "TILLOTSON") == 0) { mat_database[mat_id].eos = EosEnum::TILLOTSON; };
        if (std::strcmp(value, "SIMPLIFIED_TILLOTSON") == 0) { mat_database[mat_id].eos = EosEnum::SIMPLIFIED_TILLOTSON; };
        if (std::strcmp(value, "GRUNEISEN") == 0) { mat_database[mat_id].eos = EosEnum::GRUNEISEN; };

        fscanf(fp, "%*30c%s\n", &value);
        if (std::strcmp(value, "NONE") == 0) { mat_database[mat_id].porosity = PorosityEnum::NONE; };
        if (std::strcmp(value, "P_ALPHA") == 0) { mat_database[mat_id].porosity = PorosityEnum::P_ALPHA; };

        fscanf(fp, "%*30c%s\n", &value);
        if (std::strcmp(value, "VON_MISES") == 0) { mat_database[mat_id].yield = YieldEnum::VON_MISES; };
        if (std::strcmp(value, "LUND") == 0) { mat_database[mat_id].yield = YieldEnum::LUND; };

        fscanf(fp, "%*30c%s\n", &value);
        if (std::strcmp(value, "NONE") == 0) { mat_database[mat_id].damage = DamageEnum::NONE; };
        if (std::strcmp(value, "GRADY_KIPP") == 0) { mat_database[mat_id].damage = DamageEnum::GRADY_KIPP; };

        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].shear_modulus);
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].bulk_modulus);
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].elastic_modulus);
      };

      if (std::strcmp(title, "@TILLOTSON") == 0)
      {
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].tillotson.a);
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].tillotson.b);
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].tillotson.alpha);
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].tillotson.beta);
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].tillotson.A);
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].tillotson.B);
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].tillotson.E0);
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].tillotson.Ecv);
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].tillotson.Eiv);

        // Only used for p-alpha update
        mat_database[mat_id].sim_tillotson.A = mat_database[mat_id].tillotson.A;
        mat_database[mat_id].sim_tillotson.c = mat_database[mat_id].tillotson.a + mat_database[mat_id].tillotson.b;
      };

      if (std::strcmp(title, "@SIMPLIFIED_TILLOTSON") == 0)
      {
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].sim_tillotson.c);
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].sim_tillotson.A);
      };

      if (std::strcmp(title, "@GRUNEISEN") == 0)
      {
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].gruneisen.s);
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].gruneisen.c0);
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].gruneisen.gama0);
      };

      if (std::strcmp(title, "@P_ALPHA") == 0)
      {
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].p_alpha.alpha0);
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].p_alpha.rho_solid);
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].p_alpha.pe);
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].p_alpha.ps);
      };

      if (std::strcmp(title, "@VON_MISES") == 0)
      {
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].von_mises.Y0);
      };

      if (std::strcmp(title, "@LUND") == 0)
      {
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].lund.Yi0);
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].lund.Yd0);
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].lund.Ym);
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].lund.mui);
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].lund.mud);
      };

      if (std::strcmp(title, "@GRADY_KIPP") == 0)
      {
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].grady_kipp.m);
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].grady_kipp.k);
        fscanf(fp, "%*30c%lf\n", &mat_database[mat_id].grady_kipp.cg_ce);
        mat_database[mat_id].grady_kipp.p_bd = 1.23e9; // TODO
        mat_database[mat_id].grady_kipp.p_bp = 2.35e9;
      };
    }

  }
  fclose(fp);
}

