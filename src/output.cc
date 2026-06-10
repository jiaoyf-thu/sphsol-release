#include"output.h"


//--------------------------------------------------------------------------------------------------
// Initialize the output time list
//--------------------------------------------------------------------------------------------------
void Output::Initialize(double total_time, int total_step, const std::vector<int>& output_flag)
{
  output_step_ = 0;

  output_time_.resize(total_step + 1);

  for (int i = 0; i <= total_step; i++)
  {
    output_time_[i] = total_time * static_cast<double>(i) / static_cast<double>(total_step);
  }
  
  output_flag_.clear();

  for(int flag : output_flag)
  {
    output_flag_.push_back(flag);
  }
}


//--------------------------------------------------------------------------------------------------
// Add output time for late stage
//--------------------------------------------------------------------------------------------------
void Output::AddOutputTime(const std::vector<double>& add_time)
{
  for (double time : add_time)
  {
    output_time_.push_back(time);
  }
};


//--------------------------------------------------------------------------------------------------
// Print output times to file
//--------------------------------------------------------------------------------------------------
void Output::PrintTime()
{
  FILE* fp = fopen(kPathTime.c_str(), "w");

  fprintf(fp, "OutputTime\n");
  
  for (int i = 0; i < output_time_.size(); i++)
  {
    fprintf(fp, "%.10e\n", output_time_[i]);
  }

  fclose(fp);
}


//--------------------------------------------------------------------------------------------------
// Print particles data
//--------------------------------------------------------------------------------------------------
void Output::PrintParticles(const Particle* particles, int particles_num)
{
  // print the particles data
  std::string stepNum = std::to_string(output_step_);
  std::string outputPath = kPathOutput.data() + std::string(4 - stepNum.length(), '0') + stepNum + ".csv";

  FILE* fp = fopen(outputPath.c_str(), "w");

  fprintf(fp, "ID");
  for (int i = 0; i < output_flag_.size(); i++)
  {
    if (output_flag_[i] == 1)
    {
      switch (i)
      {
      case(0):
      {
        fprintf(fp, ",X,Y,Z");
        break;
      }
      case(1):
      {
        fprintf(fp, ",VX,VY,VZ");
        break;
      }
      case(2):
      {
        fprintf(fp, ",AX,AY,AZ");
        break;
      }
      case(3):
      {
        fprintf(fp, ",P");
        break;
      }
      case(4):
      {
        fprintf(fp, ",Pmax");
        break;
      }
      case(5):
      {
        fprintf(fp, ",S11,S12,S13,S22,S23,S33");
        break;
      }
      case(6):
      {
        fprintf(fp, ",IFLAG");
        break;
      }
      case(7):
      {
        fprintf(fp, ",MASS");
        break;
      }
      case(8):
      {
        fprintf(fp, ",RHO,RHO_0");
        break;
      }
      case(9):
      {
        fprintf(fp, ",E");
        break;
      }
      case(10):
      {
        fprintf(fp, ",ETmax");
        break;
      }
      case(11):
      {
        fprintf(fp, ",H");
        break;
      }
      case(12):
      {
        fprintf(fp, ",CS");
        break;
      }
      case(13):
      {
        fprintf(fp, ",D");
        break;
      }
      case(14):
      {
        fprintf(fp, ",PART");
        break;
      }
      case(15):
      {
        fprintf(fp, ",LM");
        break;
      }
      case(16):
      {
        fprintf(fp, ",ALPHA");
        break;
      }
      case(17):
      {
        fprintf(fp, ",R12,R13,R23");
        break;
      }
      case(18):
      {
        fprintf(fp, ",GX,GY,GZ");
        break;
      }
      default:
        break;
      }
    }
  }
  fprintf(fp, "\n");


  for (int p = 0; p < particles_num; p++)
  {
    const Particle& pp = particles[p];

    fprintf(fp, "%d", p);
    for (int i = 0; i < output_flag_.size(); i++)
    {
      if (output_flag_[i] == 1)
      {
        switch (i)
        {
        case(0):
        {
          fprintf(fp, ",%.10e,%.10e,%.10e", pp.pos(0), pp.pos(1), pp.pos(2));
          break;
        }
        case(1):
        {
          fprintf(fp, ",%.10e,%.10e,%.10e", pp.vel(0), pp.vel(1), pp.vel(2));
          break;
        }
        case(2):
        {
          fprintf(fp, ",%.10e,%.10e,%.10e", pp.dvdt(0), pp.dvdt(1), pp.dvdt(2));
          break;
        }
        case(3):
        {
          fprintf(fp, ",%.10e", pp.pressure);
          break;
        }
        case(4):
        {
          fprintf(fp, ",%.10e", pp.pressure_peak);
          break;
        }
        case(5):
        {
          fprintf(fp, ",%.10e,%.10e,%.10e,%.10e,%.10e,%.10e", pp.stress(0,0), pp.stress(0,1), pp.stress(0,2),
              pp.stress(1,1), pp.stress(1,2), pp.stress(2,2));
          break;
        }
        case(6):
        {
          fprintf(fp, ",%d", pp.iflag);
          break;
        }
        case(7):
        {
          fprintf(fp, ",%.10e", pp.mass);
          break;
        }
        case(8):
        {
          fprintf(fp, ",%.10e,%.10e", pp.rho, pp.rho0);
          break;
        }
        case(9):
        {
          fprintf(fp, ",%.10e", pp.energy);
          break;
        }
        case(10):
        {
          fprintf(fp, ",%.10e", pp.energy_thermal_peak);
          break;
        }
        case(11):
        {
          fprintf(fp, ",%.10e", pp.hsml);
          break;
        }
        case(12):
        {
          fprintf(fp, ",%.10e", pp.c_sound);
          break;
        }
        case(13):
        {
          fprintf(fp, ",%.10e", pp.damage);
          break;
        }
        case(14):
        {
          fprintf(fp, ",%d", pp.part_id);
          break;
        }
        case(15):
        {
          fprintf(fp, ",%.10e", pp.frag.Lm);
        break;
        }
        case(16):
        {
          fprintf(fp, ",%.10e", pp.alpha);
          break;
        }
        case(17):
        {
          fprintf(fp, ",%.10e,%.10e,%.10e", pp.rotate_rate(0,1), pp.rotate_rate(0,2), pp.rotate_rate(1,2));
          break;
        }
        case(18):
        {
          fprintf(fp, ",%.10e,%.10e,%.10e", pp.grav(0), pp.grav(1), pp.grav(2));
          break;
        }
        default:
          break;
        }
      }
    }
    fprintf(fp, "\n");
  }

  fclose(fp);

  output_step_ += 1;
}


//--------------------------------------------------------------------------------------------------
// Print fragments >= 0 only
//--------------------------------------------------------------------------------------------------
void Output::PrintFragments(const Particle* particles, int particles_num, std::vector<int>* fragments)
{
  std::string stepNum = std::to_string(output_step_);
  std::string outputPath = kPathOutput.data() + std::string("frag") + std::string(4 - stepNum.length(), '0') + stepNum + ".csv";

  FILE* fp = fopen(outputPath.c_str(), "w");

  fprintf(fp, "ID,X,Y,Z,VX,VY,VZ,FRAG\n");

  for (int p = 0; p < particles_num; p++)
  {
    if ((*fragments)[p] < 0)
    {
      continue;
    }
    
    const Particle& pp = particles[p];

    fprintf(fp, "%d", p);
    fprintf(fp, ",%.10e,%.10e,%.10e", pp.pos(0), pp.pos(1), pp.pos(2));
    fprintf(fp, ",%.10e,%.10e,%.10e", pp.vel(0), pp.vel(1), pp.vel(2));
    fprintf(fp, ",%d\n", (*fragments)[p]);
  }

  fclose(fp);
}


//--------------------------------------------------------------------------------------------------
// Print neighbors id, test only
//--------------------------------------------------------------------------------------------------
void Output::PrintNeighbors(std::vector<Neighbor>* neighbors, int particles_num)
{
  FILE* fp = fopen("../Output/test_neighbor_grid.csv", "w");
  for (int p = 0; p < particles_num; p++)
  {
    std::sort(neighbors[p].begin(), neighbors[p].end(), [&](Neighbor lhs, Neighbor rhs)
    {
      return lhs.id < rhs.id;
    });
    for (auto neighbor : neighbors[p])
    {
      fprintf(fp, "%d,", neighbor.id);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  char cc = getchar();
}


//--------------------------------------------------------------------------------------------------
// Print information to terminal and logger
//--------------------------------------------------------------------------------------------------
void Output::PrintLogger(const char* format, ...)
{
#ifdef MY_OPENMP
#pragma omp critical
#endif
  {
    FILE* logger_;
    logger_ = fopen(kPathLogger.c_str(), "a");

    va_list vaList;
    
    va_start(vaList, format);
    vprintf(format, vaList);
    va_end(vaList);

    va_start(vaList, format);
    vfprintf(logger_, format, vaList);
    va_end(vaList);

    fclose(logger_);
  }
}


//--------------------------------------------------------------------------------------------------
// Print logo
//--------------------------------------------------------------------------------------------------
void Output::PrintMyLogo()
{
  Vec3i current_date = GetCurrentDate();

  PrintLogger("\n------------------------------  Current Date: %4d-%02d-%02d  ------------------------------\n",
      current_date(0), current_date(1), current_date(2));

  PrintLogger("      _______  __ __________  __ \n");
  PrintLogger("     / __/ _ \\/ // / __/ __ \\/ /\n");
  PrintLogger("    _\\ \\/ ___/ _  /\\ \\/ /_/ / /__\n");
  PrintLogger("   /___/_/  /_//_/___/\\____/____/   %s\n\n", kVersion.c_str());
  PrintLogger("----------------------------------------------------------------------------------------\n\n");

};

