#ifndef SPHSOL_MATERIAL_H_
#define SPHSOL_MATERIAL_H_

#include"global.h"


enum class EosEnum
{
  TILLOTSON,
  SIMPLIFIED_TILLOTSON,
  GRUNEISEN,
  ELASTIC
};
enum class PorosityEnum
{
  NONE,
  P_ALPHA
};
enum class YieldEnum
{
  VON_MISES,
  LUND
};
enum class DamageEnum
{
  NONE,
  GRADY_KIPP
};


struct TillotsonPara
{
  double a, b, alpha, beta;
  double A, B;
  double E0, Eiv, Ecv;
};
struct SimplifiedTillotsonPara
{
  double c, A;
};
struct GruneisenPara
{
  double s;
  double c0, gama0;
};
struct PressureAlphaPara
{
  double alpha0,rho_solid;
  double ps,pe;
};
struct VonMisesPara
{
  double Y0;
};
struct LundPara
{
  double Yi0, Yd0, Ym;
  double mui, mud;
};
struct GradyKippPara
{
  double m, k;
  double cg_ce;
  double p_bd,p_bp;
};


//--------------------------------------------------------------------------------------------------
// Material data
//--------------------------------------------------------------------------------------------------
struct Material
{
  std::string  mat_name;

  // Commonly used modulus
  double shear_modulus, bulk_modulus, elastic_modulus;

  // Model selection
  EosEnum       eos;
  PorosityEnum  porosity;
  YieldEnum     yield;
  DamageEnum    damage;

  // Model parameter
  TillotsonPara           tillotson;
  SimplifiedTillotsonPara sim_tillotson;
  GruneisenPara           gruneisen;
  PressureAlphaPara       p_alpha;
  VonMisesPara            von_mises;
  LundPara                lund;
  GradyKippPara           grady_kipp;
};


// Material database
extern std::vector<Material> mat_database;


// // Initialize the material database
// void LoadMatDatabase();


#endif
