#ifndef SPHSOL_PHYSICS_H_
#define SPHSOL_PHYSICS_H_

#include"math.h"
#include"particle.h"


//--------------------------------------------------------------------------------------------------
// Equation of state
//--------------------------------------------------------------------------------------------------
void EquatonOfState(Particle* particle, const Material& mat);
double EquatonOfStateCold(const Particle* particle, const Material& mat);


//--------------------------------------------------------------------------------------------------
// Strength model
//--------------------------------------------------------------------------------------------------
double StrengthModel(const Material& mat, double damage, double J2, double pressure = 0.0);


//--------------------------------------------------------------------------------------------------
// Shear damage
//--------------------------------------------------------------------------------------------------
double ShearDamage(double plastic_strain, double pressure, const Material& mat);


#endif
