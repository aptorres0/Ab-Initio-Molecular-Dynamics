#ifndef _CONSTANTS_H
#define _CONSTANTS_H

#include <cmath>

namespace Ar {
/**
 * define characteristic quantities for the system
 */

// Characteristic length/minimum interacting distance for Ar atoms:
const double gLength = 0.34; // [nm]
const double gVolume = gLength * gLength * gLength;

// Temperature of the system:
const double gTemperatureArgon = 120.0; // [K]

// Boltzmann constant:
const double gBoltzmannConstant = 13.80649; // [yJ/K] or 10^-24 [J/K]

// Characteristic energy:
const double gEnergy =
    gBoltzmannConstant * gTemperatureArgon; // [yJ] or 10^-24 [J]

// Mass of an Argon atom:
const double gMass = 39.9; // [u] (atomic units 1[u]=1.661e-27[kg])

// Number density of Argon atoms (i.e. particles per units volume):
const double gNumberDensity_r = 0.8;                        // [1]
const double gNumberDensity = gNumberDensity_r / (gVolume); // [nm^-3]

// NOTE: need to figure out what units the below 3 are in

// Characteristic time:
const double gTime = sqrt(gMass * gLength * gLength / gEnergy); //

// Characteristic velocity:
const double gVelocity = gLength / gTime;

// Characteristic acceleration:
const double gAcceleration = gLength / (gTime * gTime);

// Characteristic momentum:
const double gMomentum = sqrt(gMass * gEnergy);

} // namespace Ar

#endif