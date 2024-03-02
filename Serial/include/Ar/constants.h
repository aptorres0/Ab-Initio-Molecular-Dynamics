#ifndef _CONSTANTS_H
#define _CONSTANTS_H

#include <cmath>

namespace Ar {
/**
 * define characteristic quantities for the system
 */

// Conversion factors:
const double gElectronVolt2Joule = 1.602176634e-19;  // [J/eV]
const double gAtomicMassUnit2Kg = 1.66053906660e-27; // [kg/u]
const double gNanoMeter2Meter = 1.0e-9;              // [m/nm]

// Characteristic length/minimum interacting distance for Ar atoms:
const double gLength = 0.34;                        // [nm]
const double gLength2 = gLength * gLength;          // [nm^2]
const double gVolume = gLength * gLength * gLength; // [nm^3]

// Temperature of the system:
const double gTemperature = 120.0; // [K]

// Boltzmann constant:
const double gBoltzmannConstant = 8.617333262e-5; // [eV/K]

// Characteristic energy:
const double gEnergy = gBoltzmannConstant * gTemperature; // [eV]

// Mass of an Argon atom:
const double gMass = 39.9; // [u] (atomic units 1[u]=1.661e-27[kg])

// Number density of Argon atoms (i.e. particles per units volume):
const double gNumberDensity_r = 0.8;                        // [1]
const double gNumberDensity = gNumberDensity_r / (gVolume); // [nm^-3]

// Characteristic time:
const double gTime =
    sqrt(gMass * gLength * gLength / gEnergy) * 0.101805; // [ps]

// Characteristic velocity:
const double gVelocity = gLength / gTime;

// Characteristic acceleration:
const double gAcceleration = gLength / (gTime * gTime);

// Characteristic momentum:
const double gMomentum = sqrt(gMass * gEnergy);

// Pressure
const double gPressure = (gEnergy / gVolume) * (160.218); // [MPa]

// Specific heat capacity
// const double gSpecificHeatCapacity = 3.0 * gBoltzmannConstant; // [eV/K]
const double gSpecificHeatCapacity =
    (gBoltzmannConstant / gMass) * 893383; // [J/kg/K]

} // namespace Ar

#endif