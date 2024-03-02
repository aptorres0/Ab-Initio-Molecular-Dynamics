/**
 * Sim.h
 *
 * Author: Alexander Paul Torres
 * Date: 14 APR 23
 * Class: PHYS 7411 - Computational Physics I
 *
 * Description:
 * Header file with declaration of the Sim class used to simulate the time
 * evolution of a system of 108 Ar atoms in the Lennard-Jones potential and
 * compute thermodinamic properties. Implements the Velocity-Verlet algorithm to
 * integrate the equations of motion. See Sim.cxx for implementation.
 *
 */
#ifndef _SIM_H
#define _SIM_H

#include "Ar/Atom.h"
#include "MyRNG.h" // Random number generator used by thermostats
#include "Vector3.h"

typedef struct RadialDistributionFunction {
  double bin[500];
  int nBins = 500;
  double binWidth = 0.01;
} Histo;

class Sim {
private:
  Ar::Atom *fAtomHead; // Pointer to head of linked list of atoms
  // Sim parameters
  int fNumAtoms = 108;          // number of atoms
  double h = 0.01;              // time step [r.u.]
  double L = cbrt(108.0 / 0.8); // length of simulation cell [r.u.]
  double CutoffRadiusSq = 6.25; // cutoff distance (2.5) squared
  // System values and counters
  double fKineticEnergy = 0;   // total kinetic energy
  double fPotentialEnergy = 0; // total potential energy
  Vector3D fNetMomentum;       // net momentum
  Vector3D rij;                // separation vector
  double fVirial = 0;          // virial
  double fVirialSum = 0;       // sum of virials
  int fVirialCount = 0;        // number of virial sums
  // Radial distribution parameters
  double dr = 0.01;                // resolution of histogram [r.u.]
  int fHistoFillCount = 0;         // count number of times histogram is filled
  double fArRadialDist[500] = {0}; // Ar radial distribution
  int nBins = 500;                 // number of bins in histogram
  // Mean-Squared Displacement (MSD) variables
  Ar::Atom *fAtomMSD; // Pointer to atom for MSD calculation
  Vector3D fR0[108];  // initial positions for MSD calculation

public:
  Sim() : fAtomHead{new Ar::Atom()}, fAtomMSD{new Ar::Atom()} {}
  ~Sim() { // avoid memory leaks
    Ar::Atom *atom;
    while (fAtomHead) {
      atom = fAtomHead;
      fAtomHead = fAtomHead->next;
      delete atom;
    }
    while (fAtomMSD) {
      atom = fAtomMSD;
      fAtomMSD = fAtomMSD->next;
      delete atom;
    }
  }

  // *******************************************************
  // Initializers
  // *******************************************************
  // Initialize velocities and positions of atoms. Optionally specify timestep
  void InitSim(Vector3D _R0[108], Vector3D _V0[108], double _TimeStep = 0);
  // Initialize atom positions for MSD calculation
  void InitMSD();
  void ResetPressure() {
    fVirialSum = 0;
    fVirialCount = 0;
  }

  // *******************************************************
  // Getters
  // *******************************************************
  double GetKineticEnergy() { return fKineticEnergy; }
  double GetPotentialEnergy() { return fPotentialEnergy; }
  Vector3D GetNetMomentum() { return fNetMomentum; }
  double GetTemperature() { return fKineticEnergy / 160.5; }
  double GetThermalWavelength() { return 0.9348779729 / sqrt(fKineticEnergy); }
  double GetVirial() { return fVirial; }
  double GetPressure() {
    return (108.0 * GetTemperature() + GetAveVirial() / 3.0) / (L * L * L);
  }
  double GetAveVirial() { return fVirialSum / fVirialCount; }
  double GetMeanSquaredDisplacement();
  Vector3D GetVelocityDistribution();
  Histo GetRadialDistributionHisto();
  double GetTimeStepSize() { return h; }
  Vector3D GetPositionOfAtom(int i);

  // *******************************************************
  // Steppers / Integrators
  // *******************************************************
  void VelocityVerletStep();
  void StepMSD();

  // *******************************************************
  // Utilities / Helper Functions
  // *******************************************************
  void ComputeForces();
  Vector3D MinimumImage(Vector3D rij);
  inline void ApplyPeriodicBC(Vector3D &r);
  void FillRadialDistributionHisto();

  // *******************************************************
  // Thermostats
  // *******************************************************
  void ScaleVelocities();    // set temperature to 120 K by fixing KE
  void GaussianThermostat(); // scale velocities w/ additional constraint force
  void AndersenThermostat();
  void LoweAnderson(); // conserves linear momentum (and angular?)
  void BerendsenThermostat(double tau); // weak coupling thermostat that allows
                                        // for small temperature fluctuations
  int CheckIfTempStable(double eps, int steps); // returns 1 if stable else 0
};

#endif