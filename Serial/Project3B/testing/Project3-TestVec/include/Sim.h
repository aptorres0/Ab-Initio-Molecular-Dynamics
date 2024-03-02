#ifndef _SIM_H
#define _SIM_H

#include "Ar/Ensemble.h"
#include "vec.h"

typedef struct {
  int NumberOfDimensions;
  int NumberOfAtoms;
  Vector3D SimulationCellBounds;
  double TimeStepSize;
  double (*PotentialEnergyFunction)(const Vector3D &);
  Vector3D (*ForceFunction)(const Vector3D &);
} SimArgs;

/**
 * Class to implement time evolution of system. Actual updating of positions and
 * velocities is done in the Sim::VelocityVerletStep function and updating of
 * accelerations is done in Sim::ComputeForces function.
 */
class Sim {
private:
  double (*U_ij)(const Vector3D &);   // Pointer to potential energy function
  Vector3D (*F_ij)(const Vector3D &); // Pointer to force function
  int fDIM;                           // dimensions of system
  int fNumAtoms;                      // number of atoms
  Vector3D fSimBox;                   // lengths of each side of simulation box
  Ar::Ensemble fEnsemble; // collection of fNumAtoms with vectors for position,
                          // velocity, and acceletation
  Vector3D rij{fDIM}; // separation vector declared once to avoid allocations
                      // and deletions
  Vector3D _force{
      fDIM}; // force vector declared once to avoid allocations and deletions
  double h;  // Time step size

public:
  Sim() : fDIM{3}, fNumAtoms{108}, h{0.01} {}

  void InitSim(SimArgs args);

  void InitSim(int _DIM, int _NumAtoms, double (*_U)(const Vector3D &),
               Vector3D (*_F)(const Vector3D &), Vector3D _SimBox, double _h);

  // Initialize position of atoms to FCC crystal lattice
  void InitPositions();

  /**
   * Steps to initializing the velocities:
   *  1. Set to random values from gaussian centered at 0 with unitless
   * variance of 1
   *  2. Apply correction to set system net momentum to 0
   *  3. Apply correction to set system temperature to 120 K (or 1 in
   * characteristic units)
   */
  void InitVelocities();

  // Compute net momentum of the system
  Vector3D GetNetMomentum() const;

  // Compute temperature of the system
  double GetTemperature() const;

  void ApplyPeriodicBC();

  // copysign function for my vector class
  Vector3D copysign_v(const Vector3D &_lhs, const Vector3D &_rhs);

  // check if periodic BC is necessary
  bool IsInSimBox(const Vector3D &check);

  const Vector3D &MinimumImage(const Vector3D &ri, const Vector3D &rj);

  void ComputeForces();

  void VelocityVerletStep();

  void PrintSample();

  double GetPotentialEnergy();

  double GetKineticEnergy();

  double GetInternalEnergy();
};

#endif