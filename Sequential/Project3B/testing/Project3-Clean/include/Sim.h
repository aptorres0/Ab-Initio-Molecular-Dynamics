#ifndef _SIM_H
#define _SIM_H

#include "Ar/Ensemble.h"
#include "vec.h"

typedef struct {
  int NumberOfDimensions;
  int NumberOfAtoms;
  VectorD SimulationCellBounds;
  double TimeStepSize;
  double (*PotentialEnergyFunction)(const VectorD &);
  VectorD (*ForceFunction)(const VectorD &);
} SimArgs;

/**
 * Class to implement time evolution of system. Actual updating of positions and
 * velocities is done in the Sim::VelocityVerletStep function and updating of
 * accelerations is done in Sim::ComputeForces function.
 */
class Sim {
private:
  double (*U_ij)(const VectorD &);  // Pointer to potential energy function
  VectorD (*F_ij)(const VectorD &); // Pointer to force function
  int fDIM;                         // dimensions of system
  int fNumAtoms;                    // number of atoms
  VectorD fSimBox;                  // lengths of each side of simulation box
  Ar::Ensemble fEnsemble; // collection of fNumAtoms with vectors for position,
                          // velocity, and acceletation
  VectorD rij{fDIM}; // separation vector declared once to avoid allocations and
                     // deletions
  VectorD _force{
      fDIM}; // force vector declared once to avoid allocations and deletions
  double h;  // Time step size

public:
  Sim() : fDIM{3}, fNumAtoms{108}, h{0.01} {}

  void InitSim(SimArgs args);

  void InitSim(int _DIM, int _NumAtoms, double (*_U)(const VectorD &),
               VectorD (*_F)(const VectorD &), VectorD _SimBox, double _h);

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
  VectorD GetNetMomentum() const;

  // Compute temperature of the system
  double GetTemperature() const;

  void ApplyPeriodicBC();

  // copysign function for my vector class
  VectorD copysign_v(const VectorD &_lhs, const VectorD &_rhs);

  // check if periodic BC is necessary
  bool IsInSimBox(const VectorD &check);

  bool ImageIsCloser(const VectorD &check);

  const VectorD &MinimumImage(const VectorD &ri, const VectorD &rj);

  void ComputeForces();

  void VelocityVerletStep();

  void PrintSample();

  double GetPotentialEnergy();

  double GetKineticEnergy();

  double GetInternalEnergy();
};

#endif