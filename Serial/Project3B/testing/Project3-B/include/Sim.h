#ifndef _SIM_H
#define _SIM_H

#include "Ar/Ensemble.h"
#include "Vector3.h"

typedef struct {
  Vector3F SimulationCellBounds;
  float TimeStepSize;
  float (*PotentialEnergyFunction)(const Vector3F &);
  Vector3F (*ForceFunction)(const Vector3F &);
} SimArgs;

/**
 * Class to implement time evolution of system. Actual updating of positions and
 * velocities is done in the Sim::VelocityVerletStep function and updating of
 * accelerations is done in Sim::ComputeForces function.
 */
class Sim {
private:
  int fNumAtoms = 108;                // number of atoms
  float (*U_ij)(const Vector3F &);    // Pointer to potential energy function
  Vector3F (*F_ij)(const Vector3F &); // Pointer to force function
  Vector3F fSimBox;                   // lengths of each side of simulation box
  Vector3F fSimBoxHalf;               // half of the lengths of each side of
                                      // simulation box
  Ar::Ensemble fEnsemble; // collection of fNumAtoms with vectors for position,
                          // velocity, and acceletation
  Vector3F rij{3}; // separation vector declared once to avoid allocations
                   // and deletions
  Vector3F _force{
      3};  // force vector declared once to avoid allocations and deletions
  float h; // Time step size = 1% of the period of oscillation for the
           // harmonic oscillator in the lennard jones potential
  float LHalf;

public:
  Sim() {}

  void InitSim(SimArgs args);

  Ar::Ensemble &GetEnsemble() { return fEnsemble; }

  void ParallelComputeForces();

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
  Vector3F GetNetMomentum() const;

  // Compute temperature of the system
  float GetTemperature() const;

  void ApplyPeriodicBC();

  // copysign function for my vector class
  Vector3F copysign_vector3(const Vector3F &_lhs, const Vector3F &_rhs);

  // check if periodic BC is necessary
  bool IsInSimBox(const Vector3F &check);

  bool ImageIsCloser(const Vector3F &);

  const Vector3F &MinimumImage(const Vector3F &ri, const Vector3F &rj);

  void ComputeForces();

  void VelocityVerletStep();

  void PrintSample();

  float GetPotentialEnergy();

  float GetKineticEnergy();

  float GetInternalEnergy();

  void CorrectTemperature();

  Vector3F GetMeanVelocity() const;

  float Getfvelo(int j) const;
  float Getfx() const { return Getfvelo(0); }
  float Getfy() const { return Getfvelo(1); }
  float Getfz() const { return Getfvelo(2); }
  Vector3F GetVelocityDistribution() const {
    return Vector3F(Getfx(), Getfy(), Getfz());
  }

  float GetAtomX(int i) const { return fEnsemble.r(i).x(); }
  float GetAtomY(int i) const { return fEnsemble.r(i).y(); }
  float GetAtomZ(int i) const { return fEnsemble.r(i).z(); }
  Vector3F GetPos(int i) const { return fEnsemble.r(i); }

  // static void *ComputeForcesThreaded(void *_args);
};

#endif