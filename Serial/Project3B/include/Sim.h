#ifndef _SIM_H
#define _SIM_H

#include "Particle.h"
#include "Vector3.h"

typedef struct {
  double BoxLength;
  double TimeStepSize;
  double (*PotentialEnergyFunction)(const Vector3D &);
  Vector3D (*ForceFunction)(const Vector3D &);
} SimArgs;

class Sim {
private:
  int fNumAtoms = 108;                // number of atoms
  double (*U_ij)(const Vector3D &);   // Pointer to potential energy function
  Vector3D (*F_ij)(const Vector3D &); // Pointer to force function
  Particle *fParticleHead;            // Pointer to particle
  Vector3D rij{3};   // separation vector declared once to avoid allocations
                     // and deletions
  double h;          // Time step size = 1% of the period of oscillation for the
                     // harmonic oscillator in the lennard jones potential
  double fBoxLength; // Length of the simulation box
  Vector3D _force;
  Vector3D fSimBox;
  Vector3D fSimBoxHalf;
  double fVirial = 0.0;
  double fViralSum = 0.0;
  double fPressure = 0.0;
  Vector3D fR0[108];            // used for mean square dispacements
  Particle *fParticleHeadNoPBC; // used for mean square dispacements

public:
  Sim() : fParticleHead{new Particle()} {}

  void InitSim(SimArgs args);

  // Initialize position of atoms to FCC crystal lattice
  void InitPositions();

  void PartA();

  void Question11A();

  void InitR0();

  double GetMeanSquareDisplacement();

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
  Vector3D copysign_vector3(const Vector3D &_lhs, const Vector3D &_rhs);

  // check if periodic BC is necessary
  bool IsInSimBox(const Vector3D &check);

  bool ImageIsCloser(const Vector3D &);

  const Vector3D &MinimumImage(const Vector3D &ri, const Vector3D &rj);

  void ComputeForces();

  void VelocityVerletStep();

  void ComputeForcesNoPBC();

  void VelocityVerletStepNoPBC();

  double GetPotentialEnergy();

  double GetKineticEnergy();

  double GetInternalEnergy();

  void CorrectTemperature();

  Vector3D GetMeanVelocity() const;

  double Getfvelo(int j) const;
  double Getfx() const { return Getfvelo(0); }
  double Getfy() const { return Getfvelo(1); }
  double Getfz() const { return Getfvelo(2); }
  Vector3D GetVelocityDistribution() const {
    return Vector3D(Getfx(), Getfy(), Getfz());
  }

  Vector3D GetPos(int i) const {
    if (i == 0)
      return fParticleHead->r;
    Particle *entry = fParticleHead->next;
    for (int j = 1; entry->next; ++j) {
      if (j == i) {
        return entry->r;
      }
      entry = entry->next;
    }
    return entry->next->r;
  }

  void Equilibrate(int num_steps);

  void ComputeVirial();

  void UpdatePressure(int num_steps);

  double GetPressure() const { return fPressure; }

  void PhaseI(int num_steps);
  void PhaseII(int num_steps);
};

#endif