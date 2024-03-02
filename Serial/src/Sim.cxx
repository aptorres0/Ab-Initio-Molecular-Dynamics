/**
 * Sim.cxx
 *
 * Author: Alexander Paul Torres
 * Date: 14 APR 23
 * Class: PHYS 7411 - Computational Physics I
 *
 * Description:
 * This file contains the implementation of the Sim class member functions.
 *
 */
#include "../include/Sim.h"
#include "../include/MyRNG.h"

// *******************************************************
// Initializers
// *******************************************************
void Sim::InitSim(Vector3D _R0[108], Vector3D _V0[108], double _TimeStep) {
  if (_TimeStep)
    h = _TimeStep;
  Vector3D A0 = Vector3D(0.0); // initialize acceleration to 0
  *fAtomHead = Ar::Atom(_R0[0], _V0[0], A0, 0);
  Ar::Atom *NewAtom = new Ar::Atom(_R0[1], _V0[1], A0, 1);
  fAtomHead->next = NewAtom;
  for (int i = 2; i < fNumAtoms; ++i) {
    NewAtom = (NewAtom->next = new Ar::Atom(_R0[i], _V0[i], A0, i));
    fNetMomentum += _V0[i];
    fKineticEnergy += _V0[i].sq();
  }
  fNetMomentum = _V0[0] + _V0[1];
  fKineticEnergy += _V0[0].sq() + _V0[1].sq();
  fKineticEnergy *= 0.5;
  ComputeForces(); // initialize acceleration vectors of atoms
}

void Sim::InitMSD() {
  // Initialize head of MSD linked list
  *fAtomMSD = Ar::Atom(fAtomHead->r, fAtomHead->v, fAtomHead->a, fAtomHead->ID);
  Ar::Atom *input = fAtomHead;
  fR0[input->ID] = input->r;
  // Loop over atoms and add to MSD linked list and set initial positions
  input = input->next;
  Ar::Atom *NextAtom = new Ar::Atom(input->r, input->v, input->a, input->ID);
  fAtomMSD->next = NextAtom;
  fR0[input->ID] = input->r;
  while ((input = input->next)) {
    fR0[input->ID] = input->r;
    NextAtom = (NextAtom->next =
                    new Ar::Atom(input->r, input->v, input->a, input->ID));
  }
}

// *******************************************************
// Getters
// *******************************************************
double Sim::GetMeanSquaredDisplacement() {
  double result = 0.0;
  for (Ar::Atom *atom = fAtomMSD; atom; atom = atom->next) {
    result += (atom->r - fR0[atom->ID]).sq();
  }
  result /= fNumAtoms;
  return result;
}

Histo Sim::GetRadialDistributionHisto() {
  double n_ig, r;
  Histo histo;
  for (int i = 0; i < nBins; i++) {
    r = i * dr;
    n_ig = 0.8 * 4 * M_PI * (pow(r + dr, 3.0) - pow(r, 3.0)) / 3.0;
    // normalize histogram to get time averages and fill with ratio of Ar to ig
    histo.bin[i] = fArRadialDist[i] / (n_ig * fNumAtoms * fHistoFillCount);
  }
  return histo;
}

Vector3D Sim::GetVelocityDistribution() {
  Vector3D mean = GetNetMomentum() / 108.0;
  Vector3D numerator(0.0), denominator(0.0), tmp;
  for (Ar::Atom *atom = fAtomHead; atom; atom = atom->next) {
    tmp = atom->v - mean;
    tmp.sq_inplace();
    denominator += tmp;
    tmp.sq_inplace();
    numerator += tmp;
  }
  numerator /= fNumAtoms;
  denominator /= fNumAtoms;
  denominator.sq_inplace();
  return Vector3D(numerator[0] / denominator[0], numerator[1] / denominator[1],
                  numerator[2] / denominator[2]);
}

Vector3D Sim::GetPositionOfAtom(int i) {
  Ar::Atom *atom = fAtomHead;
  while (atom->ID != i) {
    atom = atom->next;
  }
  if (atom->ID != i) {
    std::cout << "!!! ERROR: Atom " << i << " not found in linked list. !!!\n";
    exit(1);
  }
  return Vector3D(atom->r); // return copy of position
}

// *******************************************************
// Steppers / Integrators
// *******************************************************
void Sim::VelocityVerletStep() {
  for (Ar::Atom *atom = fAtomHead; atom; atom = atom->next) {
    atom->v += atom->a * 0.5 * h;
    atom->r += h * atom->v; // error ~ O(h^3)
    ApplyPeriodicBC(atom->r);
    // zero out acceleration before updating
    atom->a = Vector3D(0.0);
  }
  ComputeForces();
  fNetMomentum = Vector3D(0.0);
  fKineticEnergy = 0.0;
  for (Ar::Atom *atom = fAtomHead; atom; atom = atom->next) {
    atom->v += atom->a * 0.5 * h;
    fNetMomentum += atom->v;
    fKineticEnergy += atom->v.sq();
  }
  fKineticEnergy *= 0.5;
}

void Sim::StepMSD() {
  for (Ar::Atom *atom = fAtomMSD; atom; atom = atom->next) {
    atom->r += h * atom->v + 0.5 * h * h * atom->a;
    atom->v += 0.5 * h * atom->a;
    // zero out acceleration before updating
    atom->a = Vector3D(0.0);
  }
  // update accelerations and velocities
  for (Ar::Atom *atom_i = fAtomMSD; atom_i->next; atom_i = atom_i->next) {
    for (Ar::Atom *atom_j = atom_i->next; atom_j; atom_j = atom_j->next) {
      rij = atom_i->r - atom_j->r;
      double rij2i = 1.0 / rij.sq();
      double rij6i = rij2i * rij2i * rij2i;
      double funcVal = 48.0 * rij6i * rij2i * (rij6i - 0.5);
      atom_i->a += rij * funcVal;
      atom_j->a -= rij * funcVal;
    }
    atom_i->v += 0.5 * h * atom_i->a;
  }
}

// *******************************************************
// Helper Functions
// *******************************************************
void Sim::ComputeForces() {
  // compute forces, virial, and system potential energy
  fPotentialEnergy = 0.0;
  fVirial = 0.0;
  double rij2i, rij6i;
  for (Ar::Atom *atom_i = fAtomHead; atom_i->next; atom_i = atom_i->next) {
    for (Ar::Atom *atom_j = atom_i->next; atom_j; atom_j = atom_j->next) {
      rij = MinimumImage(atom_i->r - atom_j->r);
      // NOTE: This improves speed but introduces additional error
      // if (rij.sq() < CutoffRadiusSq) {
      // compute force
      rij2i = 1.0 / rij.sq();
      rij6i = rij2i * rij2i * rij2i;
      rij *= 48.0 * rij6i * rij2i * (rij6i - 0.5); // converting to a force
      atom_i->a += rij;
      atom_j->a -= rij;
      // compute potential energy
      fPotentialEnergy += 4.0 * rij6i * (rij6i - 1.0);
      // compute virial
      fVirial += 24.0 * rij6i * (2.0 * rij6i - 1.0);
    }
    //}
  }
  fVirialSum += fVirial;
  fVirialCount += 1;
}

Vector3D Sim::MinimumImage(Vector3D rij) {
  rij.x() -= (copysign(0.5 * L, rij.x() - 0.5 * L) +
              copysign(0.5 * L, rij.x() + 0.5 * L));
  rij.y() -= (copysign(0.5 * L, rij.y() - 0.5 * L) +
              copysign(0.5 * L, rij.y() + 0.5 * L));
  rij.z() -= (copysign(0.5 * L, rij.z() - 0.5 * L) +
              copysign(0.5 * L, rij.z() + 0.5 * L));
  return rij;
}

inline void Sim::ApplyPeriodicBC(Vector3D &r) {
  r.x() -= (copysign(0.5 * L, r.x() - L) + copysign(0.5 * L, r.x()));
  r.y() -= (copysign(0.5 * L, r.y() - L) + copysign(0.5 * L, r.y()));
  r.z() -= (copysign(0.5 * L, r.z() - L) + copysign(0.5 * L, r.z()));
}

void Sim::FillRadialDistributionHisto() {
  int bin;
  for (Ar::Atom *atom = fAtomHead; atom->next; atom = atom->next) {
    for (Ar::Atom *atom2 = atom->next; atom2; atom2 = atom2->next) {
      double r = MinimumImage(atom->r - atom2->r).mag();
      bin = int(r / dr);
      fArRadialDist[bin] += 2;
    }
  }
  fHistoFillCount++;
}

// *******************************************************
// Thermostats
// *******************************************************
// Strong coupling of the heat bath to the system.
void Sim::ScaleVelocities() {
  for (Ar::Atom *atom = fAtomHead; atom; atom = atom->next) {
    atom->r += h * atom->v + 0.5 * h * h * atom->a; // error ~ O(h^3)
    ApplyPeriodicBC(atom->r);
    atom->v += atom->a * 0.5 * h;
    // zero out before updating
    atom->a = Vector3D(0.0);
  }
  ComputeForces();
  // NOTE: Updating the kinetic energy here is required to get the right
  // temperature correction factor
  fKineticEnergy = 0.0;
  for (Ar::Atom *atom = fAtomHead; atom; atom = atom->next) {
    atom->v += atom->a * 0.5 * h;
    fKineticEnergy += atom->v.sq();
  }
  fKineticEnergy *= 0.5;
  double CorrectionFactor = 1.0 / sqrt(GetTemperature());
  fKineticEnergy = 0.0;
  fNetMomentum = Vector3D(0.0);
  for (Ar::Atom *atom = fAtomHead; atom; atom = atom->next) {
    atom->v *= CorrectionFactor;
    fKineticEnergy += atom->v.sq();
    fNetMomentum += atom->v;
  }
  fKineticEnergy *= 0.5;
}

void Sim::GaussianThermostat() {
  Vector3D a_old[108];
  Vector3D p_old[108];
  for (Ar::Atom *atom = fAtomHead; atom; atom = atom->next) {
    a_old[atom->ID] = atom->a;
    p_old[atom->ID] = atom->v;
    atom->v += atom->a * 0.5 * h;
    atom->r += h * atom->v; // error ~ O(h^3)
    ApplyPeriodicBC(atom->r);
    // zero out before updating
    atom->a = Vector3D(0.0);
  }
  ComputeForces();
  // compute alpha term
  double alpha = 0.0, NetVelo2 = 0.0;
  Vector3D AveA(0.0);
  for (Ar::Atom *atom = fAtomHead; atom; atom = atom->next) {
    AveA = (atom->a + a_old[atom->ID]) * 0.5;
    alpha -= dot(atom->v, AveA);
    NetVelo2 += atom->v.sq();
  }
  alpha /= NetVelo2;
  fKineticEnergy = 0.0;
  fNetMomentum = Vector3D(0.0);
  // apply correction factor and compute resulting momentum and kinetic energy
  for (Ar::Atom *atom = fAtomHead; atom; atom = atom->next) {
    atom->v -= alpha * p_old[atom->ID];
    fKineticEnergy += atom->v.sq();
    fNetMomentum += atom->v;
  }
  fKineticEnergy *= 0.5;
}

void Sim::AndersenThermostat() {
  init_random();
  double gamma = 1.0; // heat bath stochastic collision frequency
  Ar::Atom *atom;
  for (atom = fAtomHead; atom; atom = atom->next) {
    atom->r += h * atom->v + 0.5 * h * h * atom->a; // error ~ O(h^3)
    ApplyPeriodicBC(atom->r);
    atom->v += atom->a * 0.5 * h;
  }
  ComputeForces();
  for (atom = fAtomHead; atom; atom = atom->next) {
    atom->v += atom->a * 0.5 * h;
    atom->v += atom->a * 0.5 * h;
    // von neumann rejection method
    if (rndu() >= gamma * h) {
      // generate random velocity vector from gaussian
      atom->v = Vector3D(rnd(), rnd(), rnd());
    }
  }
  save_random();
}

// conserves linear momentum (and angular?)
void Sim::LoweAnderson() {
  init_random();
  double gamma = 20.0; // heat bath stochastic collision frequency
  Ar::Atom *atom;
  for (atom = fAtomHead; atom; atom = atom->next) {
    atom->r += h * atom->v + 0.5 * h * h * atom->a; // error ~ O(h^3)
    ApplyPeriodicBC(atom->r);
    atom->v += atom->a * 0.5 * h;
  }
  ComputeForces();
  for (atom = fAtomHead; atom; atom = atom->next) {
    Ar::Atom *atom2;
    for (atom2 = atom->next; atom2; atom2 = atom2->next) {
      // von neumann rejection method
      if (rndu() < gamma * h && 2.0 > (atom->r - atom2->r).sq()) {
        double lambda = rnd() * sqrt(2.0); // heat bath stochastic variable
        Vector3D vij = atom->v - atom2->v;
        Vector3D uij = (atom->r - atom2->r).unit_vec();
        atom->v += (0.5 * (lambda - dot(vij, uij)) * uij);
        atom2->v -= (0.5 * (lambda - dot(vij, uij)) * uij);
      }
    }
    atom->v += atom->a * 0.5 * h;
  }
  save_random();
}

void Sim::BerendsenThermostat(double tau = 0.1) {
  double CorrectionFactor = sqrt(1.0 + h / tau * (1. / GetTemperature() - 1.0));
  fKineticEnergy = 0.0;
  fNetMomentum = Vector3D(0.0);
  for (Ar::Atom *atom = fAtomHead; atom; atom = atom->next) {
    atom->v *= CorrectionFactor;
    fKineticEnergy += atom->v.sq();
    fNetMomentum += atom->v;
  }
  fKineticEnergy *= 0.5;
}

int Sim::CheckIfTempStable(double eps, int steps) {
  double Temp = GetTemperature();
  Ar::Atom atoms[108];
  Ar::Atom *in;
  for (in = fAtomHead; in; in = in->next) {
    atoms[in->ID] = *in;
  }
  double Tlocal = 0.0;
  for (int istep = 0; istep <= steps; istep++) {
    double Tstep = 0.0;
    for (int i = 0; i < 108; ++i) {
      atoms[i].r += h * atoms[i].v + 0.5 * h * h * atoms[i].a; // error ~ O(h^3)
      ApplyPeriodicBC(atoms[i].r);
      atoms[i].v += atoms[i].a * 0.5 * h;
    }
    for (int i = 0; i < 108; i++) {
      atoms[i].a = Vector3D(0.0);
    }
    for (int iatm = 0; iatm < 108 - 1; iatm++) {
      for (int jatm = iatm + 1; jatm < 108; jatm++) {
        rij = MinimumImage(atoms[iatm].r - atoms[jatm].r);
        // compute squared distance and skip if greater than cutoff
        // NOTE: This improves speed but introduces additional error
        double rij2 = rij.sq();
        if (rij2 < CutoffRadiusSq) {
          double rij2i = 1.0 / rij2;
          double rij6i = rij2i * rij2i * rij2i;
          double funcVal = 48.0 * rij6i * (rij6i - 0.5) * rij2i;
          atoms[iatm].a += rij * funcVal;
          atoms[jatm].a -= rij * funcVal;
        }
      }
      atoms[iatm].v += atoms[iatm].a * 0.5 * h;
      Tstep += atoms[iatm].v.sq();
    }
    Tstep /= (3.0 * (108.0 - 1.0));
    if (istep == steps) {
      Tlocal = Tstep;
    }
  }
  if (fabs(Temp - Tlocal) < eps) {
    std::cout << "Temperature is stable to within "
              << 100 * (Temp - Tlocal) / Temp << "% after " << steps
              << " steps.\n";
    return 1;
  }
  return 0;
}