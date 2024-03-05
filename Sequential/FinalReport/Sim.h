#include "Ar/Atom.h"
#include "MyRNG.h"
#include "Vector3.h"
#include <algorithm> // std::copy
#include <iterator>  // std::begin, std::end

typedef struct RadialDistributionFunction {
  double bin[400];
  int nBins = 400;
  double binWidth = 0.01;
} Histo;

class NewSim {
private:
  int fNumAtoms = 108;            // number of atoms
  Ar::Atom *fAtomHead;            // Pointer to head of linked list of atoms
  Ar::Atom *fAtomMSD;             // Pointer to atom for MSD calculation
  Vector3D fR0[108];              // initial positions for MSD calculation
  Vector3D rij;                   // separation vector
  Vector3D fSimCell;              // Simulation cell
  double L = cbrt(108.0 / 0.8);   // length of simulation cell [r.u.]
  double fPotentialEnergy = 0;    // total potential energy
  double fVirial = 0;             // virial
  double fVirialSum = 0;          // sum of virials
  int fVirialCount = 0;           // number of virial sums
  double h = 0.0008311612;        // time step [r.u.]
  double fHistoResolution = 0.01; // resolution of histogram [r.u.]
  static const int nBins = 400;   // number of bins in histogram
  double fArRadialDist[nBins]; // histogram for Ar radial distribution function
  double fIGRadialDist[nBins]; // histogram for radial distribution of ideal gas

public:
  NewSim() : fAtomHead{new Ar::Atom()}, fAtomMSD{new Ar::Atom()} {}
  ~NewSim() { // avoid memory leaks
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

  double GetTimeStepSize() { return h; }
  double GetThermalWavelength() {
    return 0.0294393 * sqrt(2 * M_PI / GetTemperature()); // [r.u.]
  }
  double GetPotentialEnergy() { return fPotentialEnergy; }
  double GetVirial() { return fVirial; }
  void ResetPressure() {
    fVirialSum = 0;
    fVirialCount = 0;
  }
  void BookRadialDistributionHisto() {
    for (int i = 0; i < nBins; i++) {
      fArRadialDist[i] = 0;
      fIGRadialDist[i] = 0;
    }
  }
  void FillRadialDistributionHisto() {
    double dr = fHistoResolution;
    int bin;
    Ar::Atom *atom;
    for (atom = fAtomHead; atom->next; atom = atom->next) {
      for (Ar::Atom *atom2 = atom->next; atom2; atom2 = atom2->next) {
        double r = MinimumImage(atom->r - atom2->r).mag();
        bin = int(r / dr);
        fArRadialDist[bin] += 2;
      }
    }
  }
  Histo GetRadialDistributionHisto() {
    double dr = fHistoResolution;
    double n_ig;
    // fill the i.g. histogram
    for (int i = 0; i < nBins; ++i) {
      double r = i * dr;
      n_ig = 0.8 * 4 * M_PI * (pow(r + dr, 3.0) - pow(r, 3.0)) / 3.0;
      fIGRadialDist[i] = n_ig;
    }
    Histo histo;
    for (int i = 0; i < nBins; i++) {
      // normalize histogram to get time averages
      fArRadialDist[i] /= (fNumAtoms * h * fVirialCount);
      // fill histogram with ratio of Ar to ideal gas
      histo.bin[i] = fArRadialDist[i] / fIGRadialDist[i];
    }
    return histo;
  }
  double GetPressure() {
    return (108.0 * GetTemperature() + GetVirialAverage() / 3.0) / (L * L * L);
  }
  double GetVirialAverage() { return fVirialSum / fVirialCount; }
  double GetKineticEnergy() {
    double result = 0.0;
    Ar::Atom *atom;
    for (atom = fAtomHead; atom; atom = atom->next) {
      result += atom->v.sq();
    }
    return 0.5 * result;
  }
  Vector3D GetNetMomentum() {
    Vector3D result(0.0);
    Ar::Atom *atom;
    for (atom = fAtomHead; atom; atom = atom->next) {
      result += atom->v;
    }
    return result;
  }
  double GetTemperature() {
    return 2.0 * GetKineticEnergy() / (3.0 * (fNumAtoms - 1.0));
  }
  Vector3D GetMeanVelocity() {
    Vector3D result(0.0);
    Ar::Atom *atom;
    for (atom = fAtomHead; atom; atom = atom->next) {
      result += atom->v;
    }
    return result / 108.0;
  }
  Vector3D GetVelocityDistribution() {
    Vector3D mean = GetMeanVelocity();
    Vector3D numerator(0.0), denominator(0.0), tmp;
    Ar::Atom *atom;
    for (atom = fAtomHead; atom; atom = atom->next) {
      tmp = atom->v - mean;
      tmp.sq_inplace();
      denominator += tmp;
      tmp.sq_inplace();
      numerator += tmp;
    }
    numerator /= fNumAtoms;
    denominator /= fNumAtoms;
    denominator.sq_inplace();
    return Vector3D(numerator[0] / denominator[0],
                    numerator[1] / denominator[1],
                    numerator[2] / denominator[2]);
  }

  void InitSim(Vector3D _R0[108], Vector3D _V0[108]) {
    fSimCell = Vector3D(L);

    *fAtomHead = Ar::Atom(_R0[0], _V0[0], 0);
    Ar::Atom *NewAtom = new Ar::Atom(_R0[1], _V0[1], 1);
    fAtomHead->next = NewAtom;
    for (int i = 2; i < fNumAtoms; ++i) {
      NewAtom->next = new Ar::Atom(_R0[i], _V0[i], i);
      NewAtom = NewAtom->next;
    }
    ComputeForces();
  }

  void InitMSD() {
    *fAtomMSD = Ar::Atom(fAtomHead->r, fAtomHead->v, fAtomHead->ID);
    Ar::Atom *NewAtom = new Ar::Atom(fAtomHead->next->r, fAtomHead->next->v,
                                     fAtomHead->next->ID);
    fAtomMSD->next = NewAtom;
    Ar::Atom *input = fAtomHead->next->next;
    fR0[0] = fAtomHead->r;
    fR0[1] = fAtomHead->next->r;
    int i = 2;
    while (input) {
      fR0[i] = input->r;
      ++i;
      NewAtom->next = new Ar::Atom(input->r, input->v, input->ID);
      NewAtom = NewAtom->next;
      input = input->next;
    }
  }

  void StepMSD() {
    Ar::Atom *atom_i;
    for (atom_i = fAtomMSD; atom_i; atom_i = atom_i->next) {
      atom_i->a = Vector3D(0.0);
    }
    Ar::Atom *atom_j;
    for (atom_i = fAtomMSD; atom_i->next; atom_i = atom_i->next) {
      for (atom_j = atom_i->next; atom_j; atom_j = atom_j->next) {
        rij = atom_i->r - atom_j->r;
        double rij2 = rij.sq();
        double rij2i = 1.0 / rij2;
        double rij6i = rij2i * rij2i * rij2i;
        double funcVal = 48.0 * rij6i * (rij6i - 0.5) * rij2i;
        atom_i->a += rij * funcVal;
        atom_j->a -= rij * funcVal;
      }
      atom_i->r += h * atom_i->v + 0.5 * h * h * atom_i->a;
      atom_i->v += 0.5 * h * atom_i->a;
    }
    for (atom_i = fAtomMSD; atom_i->next; atom_i = atom_i->next) {
      for (atom_j = atom_i->next; atom_j; atom_j = atom_j->next) {
        rij = atom_i->r - atom_j->r;
        double rij2 = rij.sq();
        double rij2i = 1.0 / rij2;
        double rij6i = rij2i * rij2i * rij2i;
        double funcVal = 48.0 * rij6i * (rij6i - 0.5) * rij2i;
        atom_i->a += rij * funcVal;
        atom_j->a -= rij * funcVal;
      }
      atom_i->v += 0.5 * h * atom_i->a;
    }
  }

  double GetMeanSquaredDisplacement() {
    double result = 0.0;
    Ar::Atom *atom_i;
    for (atom_i = fAtomMSD; atom_i; atom_i = atom_i->next) {
      result += (atom_i->r - fR0[atom_i->ID]).sq();
    }
    result /= 108.0;
    return result;
  }

  void ComputeForces() {
    // zero them out
    fPotentialEnergy = 0.0;
    fVirial = 0.0;
    Ar::Atom *atom_i;
    Ar::Atom *atom_j;
    for (atom_i = fAtomHead; atom_i != nullptr; atom_i = atom_i->next) {
      atom_i->a = Vector3D(0.0);
    }
    // compute forces, virial, and system potential energy
    for (atom_i = fAtomHead; atom_i != nullptr; atom_i = atom_i->next) {
      for (atom_j = atom_i->next; atom_j != nullptr; atom_j = atom_j->next) {
        rij = MinimumImage(atom_i->r - atom_j->r);
        // compute squared distance and skip if greater than cutoff
        double rij2 = rij.sq();
        // NOTE: This improves speed but introduces additional error
        // if (rij2 < CutoffRadiusSq) {
        // compute force
        double rij2i = 1.0 / rij2;
        double rij6i = rij2i * rij2i * rij2i;
        double funcVal = 48.0 * rij6i * (rij6i - 0.5) * rij2i;
        atom_i->a += rij * funcVal;
        atom_j->a -= rij * funcVal;
        // compute potential energy
        fPotentialEnergy += 4.0 * rij6i * (rij6i - 1.0);
        // compute virial
        fVirial += funcVal * rij2;
        fVirialSum += fVirial;
        fVirialCount += 1;
        //}
      }
    }
  }

  inline void ApplyPeriodicBC(Vector3D &r) {
    r.x() -= (copysign(0.5 * L, r.x() - L) + copysign(0.5 * L, r.x()));
    r.y() -= (copysign(0.5 * L, r.y() - L) + copysign(0.5 * L, r.y()));
    r.z() -= (copysign(0.5 * L, r.z() - L) + copysign(0.5 * L, r.z()));
  }

  Vector3D MinimumImage(Vector3D rij) {
    rij.x() -= (copysign(0.5 * L, rij.x() - 0.5 * L) +
                copysign(0.5 * L, rij.x() + 0.5 * L));
    rij.y() -= (copysign(0.5 * L, rij.y() - 0.5 * L) +
                copysign(0.5 * L, rij.y() + 0.5 * L));
    rij.z() -= (copysign(0.5 * L, rij.z() - 0.5 * L) +
                copysign(0.5 * L, rij.z() + 0.5 * L));
    return rij;
  }

  void VelocityVerletStep() {
    // ComputeForces();
    Ar::Atom *atom;
    for (atom = fAtomHead; atom != nullptr; atom = atom->next) {
      // update positions
      atom->r += h * atom->v + 0.5 * h * h * atom->a; // error ~ O(h^3)
      // apply periodic boundary conditions
      ApplyPeriodicBC(atom->r);
      // update velocities
      atom->v += atom->a * 0.5 * h;
    }
    // compute forces
    ComputeForces();
    for (atom = fAtomHead; atom != nullptr; atom = atom->next) {
      // update velocities
      atom->v += atom->a * 0.5 * h;
    }
  }

  // *******************************************************
  // Thermostats
  // *******************************************************

  void RescaleTemperature() {
    double CorrectionFactor = sqrt(1.0 / GetTemperature());
    Ar::Atom *atom;
    for (atom = fAtomHead; atom != nullptr; atom = atom->next) {
      atom->v *= CorrectionFactor;
    }
  }

  void AndersenThermostat() {
    // init rng
    init_random();

    double gamma = 1.0; // heat bath stochastic collision frequency

    // loop over atoms
    Ar::Atom *atom;
    for (atom = fAtomHead; atom; atom = atom->next) {
      // update positions
      atom->r += h * atom->v + 0.5 * h * h * atom->a; // error ~ O(h^3)
      // apply periodic boundary conditions
      ApplyPeriodicBC(atom->r);
      // update velocities
      atom->v += atom->a * 0.5 * h;
    }
    ComputeForces();
    for (atom = fAtomHead; atom; atom = atom->next) {
      // update velocities
      atom->v += atom->a * 0.5 * h;
      // update velocities
      atom->v += atom->a * 0.5 * h;
      // von neumann rejection method
      if (rndu() >= gamma * h) {
        // generate random velocity vector from gaussian
        atom->v = Vector3D(rnd(), rnd(), rnd());
      }
    }
    // save rng state
    save_random();
  }

  // conserves linear momentum (and angular?)
  void LoweAnderson() {
    // init rng
    init_random();

    double gamma = 20.0; // heat bath stochastic collision frequency

    // loop over atoms
    Ar::Atom *atom;
    // Pre-collision with heat bath step
    for (atom = fAtomHead; atom; atom = atom->next) {
      // update positions
      atom->r += h * atom->v + 0.5 * h * h * atom->a; // error ~ O(h^3)
      // apply periodic boundary conditions
      ApplyPeriodicBC(atom->r);
      // update velocities
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
      // update velocities
      atom->v += atom->a * 0.5 * h;
    }
    // save rng state
    save_random();
  }

  void RescaleVelocities() {
    // loop over atoms
    Ar::Atom *atom;
    // Pre-collision with heat bath step
    for (atom = fAtomHead; atom; atom = atom->next) {
      // update positions
      atom->r += h * atom->v + 0.5 * h * h * atom->a; // error ~ O(h^3)
      // apply periodic boundary conditions
      ApplyPeriodicBC(atom->r);
      // update velocities
      atom->v += atom->a * 0.5 * h;
    }
    ComputeForces();
    for (atom = fAtomHead; atom; atom = atom->next) {
      // update velocities
      atom->v += atom->a * 0.5 * h;
    }
    // rescale velocities
    RescaleTemperature();
  }

  int CheckIfTempStable(double eps, int steps) {
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
        // update positions
        atoms[i].r +=
            h * atoms[i].v + 0.5 * h * h * atoms[i].a; // error ~ O(h^3)
        // apply periodic boundary conditions
        ApplyPeriodicBC(atoms[i].r);
        // update velocities
        atoms[i].v += atoms[i].a * 0.5 * h;
      }
      // compute forces
      for (int i = 0; i < 108; i++) {
        atoms[i].a = Vector3D(0.0);
      }
      // compute forces, virial, and system potential energy
      for (int iatm = 0; iatm < 108 - 1; iatm++) {
        for (int jatm = iatm + 1; jatm < 108; jatm++) {
          rij = MinimumImage(atoms[iatm].r - atoms[jatm].r);
          // compute squared distance and skip if greater than cutoff
          double rij2 = rij.sq();
          // NOTE: This improves speed but introduces additional error
          // if (rij2 < CutoffRadiusSq) {
          // compute force
          double rij2i = 1.0 / rij2;
          double rij6i = rij2i * rij2i * rij2i;
          double funcVal = 48.0 * rij6i * (rij6i - 0.5) * rij2i;
          atoms[iatm].a += rij * funcVal;
          atoms[jatm].a -= rij * funcVal;
        }
        // update velocities
        atoms[iatm].v += atoms[iatm].a * 0.5 * h;
        // compute temperature
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
};
