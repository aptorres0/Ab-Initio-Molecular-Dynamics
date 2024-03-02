#include "../include/Sim.h"

Vector3F gFCC[108];
Vector3F gV0[108];
Vector3F FCC(int i) { return gFCC[i]; }
Vector3F V0(int i) { return gV0[i]; }
void InitFCC(float l) {
  // Generate position vectors for the 3D Volume of 9 unit cells
  for (int i = 0; i < 3; i++) { // loop over cells in x
    float xmin = i * l;
    for (int j = 0; j < 3; j++) { // loop over cells in y
      float ymin = j * l;
      for (int k = 0; k < 3; k++) { // loop over cells in z
        float zmin = k * l;
        gFCC[i * 36 + j * 12 + k * 4] = {xmin, ymin, zmin};
        gFCC[i * 36 + j * 12 + k * 4 + 1] = {xmin, 0.5f * l + ymin,
                                             0.5f * l + zmin};
        gFCC[i * 36 + j * 12 + k * 4 + 2] = {0.5f * l + xmin, ymin,
                                             0.5f * l + zmin};
        gFCC[i * 36 + j * 12 + k * 4 + 3] = {0.5f * l + xmin, 0.5f * l + ymin,
                                             zmin};
      }
    }
  }
}
void InitV0() {
  /**
   * Steps to initializing the velocities:
   *  1. Set to random values from gaussian centered at 0 with unitless
   * variance of 1
   *  2. Apply correction to set system net momentum to 0
   *  3. Apply correction to set system temperature to 120 K (or 1 in
   * characteristic units)
   */

  // 1. set velocities from gaussian distribution
  Vector3F NetMomentum = {0.0f, 0.0f, 0.0f};
  for (int i = 0; i < 108; ++i) {
    gV0[i].set_random();
    NetMomentum += gV0[i];
  }

  // 2. correct for momentum = 0 by subtracting the mean from each value and
  // compute the resulting temperature
  float Temp = 0.0f;
  for (int i = 0; i < 108; ++i) {
    gV0[i] -= (NetMomentum / 108.0f);
    Temp += gV0[i].sq();
  }
  Temp /= (3.0 * 107.0);

  // 3. set temperature of the system by multiplying by a correction factor
  float T_corr = 1.0 / sqrt(Temp);
  for (int i = 0; i < 108; ++i) {
    gV0[i] *= T_corr;
  }
}

void Sim::InitSim(SimArgs args) {
  fSimBox = args.SimulationCellBounds;
  fSimBoxHalf = 0.5f * fSimBox;
  LHalf = fSimBoxHalf[0];
  h = args.TimeStepSize;
  U_ij = args.PotentialEnergyFunction;
  F_ij = args.ForceFunction;
}

// Initialize position of atoms to FCC crystal lattice
void Sim::InitPositions() { InitFCC(fSimBox[0] / 3.0f); }

void Sim::InitVelocities() {
  InitV0();
  std::cout << "Starting to add particles to list...\n";
  fParticleHead = new Particle(FCC(0), V0(0));
  Particle *NewParticle = new Particle(FCC(1), V0(1));
  fParticleHead->next = NewParticle;
  for (int i = 2; i < fNumAtoms; ++i) {
    /*NewParticle->next = new Particle(FCC(i), V0(i));
    NewParticle = NewParticle->next;*/
    NewParticle = (NewParticle->next = new Particle(FCC(i), V0(i)));
    NewParticle->next = nullptr;
  }
  /*for (int i = 0; i < fNumAtoms; ++i) {
    fParticles->AddParticle(new Particle(FCC(i), V0(i)));
  }*/
}

Vector3F Sim::GetNetMomentum() const {
  Vector3F result = {0.0, 0.0, 0.0};
  // Loop over all particles in list
  Particle *entry;
  for (entry = fParticleHead; entry; entry = entry->next) {
    result += entry->v;
  }
  /*ParticleNode *entry;
  for (entry = this->fParticles->GetFirst(); entry; entry = entry->GetNext()) {
    result += entry->GetParticle()->v;
  }*/
  return result;
}

float Sim::GetTemperature() const {
  float result = 0.0;
  Particle *entry;
  for (entry = fParticleHead; entry; entry = entry->next) {
    result += entry->v.sq();
  }
  result /= 3.0 * (fNumAtoms - 1.0);
  /*ParticleNode *entry;
  for (entry = this->fParticles->GetFirst(); entry; entry = entry->GetNext()) {
    result += entry->GetParticle()->v.sq();
  }
  result /= 3.0 * (fNumAtoms - 1.0);*/
  return result;
}

void Sim::ApplyPeriodicBC() {
  Particle *entry;
  for (entry = fParticleHead; entry; entry = entry->next) {
    if (IsInSimBox(entry->r))
      continue;
    entry->r -= (copysign_vector3(fSimBoxHalf, entry->r - fSimBox) +
                 copysign_vector3(fSimBoxHalf, entry->r));
  }
  /*ParticleNode *entry;
  for (entry = this->fParticles->GetFirst(); entry; entry = entry->GetNext()) {
    if (IsInSimBox(entry->GetParticle()->r))
      continue;
    entry->GetParticle()->r -=
        (copysign_vector3(fSimBoxHalf, entry->GetParticle()->r - fSimBox) +
         copysign_vector3(fSimBoxHalf, entry->GetParticle()->r));
  }*/
  /*for (int i = 0; i < fNumAtoms; ++i) {
    if (IsInSimBox(fEnsemble.r(i)))
      continue;
    fEnsemble.r(i) -= (copysign_vector3(fSimBoxHalf, fEnsemble.r(i) - fSimBox) +
                       copysign_vector3(fSimBoxHalf, fEnsemble.r(i)));
  }*/
}

inline Vector3F Sim::copysign_vector3(const Vector3F &_lhs,
                                      const Vector3F &_rhs) {
  /*Vector3F result{0.0, 0.0, 0.0};
  result[0] = copysign(_lhs[0], _rhs[0]);
  result[1] = copysign(_lhs[1], _rhs[1]);
  result[2] = copysign(_lhs[2], _rhs[2]);
  return result;*/
  return Vector3F(copysign(_lhs[0], _rhs[0]), copysign(_lhs[1], _rhs[1]),
                  copysign(_lhs[2], _rhs[2]));
}

bool Sim::IsInSimBox(const Vector3F &check) {
  for (int i = 0; i < 3; ++i) {
    if (check[i] > fSimBox[i] || check[i] < 0)
      return false;
  }
  return true;
}

inline bool Sim::ImageIsCloser(const Vector3F &rij) {
  return (fabs(rij[0]) > LHalf || fabs(rij[1]) > LHalf || fabs(rij[2]) > LHalf);
}

// this function gets called a lot, so it's worth optimizing
const Vector3F &Sim::MinimumImage(const Vector3F &ri, const Vector3F &rj) {
  rij[0] = ri[0] - rj[0];
  rij[1] = ri[1] - rj[1];
  rij[2] = ri[2] - rj[2];
  if (ImageIsCloser(rij)) {
    rij[0] -= copysign(fSimBoxHalf[0], rij[0] - fSimBoxHalf[0]) +
              copysign(fSimBoxHalf[0], rij[0] + fSimBoxHalf[0]);
    rij[1] -= copysign(fSimBoxHalf[1], rij[1] - fSimBoxHalf[1]) +
              copysign(fSimBoxHalf[1], rij[1] + fSimBoxHalf[1]);
    rij[2] -= copysign(fSimBoxHalf[2], rij[2] - fSimBoxHalf[2]) +
              copysign(fSimBoxHalf[2], rij[2] + fSimBoxHalf[2]);
  }
  return rij;
}

void Sim::ComputeForces() {
  // zero them out
  Particle *entry;
  for (entry = fParticleHead; entry; entry = entry->next) {
    entry->a.assign(0.0);
  }
  entry = fParticleHead;
  while (entry->next) {
    Particle *entry2 = entry->next;
    while (entry2) {
      _force.assign_to(F_ij(MinimumImage(entry->r, entry2->r)));
      entry->a += _force;
      entry2->a -= _force;
      entry2 = entry2->next;
    }
    entry = entry->next;
  }
}

void Sim::VelocityVerletStep() {
  ComputeForces();
  Particle *entry = fParticleHead;
  while (entry) {
    entry->r += h * entry->v + 0.5f * h * h * entry->a;
    entry->v += 0.5f * h * entry->a;
    entry = entry->next;
  }
  ApplyPeriodicBC();
  ComputeForces();
  entry = fParticleHead;
  while (entry) {
    entry->v += 0.5f * h * entry->a;
    entry = entry->next;
  }
  /*ComputeForces();
  ParticleNode *entry;
  for (entry = this->fParticles->GetFirst(); entry; entry = entry->GetNext()) {
    entry->GetParticle()->r +=
        h * entry->GetParticle()->v + 0.5f * h * h * entry->GetParticle()->a;
    entry->GetParticle()->v += 0.5f * h * entry->GetParticle()->a;
  }
  ApplyPeriodicBC();
  ComputeForces();
  for (entry = this->fParticles->GetFirst(); entry; entry = entry->GetNext()) {
    entry->GetParticle()->v += 0.5f * h * entry->GetParticle()->a;
  }*/
  /*ComputeForces();
  for (int i = 0; i < fNumAtoms; ++i) {
    fEnsemble.r(i) += h * fEnsemble.v(i) + 0.5f * h * h * fEnsemble.a(i);
    fEnsemble.v(i) += 0.5f * h * fEnsemble.a(i);
  }
  ApplyPeriodicBC();
  ComputeForces();
  for (int i = 0; i < fNumAtoms; ++i) {
    fEnsemble.v(i) += 0.5f * h * fEnsemble.a(i);
  }*/
}

void Sim::PrintSample() {
  std::cout << "r[0] = " << this->fParticles->GetFirst()->GetParticle()->r
            << "v[0] = " << this->fParticles->GetFirst()->GetParticle()->v
            << "a[0] = " << this->fParticles->GetFirst()->GetParticle()->a
            << '\n';
  std::cout << "r[1] = "
            << this->fParticles->GetFirst()->GetNext()->GetParticle()->r
            << "v[1] = "
            << this->fParticles->GetFirst()->GetNext()->GetParticle()->v
            << "a[1] = "
            << this->fParticles->GetFirst()->GetNext()->GetParticle()->a
            << '\n';
}

float Sim::GetPotentialEnergy() {
  float result = 0.0f;
  Particle *entry;
  for (entry = fParticleHead; entry; entry = entry->next) {
    Particle *entry2;
    for (entry2 = entry->next; entry2; entry2 = entry2->next) {
      result += U_ij(MinimumImage(entry->r, entry2->r));
    }
  }
  return result;
  /*ParticleNode *entry;
  for (entry = this->fParticles->GetFirst(); entry; entry = entry->GetNext()) {
    ParticleNode *entry2;
    for (entry2 = entry->GetNext(); entry2; entry2 = entry2->GetNext()) {
      result +=
          U_ij(MinimumImage(entry->GetParticle()->r, entry2->GetParticle()->r));
    }
  }
  return result;*/
  /*float result{};
  for (int i = 0; i < fNumAtoms; ++i) {
    for (int j = i + 1; j < fNumAtoms; ++j) {
      result += U_ij(MinimumImage(fEnsemble.r(i), fEnsemble.r(j)));
    }
  }
  return result;*/
}

float Sim::GetKineticEnergy() {
  float result = 0.0f;
  Particle *entry;
  for (entry = fParticleHead; entry; entry = entry->next) {
    result += entry->v.sq();
  }
  return 0.5f * result;
  /*ParticleNode *entry;
  for (entry = this->fParticles->GetFirst(); entry; entry = entry->GetNext()) {
    result += entry->GetParticle()->v.sq();
  }
  return result;*/
  /*float result{};
  for (int i = 0; i < fNumAtoms; ++i) {
    result += fEnsemble.v(i).sq();
  }
  return result / 2.0;*/
}

float Sim::GetInternalEnergy() {
  return GetKineticEnergy() + GetPotentialEnergy();
}

void Sim::CorrectTemperature() {
  float T_corr = 1.0 / sqrt(GetTemperature());
  Particle *entry;
  for (entry = fParticleHead; entry; entry = entry->next) {
    entry->v *= T_corr;
  }
  /*ParticleNode *entry;
  for (entry = this->fParticles->GetFirst(); entry; entry = entry->GetNext()) {
    entry->GetParticle()->v *= T_corr;
  }*/
  /*for (int i = 0; i < fNumAtoms; ++i) {
    fEnsemble.v(i) *= T_corr;
  }*/
}

Vector3F Sim::GetMeanVelocity() const {
  Vector3F result = {0.0, 0.0, 0.0};
  Particle *entry;
  for (entry = fParticleHead; entry; entry = entry->next) {
    result += entry->v;
  }
  /*ParticleNode *entry;
  for (entry = this->fParticles->GetFirst(); entry; entry = entry->GetNext()) {
    result += entry->GetParticle()->v;
  }*/
  result /= 108.0f;
  return result;
  /*Vector3F result{};
  for (int i = 0; i < fNumAtoms; ++i) {
    result += fEnsemble.v(i);
  }
  result /= 108.0f;
  return result;*/
}

float Sim::Getfvelo(int j) const {
  float numerator = 0.0;
  float denominator = 0.0;
  Vector3F mean = GetMeanVelocity();
  Particle *entry;
  for (entry = fParticleHead; entry; entry = entry->next) {
    float tmpn = entry->v[j] - mean[j];
    numerator += tmpn * tmpn * tmpn * tmpn;
    denominator += tmpn * tmpn;
  }
  /*ParticleNode *entry;
  for (entry = this->fParticles->GetFirst(); entry; entry = entry->GetNext()) {
    float tmpn = entry->GetParticle()->v[j] - mean[j];
    numerator += tmpn * tmpn * tmpn * tmpn;
    denominator += tmpn * tmpn;
  }*/
  numerator /= 108.0f;
  denominator /= 108.0f;
  return numerator / (denominator * denominator);
  /*float numerator = 0.0;
  float denominator = 0.0;
  Vector3F mean = GetMeanVelocity();
  for (int i = 0; i < fNumAtoms; ++i) {
    float tmpn = fEnsemble.v(i)[j] - mean[j];
    numerator += tmpn * tmpn * tmpn * tmpn;
    denominator += tmpn * tmpn;
  }
  numerator /= 108.0f;
  denominator /= 108.0f;
  return numerator / (denominator * denominator);*/
}