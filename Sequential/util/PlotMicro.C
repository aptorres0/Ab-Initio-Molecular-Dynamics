
{
  auto gE =
      new TGraph("../build/MicrocanonicalEnsembleSim.txt", "%lg %lg %*lg %*lg");
  gE->SetTitle("Total Energy of the System vs Time;Time [ps];Energy [eV]");
  auto gKE =
      new TGraph("../build/MicrocanonicalEnsembleSim.txt", "%lg %*lg %lg %*lg");
  gKE->SetTitle("Total Kinetic Energy vs Time;Time [ps];Energy [eV]");
  auto gU =
      new TGraph("../build/MicrocanonicalEnsembleSim.txt", "%lg %*lg %*lg %lg");
  gU->SetTitle("Total Potential Energy vs Time;Time [ps];Energy [eV]");

  /*auto cE = new TCanvas();
  gE->Draw();
  auto cKE = new TCanvas();
  gKE->Draw();
  auto cU = new TCanvas();
  gU->Draw();*/

  auto c = new TCanvas();
  c->Divide(1, 3);
  c->cd(1);
  gE->Draw();
  c->cd(2);
  gKE->Draw();
  c->cd(3);
  gU->Draw();
}