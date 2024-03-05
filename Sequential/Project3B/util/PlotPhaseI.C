
{
  auto gE =
      new TGraph("../build/CanoncialEnsembleI.txt", "%lg %lg %*lg %*lg %*lg");
  gE->SetTitle("Total Energy of the System vs Time;Time [ps];Energy [eV]");
  auto gKE = new TGraph("../build/CanoncialEnsembleI.txt",
                        "%lg %*lg %lg %*lg %*lg %*lg");
  gKE->SetTitle("Total Kinetic Energy vs Time;Time [ps];Energy [eV]");
  auto gU = new TGraph("../build/CanoncialEnsembleI.txt",
                       "%lg %*lg %*lg %lg %*lg %*lg");
  gU->SetTitle("Total Potential Energy vs Time;Time [ps];Energy [eV]");
  auto gT = new TGraph("../build/CanoncialEnsembleI.txt",
                       "%lg %*lg %*lg %*lg %lg %*lg");
  gT->SetTitle("Temperature vs Time;Time [ps];Temperature [K]");
  auto gP = new TGraph("../build/CanoncialEnsembleI.txt",
                       "%lg %*lg %*lg %*lg %*lg %lg");
  gP->SetTitle("Pressure vs Time;Time [ps];Pressure [GPa]");

  /*auto cE = new TCanvas();
  gE->Draw();
  auto cKE = new TCanvas();
  gKE->Draw();
  auto cU = new TCanvas();
  gU->Draw();
  auto cT = new TCanvas();
  gT->Draw();
  auto cP = new TCanvas();
  gP->Draw();*/

  auto c = new TCanvas();
  c->Divide(1, 3);
  c->cd(1);
  gE->Draw();
  c->cd(2);
  gKE->Draw();
  c->cd(3);
  gU->Draw();

  auto c2 = new TCanvas();
  c2->Divide(1, 2);
  c2->cd(1);
  gT->Draw();
  c2->cd(2);
  gP->Draw();
}