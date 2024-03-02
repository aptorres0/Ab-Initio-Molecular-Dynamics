
{
  auto gE = new TGraph("../build/CanoncialEnsembleII.txt",
                       "%lg %lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg");
  gE->SetTitle("Time Ave. Energy per Atom vs Time;Time [ps];Energy [eV]");
  auto gKE = new TGraph("../build/CanoncialEnsembleII.txt",
                        "%lg %*lg %lg %*lg %*lg %*lg %*lg %*lg %*lg");
  gKE->SetTitle(
      "Time Ave. Kinetic Energy per Atom vs Time;Time [ps];Energy [eV]");
  auto gU = new TGraph("../build/CanoncialEnsembleII.txt",
                       "%lg %*lg %*lg %lg %*lg %*lg %*lg %*lg %*lg");
  gU->SetTitle(
      "Time Ave. Potential Energy per Atom vs Time;Time [ps];Energy [eV]");
  auto gT = new TGraph("../build/CanoncialEnsembleII.txt",
                       "%lg %*lg %*lg %*lg %lg %*lg %*lg %*lg %*lg");
  gT->SetTitle("Time Ave. Temperature vs Time;Time [ps];Temperature [K]");
  auto gDT = new TGraph("../build/CanoncialEnsembleII.txt",
                        "%lg %*lg %*lg %*lg %*lg %lg %*lg %*lg %*lg");
  gDT->SetTitle(
      "Temperature Fluctuations vs Time;Time [ps];Temp. Fluctuations");
  auto gP = new TGraph("../build/CanoncialEnsembleII.txt",
                       "%lg %*lg %*lg %*lg %*lg %*lg %lg %*lg %*lg");
  gP->SetTitle("Time Ave. Pressure vs Time;Time [ps];Pressure [kPa]");
  auto gCv = new TGraph("../build/CanoncialEnsembleII.txt",
                        "%lg %*lg %*lg %*lg %*lg %*lg %*lg %lg %*lg");
  gCv->SetTitle(
      "Specific Heat Capacity vs Time;Time [ps];#mbox{C}_{V} [kJ/kg K]");
  auto gR2 = new TGraph("../build/CanoncialEnsembleII.txt",
                        "%lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %lg");
  gR2->SetTitle("Mean Square Displacement vs Time;Time [ps]; Mean Square "
                "Displacement [#mbox{nm}^{2}]");

  auto c = new TCanvas();
  c->Divide(1, 3);
  c->cd(1);
  gE->Draw();
  c->cd(2);
  gKE->Draw();
  c->cd(3);
  gU->Draw();

  auto c2 = new TCanvas();
  c2->Divide(1, 3);
  c2->cd(1);
  gT->Draw();
  c2->cd(2);
  gDT->Draw();
  c2->cd(3);
  gP->Draw();

  auto c3 = new TCanvas();
  c3->Divide(1, 2);
  c3->cd(1);
  gCv->Draw();
  c3->cd(2);
  gR2->Draw();
}