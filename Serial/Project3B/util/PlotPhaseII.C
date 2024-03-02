
{
  auto gE = new TGraph("../build/CanoncialEnsembleII.txt",
                       "%lg %lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg");
  gE->SetTitle("Total Energy of the System vs Time;Time [ps];Energy [eV]");
  auto gKE = new TGraph("../build/CanoncialEnsembleII.txt",
                        "%lg %*lg %lg %*lg %*lg %*lg %*lg %*lg %*lg");
  gKE->SetTitle("Total Kinetic Energy vs Time;Time [ps];Energy [eV]");
  auto gU = new TGraph("../build/CanoncialEnsembleII.txt",
                       "%lg %*lg %*lg %lg %*lg %*lg %*lg %*lg %*lg");
  gU->SetTitle("Total Potential Energy vs Time Step;Time [ps];Energy [eV]");
  auto gT = new TGraph("../build/CanoncialEnsembleII.txt",
                       "%lg %*lg %*lg %*lg %lg %*lg %*lg %*lg %*lg");
  gT->SetTitle("Temperature vs Time Step;Time [ps];Temperature [K]");
  auto gDT = new TGraph("../build/CanoncialEnsembleII.txt",
                        "%lg %*lg %*lg %*lg %*lg %lg %*lg %*lg %*lg");
  gDT->SetTitle(
      "Temperature Fluctuations vs Time Step;Time [ps];Temp. Fluctuations");
  auto gP = new TGraph("../build/CanoncialEnsembleII.txt",
                       "%lg %*lg %*lg %*lg %*lg %*lg %lg %*lg %*lg");
  gP->SetTitle("Pressure vs Time Step;Time [ps];Pressure [GPa]");
  auto gCv = new TGraph("../build/CanoncialEnsembleII.txt",
                        "%lg %*lg %*lg %*lg %*lg %*lg %*lg %lg %*lg");
  gCv->SetTitle("Heat Capacity vs Time Step;Time [ps];Heat Capacity [r.u.]");
  auto gR2 = new TGraph("../build/CanoncialEnsembleII.txt",
                        "%lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %lg");
  gR2->SetTitle("Mean Square Displacement vs Time Step;Time [ps]; Mean Square "
                "Displacement [#mbox{nm}^{2}]");

  /*auto func = new TF1("func", "6*x*[0]");
  func->SetParName(0, "#D");
  gR2->Fit(func);
  gR2->PaintStats(func);*/

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

  /*auto cE = new TCanvas();
  gE->Draw();
  auto cKE = new TCanvas();
  gKE->Draw();
  auto cU = new TCanvas();
  gU->Draw();
  auto cT = new TCanvas();
  gT->Draw();
  auto cDT = new TCanvas();
  gDT->Draw();
  auto cP = new TCanvas();
  gP->Draw();
  auto cCv = new TCanvas();
  gCv->Draw();
  auto cR2 = new TCanvas();
  gR2->Draw();*/
}