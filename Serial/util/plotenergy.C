
{
  auto gE = new TGraph("../build/CanoncialEnsembleI.txt",
                       "%lg %lg %*lg %*lg %*lg %*lg");
  gE->SetTitle("Total Energy of the System vs Time Step;Time Step");
  auto gKE = new TGraph("../build/CanoncialEnsembleI.txt",
                        "%lg %*lg %lg %*lg %*lg %*lg");
  gKE->SetTitle("Total Kinetic Energy vs Time Step;Time Step");
  auto gU = new TGraph("../build/CanoncialEnsembleI.txt",
                       "%lg %*lg %*lg %lg %*lg %*lg");
  gU->SetTitle("Total Potential Energy vs Time Step;Time Step");
  auto gT = new TGraph("../build/CanoncialEnsembleI.txt",
                       "%lg %*lg %*lg %*lg %lg %*lg");
  gT->SetTitle("Temperature vs Time Step;Time Step;Temperature");
  auto gP = new TGraph("../build/CanoncialEnsembleI.txt",
                       "%lg %*lg %*lg %*lg %*lg %lg");
  gP->SetTitle("Pressure vs Time Step;Time Step;Pressure");

  auto c = new TCanvas();
  c->Divide(3, 2, 0, 0);
  c->cd(1);
  gE->Draw();
  c->cd(2);
  gKE->Draw();
  c->cd(3);
  gU->Draw();
  c->cd(4);
  gT->Draw();
  c->cd(5);
  gP->Draw();
}