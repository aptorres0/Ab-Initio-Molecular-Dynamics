
{
  auto gE = new TGraph("../build/energies.txt", "%lg %lg %*lg %*lg");
  gE->SetTitle("Total Energy of the System vs Time Step;Time Step");
  gE->GetXaxis()->SetTitle("Time Step");
  gE->GetXaxis()->CenterTitle(true);
  auto gT = new TGraph("../build/energies.txt", "%lg %*lg %lg %*lg");
  gT->SetTitle("Total Kinetic Energy vs Time Step;Time Step");
  gT->GetXaxis()->SetTitle("Time Step");
  gT->GetXaxis()->CenterTitle(true);
  auto gU = new TGraph("../build/energies.txt", "%lg %*lg %*lg %lg");
  gU->SetTitle("Total Potential Energy vs Time Step;Time Step");
  gU->GetXaxis()->SetTitle("Time Step");
  gU->GetXaxis()->CenterTitle(true);

  auto c = new TCanvas();
  c->Divide(2, 2, 0, 0);
  c->cd(1);
  gE->Draw();
  c->cd(2);
  gT->Draw();
  c->cd(3);
  gU->Draw();

  auto cE = new TCanvas();
  gE->Draw("AL");
}