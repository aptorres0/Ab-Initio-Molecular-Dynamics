
{
  auto gE = new TGraph("../build/energies.txt", "%lg %lg %*lg %*lg");
  gE->SetTitle("E;Time Step");
  auto gT = new TGraph("../build/energies.txt", "%lg %*lg %lg %*lg");
  gT->SetTitle("T;Time Step");
  auto gU = new TGraph("../build/energies.txt", "%lg %*lg %*lg %lg");
  gU->SetTitle("U;Time Step");

  auto c = new TCanvas();
  c->Divide(2, 2, 0, 0);
  c->cd(1);
  gE->Draw();
  c->cd(2);
  gT->Draw();
  c->cd(3);
  gU->Draw();
}