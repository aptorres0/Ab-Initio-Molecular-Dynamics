{
  auto canvas = new TCanvas("canvas");
  auto mg = new TMultiGraph("mg", "mg");

  auto gE = new TGraph("../build/energy.txt", "%lg %lg %*lg %*lg");
  gE->SetTitle("E");
  auto gT = new TGraph("../build/energy.txt", "%lg %*lg %lg %*lg");
  gT->SetTitle("T");
  auto gU = new TGraph("../build/energy.txt", "%lg %*lg %*lg %lg");
  gU->SetTitle("U");

  canvas->Divide(2, 2, 0, 0);
  canvas->cd(1);
  gE->Draw();
  canvas->cd(2);
  gT->Draw();
  canvas->cd(3);
  gU->Draw();
}