{
  auto gE = new TGraph("../data/Step7.txt", "%lg %lg %*lg %*lg");
  gE->SetTitle("System Interal Energy vs Time;Time [ps];Energy [eV]");
  auto gKE = new TGraph("../data/Step7.txt", "%lg %*lg %lg %*lg");
  gKE->SetTitle("Total KE vs Time;Time [ps];KE [eV]");
  auto gPE = new TGraph("../data/Step7.txt", "%lg %*lg %*lg %lg");
  gPE->SetTitle("Total PE vs Time;Time [ps];PE [eV]");

  auto c = new TCanvas("Step 7", "Step 7");
  c->Divide(1, 3);
  c->cd(1);
  gE->Draw("AL");
  c->cd(2);
  gKE->Draw("AL");
  c->cd(3);
  gPE->Draw("AL");
}