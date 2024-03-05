{
  auto gE =
      new TGraph("../build/andersonthermostat.txt", "%lg %lg %*lg %*lg %*lg");
  gE->SetTitle("Energy vs Time;Time [r.u.];Energy [r.u.]");
  auto gKE =
      new TGraph("../build/andersonthermostat.txt", "%lg %*lg %lg %*lg %*lg");
  gKE->SetTitle("KE vs Time;Time [r.u.];KE [r.u.]");
  auto gPE =
      new TGraph("../build/andersonthermostat.txt", "%lg %*lg %*lg %lg %*lg");
  gPE->SetTitle("PE vs Time;Time [r.u.];PE [r.u.]");
  auto gTemp =
      new TGraph("../build/andersonthermostat.txt", "%lg %*lg %*lg %*lg %lg");
  gTemp->SetTitle("Temperature vs Time;Time [r.u.];Temp [r.u.]");

  auto c = new TCanvas();
  c->Divide(2, 2, 0, 0);
  c->cd(1);
  gE->Draw("AL");
  c->cd(2);
  gKE->Draw("AL");
  c->cd(3);
  gPE->Draw("AL");
  c->cd(4);
  gTemp->Draw("AL");
}