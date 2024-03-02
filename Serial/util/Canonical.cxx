{
  // Phase I
  auto gE_I =
      new TGraph("../data/CanonicalI.txt", "%lg %lg %*lg %*lg %*lg %*lg");
  gE_I->SetTitle("System Interal Energy vs Time;Time [ps];Energy [eV]");
  auto gKE_I =
      new TGraph("../data/CanonicalI.txt", "%lg %*lg %lg %*lg %*lg %*lg");
  gKE_I->SetTitle("Total KE vs Time;Time [ps];KE [eV]");
  auto gPE_I =
      new TGraph("../data/CanonicalI.txt", "%lg %*lg %*lg %lg %*lg %*lg");
  gPE_I->SetTitle("Total PE vs Time;Time [ps];PE [eV]");
  auto gTemp_I =
      new TGraph("../data/CanonicalI.txt", "%lg %*lg %*lg %*lg %lg %*lg");
  gTemp_I->SetTitle("Temperature vs Time;Time [ps];Temp [K]");
  auto gPress_I =
      new TGraph("../data/CanonicalI.txt", "%lg %*lg %*lg %*lg %*lg %lg");
  gPress_I->SetTitle("Pressure vs Time; Time [ps]; Pressure [MPa]");

  auto c_Ia = new TCanvas("Phase I Energy", "Phase I Energy");
  c_Ia->Divide(1, 3);
  c_Ia->cd(1);
  gE_I->Draw("AL");
  c_Ia->cd(2);
  gKE_I->Draw("AL");
  c_Ia->cd(3);
  gPE_I->Draw("AL");
  auto c_Ib = new TCanvas("Phase I Temp & Press", "Phase I Temp & Press");
  c_Ib->Divide(1, 2);
  c_Ib->cd(1);
  gTemp_I->Draw("AL");
  c_Ib->cd(2);
  gPress_I->Draw("AL");

  // Phase II
  auto gE_II = new TGraph("../data/CanonicalII.txt",
                          "%lg %lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg");
  gE_II->SetTitle("Time Ave. (<KE> + <PE>)/N vs Time;Time [ps];Energy [eV]");
  auto gKE_II = new TGraph("../data/CanonicalII.txt",
                           "%lg %*lg %lg %*lg %*lg %*lg %*lg %*lg %*lg");
  gKE_II->SetTitle("Time Ave <KE>/N vs Time;Time [ps];KE [eV]");
  auto gPE_II = new TGraph("../data/CanonicalII.txt",
                           "%lg %*lg %*lg %lg %*lg %*lg %*lg %*lg %*lg");
  gPE_II->SetTitle("Time Ave. <PE>/N vs Time;Time [ps];PE [eV]");
  auto gTemp_II = new TGraph("../data/CanonicalII.txt",
                             "%lg %*lg %*lg %*lg %lg %*lg %*lg %*lg %*lg");
  gTemp_II->SetTitle("Time Ave. <T> vs Time;Time [ps];Temp [K]");
  auto gDT_II = new TGraph("../data/CanonicalII.txt",
                           "%lg %*lg %*lg %*lg %*lg %lg %*lg %*lg %*lg");
  gDT_II->SetTitle("Temperature Fluctuation vs Time;Time [ps];#delta T [K]");
  auto gPress_II = new TGraph("../data/CanonicalII.txt",
                              "%lg %*lg %*lg %*lg %*lg %*lg %lg %*lg %*lg");
  gPress_II->SetTitle("Time Ave. <P> vs Time; Time [ps]; Pressure [MPa]");
  auto gCv_II = new TGraph("../data/CanonicalII.txt",
                           "%lg %*lg %*lg %*lg %*lg %*lg %*lg %lg %*lg");
  gCv_II->SetTitle("Heat Capacity vs Time; Time [ps]; Cv [J/kg/K]");
  auto gMSD_II = new TGraph("../data/CanonicalII.txt",
                            "%lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %lg");
  gMSD_II->SetTitle("Mean Square Displacement vs Time; Time [ps]; MSD [nm^2]");
  auto gRadialDist_II = new TGraph("../data/radialdist.txt", "%lg %lg");
  gRadialDist_II->SetTitle("Radial Distribution Function; r [nm]; g(r)");

  auto c_IIa = new TCanvas("Phase IIa", "Phase IIa");
  c_IIa->Divide(1, 3);
  c_IIa->cd(1);
  gE_II->Draw("AL");
  c_IIa->cd(2);
  gKE_II->Draw("AL");
  c_IIa->cd(3);
  gPE_II->Draw("AL");

  auto c_IIb = new TCanvas("Phase IIb", "Phase IIb");
  c_IIb->Divide(1, 3);
  c_IIb->cd(1);
  gTemp_II->Draw("AL");
  c_IIb->cd(2);
  gDT_II->Draw("AL");
  c_IIb->cd(3);
  gPress_II->Draw("AL");

  auto c_Cv = new TCanvas("Cv", "Cv");
  gCv_II->Draw("AL");

  auto c_MSD = new TCanvas("MSD", "MSD");
  gMSD_II->Draw("AL");

  auto c_RadialDist = new TCanvas("RadialDist", "RadialDist");
  gRadialDist_II->Draw("AL");
}