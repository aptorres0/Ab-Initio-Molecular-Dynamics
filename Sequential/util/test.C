{
  // Step 7
  auto gE_7 = new TGraph("/Users/ajpt/School/Project3/data/Step7_new.txt",
                         "%lg %lg %*lg %*lg");
  gE_7->SetTitle("Energy vs Time;Time [ps];Energy [eV]");
  auto gKE_7 = new TGraph("/Users/ajpt/School/Project3/data/Step7_new.txt",
                          "%lg %*lg %lg %*lg");
  gKE_7->SetTitle("KE vs Time;Time [ps];KE [eV]");
  auto gPE_7 = new TGraph("/Users/ajpt/School/Project3/data/Step7_new.txt",
                          "%lg %*lg %*lg %lg");
  gPE_7->SetTitle("PE vs Time;Time [ps];PE [eV]");

  auto c7 = new TCanvas("step7", "step7");
  c7->Divide(1, 3);
  c7->cd(1);
  gE_7->Draw("AL");
  c7->cd(2);
  gKE_7->Draw("AL");
  c7->cd(3);
  gPE_7->Draw("AL");

  // Microcanonical Ensemble
  auto gE_m =
      new TGraph("/Users/ajpt/School/Project3/data/Microcanonical_new.txt",
                 "%lg %lg %*lg %*lg");
  gE_m->SetTitle("Energy vs Time;Time [ps];Energy [eV]");
  auto gKE_m =
      new TGraph("/Users/ajpt/School/Project3/data/Microcanonical_new.txt",
                 "%lg %*lg %lg %*lg");
  gKE_m->SetTitle("KE vs Time;Time [ps];KE [eV]");
  auto gPE_m =
      new TGraph("/Users/ajpt/School/Project3/data/Microcanonical_new.txt",
                 "%lg %*lg %*lg %lg");
  gPE_m->SetTitle("PE vs Time;Time [ps];PE [eV]");

  auto cm = new TCanvas("micro", "micro");
  cm->Divide(1, 3);
  cm->cd(1);
  gE_m->Draw("AL");
  cm->cd(2);
  gKE_m->Draw("AL");
  cm->cd(3);
  gPE_m->Draw("AL");

  // Canonical Ensemble - ALL
  auto gE_all =
      new TGraph("/Users/ajpt/School/Project3/data/Canonical_All_new.txt",
                 "%lg %lg %*lg %*lg %*lg %*lg");
  gE_all->SetTitle("Energy vs Time;Time [ps];Energy [eV]");
  auto gKE_all =
      new TGraph("/Users/ajpt/School/Project3/data/Canonical_All_new.txt",
                 "%lg %*lg %lg %*lg %*lg %*lg");
  gKE_all->SetTitle("KE vs Time;Time [ps];KE [eV]");
  auto gPE_all =
      new TGraph("/Users/ajpt/School/Project3/data/Canonical_All_new.txt",
                 "%lg %*lg %*lg %lg %*lg %*lg");
  gPE_all->SetTitle("PE vs Time;Time [ps];PE [eV]");
  auto gT_all =
      new TGraph("/Users/ajpt/School/Project3/data/Canonical_All_new.txt",
                 "%lg %*lg %*lg %*lg %lg %*lg");
  gT_all->SetTitle("Temperature vs Time;Time [ps];Temperature [K]");
  auto gP_all =
      new TGraph("/Users/ajpt/School/Project3/data/Canonical_All_new.txt",
                 "%lg %*lg %*lg %*lg %*lg %lg");
  gP_all->SetTitle("Pressure vs Time;Time [ps];Pressure [MPa]");

  auto call_E = new TCanvas("Canonical - All E", "Canonical - All E");
  call_E->Divide(1, 3);
  call_E->cd(1);
  gE_all->Draw("AL");
  call_E->cd(2);
  gKE_all->Draw("AL");
  call_E->cd(3);
  gPE_all->Draw("AL");
  auto call_PT = new TCanvas("Canonical - All PT", "Canonical - All PT");
  call_PT->Divide(1, 2);
  call_PT->cd(1);
  gT_all->Draw("AL");
  call_PT->cd(2);
  gP_all->Draw("AL");

  // Canonical Ensemble - Phase I
  auto gE_M = new TGraph("/Users/ajpt/School/Project3/data/CanonicalI_new.txt",
                         "%lg %lg %*lg %*lg");
  gE_M->SetTitle("Energy vs Time;Time [ps];Energy [eV]");
  auto gKE_M = new TGraph("/Users/ajpt/School/Project3/data/CanonicalI_new.txt",
                          "%lg %*lg %lg %*lg");
  gKE_M->SetTitle("KE vs Time;Time [ps];KE [eV]");
  auto gPE_M = new TGraph("/Users/ajpt/School/Project3/data/CanonicalI_new.txt",
                          "%lg %*lg %*lg %lg");
  gPE_M->SetTitle("PE vs Time;Time [ps];PE [eV]");
  auto gT_M = new TGraph("/Users/ajpt/School/Project3/data/CanonicalI_new.txt",
                         "%lg %*lg %*lg %*lg %lg");
  gT_M->SetTitle("Temperature vs Time;Time [ps];Temperature [K]");
  auto gP_M = new TGraph("/Users/ajpt/School/Project3/data/CanonicalI_new.txt",
                         "%lg %*lg %*lg %*lg %*lg %lg");
  gP_M->SetTitle("Pressure vs Time;Time [ps];Pressure [MPa]");

  auto cM_E = new TCanvas("MacroE", "MacroE");
  cM_E->Divide(1, 3);
  cM_E->cd(1);
  gE_M->Draw("AL");
  cM_E->cd(2);
  gKE_M->Draw("AL");
  cM_E->cd(3);
  gPE_M->Draw("AL");
  auto cM_PT = new TCanvas("MacroPT", "MacroPT");
  cM_PT->Divide(1, 2);
  cM_PT->cd(1);
  gT_M->Draw("AL");
  cM_PT->cd(2);
  gP_M->Draw("AL");

  // Canonical Ensemble - Phase II
  auto gE_MII =
      new TGraph("/Users/ajpt/School/Project3/data/CanonicalII_new.txt",
                 "%lg %lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg");
  gE_MII->SetTitle("Time Ave. (<KE> + <PE>)/N vs Time;Time [ps];Energy [eV]");
  auto gKE_MII =
      new TGraph("/Users/ajpt/School/Project3/data/CanonicalII_new.txt",
                 "%lg %*lg %lg %*lg %*lg %*lg %*lg %*lg %*lg");
  gKE_MII->SetTitle("Time Ave. <KE>/N vs Time;Time [ps];KE [eV]");
  auto gPE_MII =
      new TGraph("/Users/ajpt/School/Project3/data/CanonicalII_new.txt",
                 "%lg %*lg %*lg %lg %*lg %*lg %*lg %*lg %*lg");
  gPE_MII->SetTitle("Time Ave. <PE>/N vs Time;Time [ps];PE [eV]");
  auto gT_MII =
      new TGraph("/Users/ajpt/School/Project3/data/CanonicalII_new.txt",
                 "%lg %*lg %*lg %*lg %lg %*lg %*lg %*lg %*lg");
  gT_MII->SetTitle("Temperature vs Time;Time [ps];Temperature [K]");
  auto gdT_MII =
      new TGraph("/Users/ajpt/School/Project3/data/CanonicalII_new.txt",
                 "%lg %*lg %*lg %*lg %*lg %lg %*lg %*lg %*lg");
  gdT_MII->SetTitle("Temperature Fluctuations vs Time;Time [ps]; #delta T [1]");
  auto gP_MII =
      new TGraph("/Users/ajpt/School/Project3/data/CanonicalII_new.txt",
                 "%lg %*lg %*lg %*lg %*lg %*lg %lg %*lg %*lg");
  gP_MII->SetTitle("Pressure vs Time;Time [ps];Pressure [MPa]");

  auto gCv_II = new TGraph("../data/CanonicalII_new.txt",
                           "%lg %*lg %*lg %*lg %*lg %*lg %*lg %lg %*lg");
  gCv_II->SetTitle("Heat Capacity vs Time; Time [ps]; #mbox{C}_{v} [J/kg/K]");
  auto gMSD_II = new TGraph("../data/CanonicalII_new.txt",
                            "%lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %lg");
  gMSD_II->SetTitle(
      "Mean Square Displacement vs Time; Time [ps]; MSD [#mbox{nm}^{2}]");
  auto gRadialDist_II = new TGraph("../data/radialdist_new.txt", "%lg %lg");
  gRadialDist_II->SetTitle("Radial Distribution Function; r [nm]; g(r)");

  auto cMII_E = new TCanvas("MacroE - II", "MacroE - II");
  cMII_E->Divide(1, 3);
  cMII_E->cd(1);
  gE_MII->Draw("AL");
  cMII_E->cd(2);
  gKE_MII->Draw("AL");
  cMII_E->cd(3);
  gPE_MII->Draw("AL");

  auto cMII_PT = new TCanvas("MacroPT - II", "MacroPT - II");
  cMII_PT->Divide(1, 3);
  cMII_PT->cd(1);
  gT_MII->Draw("AL");
  cMII_PT->cd(2);
  gdT_MII->Draw("AL");
  cMII_PT->cd(3);
  gP_M->Draw("AL");

  auto c_Cv = new TCanvas("Cv", "Cv");
  gCv_II->Draw("AL");

  auto c_MSD = new TCanvas("MSD", "MSD");
  gMSD_II->Draw("AL");

  auto c_RadialDist = new TCanvas("RadialDist", "RadialDist");
  gRadialDist_II->Draw("AL");
}