{
  TTree *phase1_tree = new TTree("phase1_tree", "Phase I");
  phase1_tree->ReadFile("../data/CanonicalII.txt", "t/D:E:KE:PE:T:dT:P:Cv:MSD",
                        ' ');

  auto h_Et = new TH2D("h_Et", "Energy vs Time", 100, 0, 100, 100, 0, 100);
  phase1_tree->Draw("E:t");
}