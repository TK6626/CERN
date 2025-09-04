ROOT::RDF::RNode applyCuts(ROOT::RDF::RNode df, Int_t track_num, Float_t low_pt_threshold, Int_t trk_charge) {
    auto df_cut1 = cut_branch_single(df, track_num, "ntrk", "==", true);
    auto df_cut2 = cut_branch_single(df_cut1, low_pt_threshold, "alltrk_pt", "<", true);
    auto df_cut3 = cut_branch_sum(df_cut2, trk_charge, "trk_q", "==", true);
    return df_cut3;
}


ROOT::RDF::RNode defineFourMomentum(ROOT::RDF::RNode df, Int_t track_num, Float_t pion_mass) {
    return df.Define("four_momentum",
        [track_num, pion_mass](const RVecI &trk_q, const RVecF &trk_pt, const RVecF &trk_eta, const RVecF &trk_phi) {
            std::array<Int_t, 2> pion_neg;
            std::array<Int_t, 2> pion_pos;
            Int_t pos_count = 0, neg_count = 0;
            std::array<TLorentzVector, 4> pion_4p;

            for (Int_t i = 0; i < track_num; ++i) {
                if (trk_q[i] == +1 && pos_count < 2) {
                    pion_pos[pos_count++] = i;
                } else if (trk_q[i] == -1 && neg_count < 2) {
                    pion_neg[neg_count++] = i;
                }
            }

            if (pos_count < 2 || neg_count < 2) {
                return std::array<TLorentzVector, 4>{};
            }

            for (Int_t i = 0; i < track_num; ++i) {
                TLorentzVector p_pion;
                p_pion.SetPtEtaPhiM(trk_pt[i], trk_eta[i], trk_phi[i], pion_mass);
                pion_4p[i] = p_pion;
            }

            TLorentzVector P_00 = pion_4p[pion_pos[0]] + pion_4p[pion_neg[0]];
            TLorentzVector P_11 = pion_4p[pion_pos[1]] + pion_4p[pion_neg[1]];
            TLorentzVector P_01 = pion_4p[pion_pos[0]] + pion_4p[pion_neg[1]];
            TLorentzVector P_10 = pion_4p[pion_pos[1]] + pion_4p[pion_neg[0]];

            return std::array<TLorentzVector, 4>{P_00, P_01, P_10, P_11};
        },
        {"trk_q", "trk_pt", "trk_eta", "trk_phi"});
}



int main() {
    // Constants
    const Float_t pion_mass = 0.139; // GeV/c^2
    const Int_t track_num = 4;
    const Float_t low_pt_threshold = 100000;
    const Int_t trk_charge = 0;

    // Load datasets
    std::vector<std::string> files_20 = {"/path/to/TOTEM20.root", "/path/to/TOTEM21.root"};
    std::vector<std::string> files_40 = {"/path/to/TOTEM40.root", "/path/to/TOTEM41.root"};
    RD df_20("tree", files_20);
    RD df_40("tree", files_40);

    // Apply cuts
    auto df_20_cut3 = applyCuts(df_20, track_num, low_pt_threshold, trk_charge);
    auto df_40_cut3 = applyCuts(df_40, track_num, low_pt_threshold, trk_charge);

    // Define four-momentum
    auto df_20_cut4 = defineFourMomentum(df_20_cut3, track_num, pion_mass);
    auto df_40_cut4 = defineFourMomentum(df_40_cut3, track_num, pion_mass);

    // Create histograms
    TH2F* histo_20 = new TH2F("histo_20", "Invariant Mass;M_{1} (MeV/c^{2});M_{2} (MeV/c^{2})", 300, 260, 1000, 300, 260, 1000);
    TH2F* histo_40 = new TH2F("histo_40", "Invariant Mass;M_{1} (MeV/c^{2});M_{2} (MeV/c^{2})", 300, 260, 1000, 300, 260, 1000);

    // Fill histograms
    fillHistogram(df_20_cut4, histo_20);
    fillHistogram(df_40_cut4, histo_40);

    // Create projections
    TH1D* y_projection_20 = createProjection(histo_20, "y_projection_20", 484, 504);
    TH1D* y_projection_40 = createProjection(histo_40, "y_projection_40", 484, 504);

    // Fit projections
    auto r_20 = fitProjection(y_projection_20, "fit_y_20", 480, 515);
    auto r_40 = fitProjection(y_projection_40, "fit_y_40", 480, 515);

    // Draw and save results
    TCanvas* canvas = new TCanvas("canvas", "Fitted Mass Peaks", 800, 600);
    y_projection_20->Draw("E1");
    y_projection_40->Draw("E1 SAME");
    canvas->SaveAs("output.pdf");

    return 0;
}