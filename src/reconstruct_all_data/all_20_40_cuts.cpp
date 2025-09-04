#include <vector>
#include <array>
#include <stdio.h>
#include <cmath>

#include "TMath.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TLegend.h"

// include headers containing useful functions and definitions
#include "../../lib/computations.h"
#include "../../lib/cut_branch.h"


// command
// !g++ src/reconstruct_all_data/all_20_40_cuts.cpp lib/*.cpp -Llib -ldict `root-config --cflags --libs` -o bin/reconstruct_all_data/all_20_40_cuts; ./bin/reconstruct_all_data/all_20_40_cuts

int main() {
    // Constants
    const Float_t pion_mass = 0.139; // GeV/c^2
    const Float_t x_min = 260;
    const Float_t x_max = 1000;
    const Int_t x_bins = 500;
    const Float_t y_min = x_min;
    const Float_t y_max = x_max;
    const Int_t y_bins = x_bins;

    // Load in all 20 and 40 combined type events, get column names, enable multithreading
    ROOT::EnableImplicitMT();

    std::vector<std::string> files_20 = {
        "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM20.root",
        "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM21.root",
        "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM22.root",
        "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM23.root"
    };

    std::vector<std::string> files_40 = {
        "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM40.root",
        "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM41.root",
        "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM42.root",
        "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM43.root"
    };

    RD df_20_df("tree", files_20);
    RD df_40_df("tree", files_40);

	std::cout << df_20_df.Describe() << std::endl;

    // Apply cuts and transformations
	
	// cut on track number 
    Int_t track_num = 4;
    VectorConstraint<Int_t> track_constraint{
        .vector_condition = [track_num](Int_t ntrk) { return ntrk == track_num; },
        .op = Operator::all
    };
    RN df_20 = cut_branch(df_20_df, "ntrk", track_constraint, true);
    RN df_40 = cut_branch(df_40_df, "ntrk", track_constraint, true);

    // cut on transverse momentum
	float_t low_pt_threshold = 100000;
    VectorConstraint<float_t> momentum_constraint{
        .vector_condition = [low_pt_threshold](float_t pt) { return pt < low_pt_threshold; },
        .op = Operator::all
    };
    df_20 = cut_branch(df_20, "alltrk_pt", momentum_constraint, true);
    df_40 = cut_branch(df_40, "alltrk_pt", momentum_constraint, true);

    // cut on track charge
	Int_t num_pos = 2, num_neg = 2;
    VectorConstraint<Int_t> charge_constraint_pos{
        .vector_condition = [](Int_t q) { return q == +1; },
        .op = Operator::equal,
        .element_count = num_pos
    };
    VectorConstraint<Int_t> charge_constraint_neg{
        .vector_condition = [](Int_t q) { return q == -1; },
        .op = Operator::equal,
        .element_count = num_neg
    };
    df_20 = cut_branch(df_20, "trk_q", charge_constraint_pos, true);
    df_20 = cut_branch(df_20, "trk_q", charge_constraint_neg, true);
    df_40 = cut_branch(df_40, "trk_q", charge_constraint_pos, true);
    df_40 = cut_branch(df_40, "trk_q", charge_constraint_neg, true);

    // cut on z interacion point
	Float_t ZPV_up = 5;
    Float_t ZPV_down = -5;
    VectorConstraint<Float_t> zpv_cut{
        .vector_condition = [ZPV_down, ZPV_up](Int_t z_pv) { return ZPV_down < z_pv && z_pv < ZPV_up; },
        .op = Operator::all
    };
    //df_20 = cut_branch(df_20, "zPV", zpv_cut, true);
    //df_40 = cut_branch(df_40, "zPV", zpv_cut, true);

	// cut on the distance to the z point
	Float_t dz_up = 1;
    Float_t dz_down = -1;
    VectorConstraint<Float_t> dz_cut{
        .vector_condition = [dz_down, dz_up](Int_t dz_val) { return dz_down < dz_val && dz_val < dz_up; },
        .op = Operator::all
    };
    //df_20 = cut_branch(df_20, "trk_dz", dz_cut, true);
    //df_40 = cut_branch(df_40, "trk_dz", dz_cut, true);

	// cut on the xy distance to the interaction point
	Float_t dxy_up = 0.5;
    Float_t dxy_down = -0.5;
    VectorConstraint<Float_t> dxy_cut{
        .vector_condition = [dxy_down, dxy_up](Int_t dxy_val) { return dxy_down < dxy_val && dxy_val < dxy_up; },
        .op = Operator::all
    };
    //df_20 = cut_branch(df_20, "trk_dxy", dxy_cut, true);
    //df_40 = cut_branch(df_40, "trk_dxy", dxy_cut, true);

    df_20 = df_20.Define("four_momentum",
        [track_num, pion_mass](const RVecI &trk_q, const RVecF &trk_pt, const RVecF &trk_eta, const RVecF &trk_phi) {
            std::array<Int_t, 2> pion_neg;
            std::array<Int_t, 2> pion_pos;
            Int_t pos_count = 0, neg_count = 0;
            ROOT::VecOps::RVec<TLorentzVector> pion_4p(track_num);

            for (Int_t i = 0; i < track_num; ++i) {
                if (trk_q[i] == +1 && pos_count < 2) {
                    pion_pos[pos_count++] = i;
                } else if (trk_q[i] == -1 && neg_count < 2) {
                    pion_neg[neg_count++] = i;
                }
            }

            if (pos_count < 2 || neg_count < 2) {
                return ROOT::VecOps::RVec<TLorentzVector>{};
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

            return ROOT::VecOps::RVec<TLorentzVector>{P_00, P_01, P_10, P_11};
        },
        {"trk_q", "trk_pt", "trk_eta", "trk_phi"}
    );

    df_40 = df_40.Define("four_momentum",
        [track_num, pion_mass](const RVecI &trk_q, const RVecF &trk_pt, const RVecF &trk_eta, const RVecF &trk_phi) {
            std::array<Int_t, 2> pion_neg;
            std::array<Int_t, 2> pion_pos;
            Int_t pos_count = 0, neg_count = 0;
            ROOT::VecOps::RVec<TLorentzVector> pion_4p(track_num);

            for (Int_t i = 0; i < track_num; ++i) {
                if (trk_q[i] == +1 && pos_count < 2) {
                    pion_pos[pos_count++] = i;
                } else if (trk_q[i] == -1 && neg_count < 2) {
                    pion_neg[neg_count++] = i;
                }
            }

            if (pos_count < 2 || neg_count < 2) {
                return ROOT::VecOps::RVec<TLorentzVector>{};
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

            return ROOT::VecOps::RVec<TLorentzVector>{P_00, P_01, P_10, P_11};
        },
        {"trk_q", "trk_pt", "trk_eta", "trk_phi"}
    );

    // Plot the invariant mass of the reconstructed particles
    TH2F* histo_20 = new TH2F("histo_20", "Invariant Mass;M_{1} (MeV/c^{2});M_{2} (MeV/c^{2})", 600, 300, 1200, 600, 300, 1200);
    df_20.Foreach(
        [&histo_20](const ROOT::VecOps::RVec<TLorentzVector>& four_momentum) {
            histo_20->Fill(four_momentum[0].M() * 1e3, four_momentum[1].M() * 1e3);
            histo_20->Fill(four_momentum[2].M() * 1e3, four_momentum[3].M() * 1e3);
        },
        {"four_momentum"}
    );

    TH2F* histo_40 = new TH2F("histo_40", "Invariant Mass;M_{1} (MeV/c^{2});M_{2} (MeV/c^{2})", 600, 300, 1200, 600, 300, 1200);
    df_40.Foreach(
        [&histo_40](const ROOT::VecOps::RVec<TLorentzVector>& four_momentum) {
            histo_40->Fill(four_momentum[0].M() * 1e3, four_momentum[1].M() * 1e3);
            histo_40->Fill(four_momentum[2].M() * 1e3, four_momentum[3].M() * 1e3);
        },
        {"four_momentum"}
    );

    TCanvas* can1 = new TCanvas("can1");
    histo_20->Draw("COLZ");
    histo_40->Draw("SAME");
    can1->Update();
    TFile* ouf = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/kaons_from_rho_y_proj_all_20_reconstructmass.root", "RECREATE");
    can1->Write();
    ouf->Close();

    // Save histograms as ROOT files
    TFile* f20 = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/histo_20.root", "RECREATE");
    histo_20->Write();
    f20->Close();

    TFile* f40 = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/histo_40.root", "RECREATE");
    histo_40->Write();
    f40->Close();

    // // Save the RDataFrames directly for further analysis
    df_20.Snapshot("tree1", "data/glueball_mass_reconstruction/data_20.root");
    df_40.Snapshot("tree1", "data/glueball_mass_reconstruction/data_40.root");

    return 0;
}
