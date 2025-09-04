#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "THStack.h"

#include "../../../../lib/custom_definitions.h"
#include "../../../../lib/cut_branch.h"
#include "../../../../lib/computations.h"
#include "../../../../lib/apply_cuts.cpp"
#include "../../../../lib/plotting_params.h"
#include "../../../../lib/admin_utils.h"


// !./bashing/run_file.sh  reconstruct_all_data/rho_reconstruction/vary_cuts_glueball/dz_glueball
void fillHistFromP4(const RVecLor& p4, TH1* hist, Float_t m_rho, Float_t mass_bound);

/** make cuts on data to 
 * see how varying dz effects the (as so far of yet nonexistant peak at 2220Mev)
 */
int main() {

	// relevent constnats
	const Float_t pion_mass = 1e-3 * m_pi_char; // change to GeV)	

	// Load in aggregated uncut data
	ROOT::EnableImplicitMT();
	TString topology = "20";
	TString file = "data/glueball_mass_reconstruction/data_" + topology + "_SNR.root";
	RDF df_df("tree", file);
	
	RN df = df_df;
	
	Int_t index = 2;
	RVecF mass_bound_list = {125, 62, 31, 15, 10};
	//RVecI num_bins = {250, 200, 200, 170, 150};
	RVecI num_bins = {80, 50, 45, 40, 35};
	RVecI lower_limit{1280, 1390, 1450, 1490, 1490};	
	Float_t mass_bound = mass_bound_list[index]; // MeV/c^2


	// dz has approx gaus distro near the mean=0, falls of sharply at dz =±3
    RVecF cuts = {100, 15, 7, 4, 1.8};
	cuts /=1e2;

	TH1F* hist = new TH1F("hist", "Invarient mass;M (MeV/c^{2});Events", num_bins[index], lower_limit[index], 2500);
	TCanvas* c = new TCanvas("c", "Canvas", 800, 600);
	TLegend* leg = new TLegend(0.7,0.4, 0.9,0.9, "Selected Cuts (dm)");
	auto stack = new THStack();
    	
	const char* trk_var = "trk_dz";

	Int_t i = 0;
	for (auto val : cuts) {
    auto hist_cut = new TH1F(TString::Format("hist_cut_%d", i), "Invariant Mass;M (MeV/c^{2});Events", num_bins[index], lower_limit[index], 2500);

    auto df_filtered = cut_branch(df, trk_var, VectorConstraint<Float_t>{
        .vector_condition = [val](Float_t dz) { return TMath::Abs(dz) < val; },
        .op = Operator::all,
        .element_count = 0
    }, false);

    df_filtered.Foreach([&](const RVecLor& p4) {
        fillHistFromP4(p4, hist_cut, m_rho, mass_bound);
    }, {"four_momentum"});

    hist_cut->SetLineColor(kBlack);
    hist_cut->SetFillColor(color_scheme[i % color_scheme.size()]);
    leg->AddEntry(hist_cut, TString::Format("%s  < %.3g", trk_var, val), "f");
    stack->Add(hist_cut);
    ++i;
	}

	stack->Draw("HIST F nostack");
	stack->SetTitle("Invariant Mass;M (MeV/c^{2});" + TString::Format("Events [%.2g MeV/c^{2}]", stack->GetHistogram()->GetXaxis()->GetBinWidth(1)));
	leg->Draw();
	c->Update();

	SaveOnlyCanvas(c, TString::Format("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/vary_cuts_glueball/dz_" + topology + "/BIG BINS %s mass_bound = %.4g MeV_c^2.root", trk_var, mass_bound), "RECREATE");

	return 0;
}


void fillHistFromP4(const RVecLor& p4, TH1* hist, Float_t m_rho, Float_t mass_bound) {
    auto m0 = p4[0].M() * 1e3;
    auto m1 = p4[1].M() * 1e3;
    auto m2 = p4[2].M() * 1e3;
    auto m3 = p4[3].M() * 1e3;

    bool cond1 = TMath::Abs(m0 - m_rho) < mass_bound && TMath::Abs(m1 - m_rho) < mass_bound;
    bool cond2 = TMath::Abs(m2 - m_rho) < mass_bound && TMath::Abs(m3 - m_rho) < mass_bound;

    if (cond1) hist->Fill((p4[0] + p4[1]).M() * 1e3);
    if (cond2) hist->Fill((p4[2] + p4[3]).M() * 1e3);
}
