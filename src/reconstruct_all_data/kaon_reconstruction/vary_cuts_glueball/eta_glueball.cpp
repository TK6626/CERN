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

// !./bashing/run_file.sh reconstruct_all_data/rho_reconstruction/vary_cuts_glueball/eta_glueball


/** make cuts on data to 
 * see how varying eta effects the (as so far of yet nonexistant peak at 2220Mev)
 */
int main() {

	// relevent constnats
	const Float_t pion_mass = 1e-3 * m_pi_char; // change to GeV)	

	// Load in aggregated uncut data
	ROOT::EnableImplicitMT();
	TString topology = "40";
	TString file = "data/glueball_mass_reconstruction/data_" + topology + "_SNR.root";
	RDF df_df("tree", file);
	
	RN df = df_df;
	
	Int_t index = 3;
	RVecF mass_bound_list = {125, 62, 31, 15, 10};
	RVecI num_bins = {250, 200, 200, 170, 150};
	RVecI lower_limit{1280, 1390, 1450, 1490, 1490};	
	Float_t mass_bound = mass_bound_list[index]; // MeV/c^2


	// eta has approx gaus distro near the mean=0, falls of sharply at eta =±3
    RVecF pt_cuts = {2.5, 2, 1.5, 1.2, 0.9};


	TH1F* hist = new TH1F("hist", "Invarient mass;M (MeV/c^{2});Events", num_bins[index], lower_limit[index], 2500);
	TCanvas* c = new TCanvas("c", "Canvas", 800, 600);
	TLegend* leg = new TLegend(0.7,0.4, 0.9,0.9, "selected cuts");
	auto stack = new THStack();
    	
	const char* trk_var = "trk_eta";

	Int_t i = 0;
	for (auto val : pt_cuts) {
    auto hist_cut = new TH1F(TString::Format("hist_cut_%d", i), "Invariant Mass;M (MeV/c^{2});Events", num_bins[index], lower_limit[index], 2500);

    auto df_filtered = cut_branch(df, trk_var, VectorConstraint<Float_t>{
        .vector_condition = [val](Float_t eta) { return TMath::Abs(eta) < val; },
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

	SaveOnlyCanvas(c, TString::Format("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/vary_cuts_glueball/eta_" + topology + "/%s mass_bound = %.4g MeV_c^2.root", trk_var, mass_bound), "RECREATE");

	return 0;
}
