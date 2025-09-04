#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "THStack.h"
#include "TLine.h"
#include "TLegend.h"



#include "../../../../lib/custom_definitions.h"
#include "../../../../lib/cut_branch.h"
#include "../../../../lib/computations.h"
#include "../../../../lib/apply_cuts.cpp"
#include "../../../../lib/plotting_params.h"
#include "../../../../lib/admin_utils.h"

//  !./bashing/run_file.sh reconstruct_all_data/rho_reconstruction/vary_cuts_glueball/vary_rho_mass


/** make cuts on data to 
 * see how varying dxy effects the (as so far of yet nonexistant peak at 2220Mev)
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
	
	Int_t index = 1;
	RVecF mass_bound_list = {125, 62, 31, 15, 10};
	RVecI num_bins = {170, 130, 120, 100, 85};
	RVecI lower_limit{1220, 1270, 1300, 1320, 1330};	
	Float_t mass_bound = mass_bound_list[index]; // MeV/c^2

	RVecI m_rho_list = {850, 820, 800, 780, 770, 760};  // MeV
	TH1F* hist = new TH1F("hist", "Invarient mass;M (MeV/c^{2});Events", num_bins[index], lower_limit[index], 2500);
	TCanvas* c = new TCanvas("c", "Canvas", 800, 600);
	TLegend* leg = new TLegend(0.7,0.4, 0.9,0.9, "selected cuts");
	auto stack = new THStack();
    	

	Int_t i = 0;
	for (auto val : m_rho_list) {
    auto hist_cut = new TH1F(TString::Format("hist_cut_%d", i), "Invariant Mass;M (MeV/c^{2});Events", num_bins[index], lower_limit[index], 2500);

    df.Foreach([&](const RVecLor& p4) {
        fillHistFromP4(p4, hist_cut, val, mass_bound);
    }, {"four_momentum"});

    hist_cut->SetLineColor(kBlack);
    hist_cut->SetFillColor(color_scheme[i % color_scheme.size()]);
    leg->AddEntry(hist_cut, TString::Format("m = %d", val), "f");
    stack->Add(hist_cut);
    ++i;
	}

	stack->Draw("HIST F nostack");
	stack->SetTitle("Invariant Mass;M (MeV/c^{2});" + TString::Format("Events [%.2g MeV/c^{2}]", stack->GetHistogram()->GetXaxis()->GetBinWidth(1)));
	TLine* line = DrawVLine(2220, stack);
	line->Draw();
	leg->Draw();
	c->Update();

	SaveOnlyCanvas(c, TString::Format("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/vary_cuts_glueball/rho_mass_bound_" + topology + "/mass_bound = %.4g MeV_c^2.root", mass_bound), "RECREATE");

	return 0;
}


