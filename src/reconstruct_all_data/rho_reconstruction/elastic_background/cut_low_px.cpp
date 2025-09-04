#include <stdio.h>
#include <stdlib.h>

#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "THStack.h"
#include "TLorentzVector.h"
#include "TLegend.h"

#include "../../../../lib/custom_definitions.h"
#include "../../../../lib/cut_branch.h"
#include "../../../../lib/computations.h"
#include "../../../../lib/apply_cuts.cpp"
#include "../../../../lib/admin_utils.h"
#include "../../../../lib/plotting_params.h"



// !./bashing/run_file.sh reconstruct_all_data/rho_reconstruction/elastic_background/cut_low_px

/** make lower bound cut on proton transverse momentum to see if has any effect on 
 * di-rho invariant mass distribution 
 */

int main() {

	// relevent constnats
	const Float_t pion_mass = 1e-3 * m_pi_char; // change to GeV)	
	const Float_t rho_mass = m_rho; // change to GeV
	const Float_t mass_bound = 125; // MeV	
	// Load in aggregated uncut data
	ROOT::EnableImplicitMT();
	TString topology = "40";
	TString file = "data/glueball_mass_reconstruction/data_" + topology + "_SNR.root";
	RDF df_df("tree", file);

	//std::cout << df_df.Describe() << std::endl;
	
	RN df = df_df;

	Float_t n_bins = 100;
	RVecF pt_cuts = {0.0, 0.25, 0.3, 0.4}; // px ~ 0.5GeV py ~ 0.25 GeV 

	TCanvas* c = new TCanvas("c", "c", 600, 800);
	TLegend* leg = new TLegend(0.7,0.45,  0.9,0.9, "Select pr_pt Cuts (MeV/c)");
	auto stack = new THStack();

	//cut on low transverse momenta
	Int_t i = 0;
	for (const Float_t val : pt_cuts) {
    	
		TH1F* h_pt = new TH1F(TString::Format("hist_cut_%d", i), "h_pt", n_bins, 1400, 2500);
		
		auto df_cut = df.Filter([val] (const Double_t px_a, const Double_t py_a, const Double_t px_b, const Double_t py_b) {
			Bool_t a = px_a*px_a > val*val;
			Bool_t b = px_b*px_b > val*val;
			return a && b;
			}, {"pr_ptx_a", "pr_pty_a", "pr_ptx_b", "pr_pty_b"}
		);

    	df_cut.Foreach([&](const RVecLor& p4) {
        	fillHistFromP4(p4, h_pt, rho_mass, mass_bound);
    	}, {"four_momentum"});
	
    	h_pt->SetLineColor(kBlack);
    	h_pt->SetFillColor(color_scheme[i % color_scheme.size()]);
    	leg->AddEntry(h_pt, TString::Format("|pr_pt| > %.3g", val * 1e3), "f");
    	stack->Add(h_pt);
    	++i;
	}
	
	stack->Draw("HIST F nostack");
	stack->SetTitle("Invariant Mass;M (MeV/c^{2});" + TString::Format("Events [%.2g MeV/c^{2}]", stack->GetHistogram()->GetXaxis()->GetBinWidth(1)));
	leg->Draw();
	c->Update();

	SaveOnlyCanvas(c, TString::Format("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/elastic_background/px/px_" + topology + "low_pr_px_cut_mass_bound = %.4g MeV_c^2.root", mass_bound), "RECREATE");



	return 0;
}
