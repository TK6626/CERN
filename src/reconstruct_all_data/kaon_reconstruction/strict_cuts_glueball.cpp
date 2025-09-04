
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TLine.h"

#include "../../../lib/cut_branch.h"
#include "../../../lib/computations.h"
#include "../../../lib/apply_cuts.cpp"
#include "../../../lib/plotting_params.h"
#include "../../../lib/admin_utils.h"

// !./bashing/run_file.sh reconstruct_all_data/kaon_reconstruction/strict_cuts_glueball


/** make cuts on data to 
 * try and remove extraneous kaon backroudn
 * due arising from secodary interactions
 */
int main() {

	// relevent constnats
	RVecF mass_bound_list = {125, 62, 31, 15, 7};
	mass_bound_list /=1e3;
	RVecI index = {0, 1, 2, 3, 4};
	RVecStr topologies = {"20", "40"};
	RVecStr kaon_types = {"*", "s"};
	
	for (const auto& K_type : kaon_types) {
    	Float_t lower_bound = 0.0;
		Float_t kaon_mass = 0.0;
    	if (K_type == "*") {kaon_mass = m_kaon_star * 1e-3; lower_bound = 700;} 
		else if (K_type == "s") {kaon_mass = m_kaon_0 * 1e-3; lower_bound = 1500;}
	for (const TString topology : topologies) {
	for (const Int_t j : index) {
	Float_t mass_bound = mass_bound_list[j];
	//Float_t kaon_mass = m_kaon_star * 1e-3;
	//TString K_type = "*";
	//TString topology = "40";

	// Load in aggregated uncut data
	
	ROOT::EnableImplicitMT();
	
	TString file = "data/glueball_mass_reconstruction/data_" + topology + "_SNR.root";
	RDF df_df("tree", file);
	RN df = CutApplier(df_df)
		.apply_ntrk_cut(4)
		.apply_eta_cut(1.5)
    	.apply_pt_cut(0.7)
    	.apply_charge_cut()
		.apply_dxy_cut(-0.2, 0.2)
		.apply_dz_cut(-2.38, 2.38)
		.apply_zPV_cut(-13.1, 13.1)
		.result();

	TCanvas* c = new TCanvas("c", "c", 600, 800);
	TH1F* hist = new TH1F("hist", "Invariant Mass;M (MeV/c^{2});Counts", 200, lower_bound, 2600);
	
	//plotting invariant mass of the KK
    df.Foreach([&hist, kaon_mass, mass_bound] (const RVecLor& p4) {
        	if (TMath::Abs(p4[0].M() - kaon_mass) < mass_bound 
			&& TMath::Abs(p4[1].M() - kaon_mass) < mass_bound) {
            	hist->Fill((p4[0] + p4[1]).M() * 1e3); // plot in MeV
       	 	}
			if (TMath::Abs(p4[2].M() - kaon_mass) < mass_bound 
			&& TMath::Abs(p4[3].M() - kaon_mass) < mass_bound) {
            	hist->Fill((p4[2] + p4[3]).M() * 1e3); // plot in MeV
        	}
		}, {"four_momentum"});

	TLine* line = DrawVLine(2220, hist);
	line->Draw();
	hist->Draw("HIST F");
	c->Update();
	SaveOnlyCanvas(c, TString::Format("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/kaon_reconstruction/strict_cuts/K_" + K_type + "_glueball_" + topology + "/mass bound = %.4g MeV_c^2.root", mass_bound * 1e3), "RECREATE"); 
	delete c;
	delete hist;
	}
	}
	}
	return 0;
}
