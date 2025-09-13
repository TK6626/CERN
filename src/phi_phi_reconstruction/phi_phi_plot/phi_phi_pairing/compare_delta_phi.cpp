#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include <cmath>

#include "../../../../lib/computations.h"
#include "../../../../lib/admin_utils.h"
#include "../../../../lib/plotting_params.h"


// !./bashing/run_file.sh phi_phi_reconstruction/phi_phi_plot/phi_phi_pairing/compare_delta_phi

int main() {

	SetPlotStyle();
	// set up transparant colours for use later
	Float_t kaon_mass = m_kaon_char / 1e3;
	Float_t phi_mass = m_phi / 1e3;
	RVecF  mass = {25, 50, 50000};
	mass /=1e3; 
	Float_t Pi = pi;

	RVecStr topology = {"20", "40"};
	ROOT::EnableImplicitMT();
	gStyle->SetOptStat(111111);	
	
	for (Float_t mass_bound : mass) {
	for (TString topo : topology) {

		TString file = TString::Format("data/phi_phi_reconstruction/uncut_SNR_"	+ topo + ".root");
		RDF df_df("tree", file);

		TCanvas* c1 = new TCanvas("c1", "c1", 600, 800);
		TH2F* hist1 = new TH2F("hist1", "", 100, 0, 3.2, 100, 0, 3.2);


		RN df1 = df_df.Filter([phi_mass, mass_bound] (const RVecLorCyl p) {
			return 	( 
						(TMath::Abs(p[0].M() - phi_mass) < mass_bound) 
						|| (TMath::Abs(p[1].M() - phi_mass) < mass_bound)
					);
			},{"phi_four_momentum"});
		
		df1.Foreach(
			[&hist1, Pi] (const RVecLorCyl p) {	
				Float_t phi1 = std::fmod(p[0].Phi() - p[1].Phi() + 3*Pi, Pi);
				Float_t phi2 = std::fmod(p[2].Phi() - p[3].Phi() + 3*Pi, Pi);	
				hist1->Fill(phi1, phi2);
			}, {"phi_four_momentum"}
		);
				
		hist1->SetTitle(TString::Format("#Delta#phi Between Possible #phi#phi Pairings;#Delta#phi (rad) Events [%.3g];#Delta#phi (rad) Events [%.3g]", hist1->GetXaxis()->GetBinWidth(1), hist1->GetYaxis()->GetBinWidth(1)));
	
		RVecDraw dat = {{hist1, "HIST"}};
		
		SaveCanvas(c1, dat, TString::Format("media/root_files/phi_phi_reconstruction/phi_phi_plot/compare_delta_phi/"+topo+"/mass_cut=%.4gMeV.root", mass_bound*1e3), "RECREATE");
		
		delete c1;
		delete hist1;
	}
	}
	return 0;
}




