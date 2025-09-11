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


// !./bashing/run_file.sh phi_phi_reconstruction/phi_phi_plot/phi_phi_pairing/delta_phi

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
	gStyle->SetOptStat(0);	
	
	for (Float_t mass_bound : mass) {
	for (TString topo : topology) {

		TString file = TString::Format("data/phi_phi_reconstruction/uncut_SNR_"	+ topo + ".root");
		RDF df_df("tree", file);
		TCanvas* c1 = new TCanvas("c1", "c1", 600, 800);
		TCanvas* c2 = new TCanvas("c2", "c3", 600, 800);
		TCanvas* c3 = new TCanvas("c3", "c3", 600, 800);
		TH1F* hist_both = new TH1F("hist_both", "", 200, 0, Pi);
		TH1F* hist1 = new TH1F("hist2", "", 800, -3.2*1e-3, 3.2*1e-3);
		TH1F* hist2 = new TH1F("hist1", "", 800, -3.2*1e-3, 3.2*1e-3);


		RN df1 = df_df.Filter([phi_mass, mass_bound] (const RVecLorCyl p) {
			return 	( 
						(TMath::Abs(p[0].M() - phi_mass) < mass_bound) 
						|| (TMath::Abs(p[1].M() - phi_mass) < mass_bound)
					);
			},{"phi_four_momentum"});
		
		df1.Foreach(
			[&hist_both, &hist1, Pi] (const RVecLorCyl p) {	
				
				
				Float_t phi1 = std::fmod(p[0].Phi() - p[1].Phi() + 3*Pi, Pi);
				Float_t phi2 = std::fmod(p[2].Phi() - p[3].Phi() + 3*Pi, Pi);	
				hist_both->Fill(phi1);
				hist_both->Fill(phi2);
			}, {"phi_four_momentum"}
		);
				
		hist_both->SetTitle(TString::Format("#Delta#phi Between Possible #phi#phi Pairings;#Delta#phi (rad);Events  [%.3g]", hist_both->GetXaxis()->GetBinWidth(1)));
	
		RVecDraw dat = {{hist_both, "HIST"}};
		
		SaveCanvas(c1, dat, TString::Format("media/root_files/phi_phi_reconstruction/phi_phi_plot/delta_phi/"+topo+"/mass_cut=%.4gMeV.root", mass_bound*1e3), "RECREATE");
		
		delete c1;
		delete c2;
		delete c3;
		delete hist_both;
		delete hist1;
		delete hist2;
	}
	}
	return 0;
}




