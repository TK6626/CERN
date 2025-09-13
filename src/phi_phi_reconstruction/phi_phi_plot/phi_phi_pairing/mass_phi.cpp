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


// !./bashing/run_file.sh phi_phi_reconstruction/phi_phi_plot/phi_phi_pairing/mass_phi

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
		TH2F* hist1 = new TH2F("hist2", "", 150, 0, 3.2, 150, 2, 5);
		TH2F* hist2 = new TH2F("hist1", "", 150, 0, 3.2, 150, 2, 5);
		TH2F* hist_both = new TH2F("hist_both", "", 150, 0, 3.2, 150, 2, 5);
		
		RN df1 = df_df;
		RN df2 = df_df;

		df1.Filter([phi_mass, mass_bound] (const RVecLorCyl p) {
			return 	( 
						(TMath::Abs(p[0].M() - phi_mass) < mass_bound) 
						|| (TMath::Abs(p[1].M() - phi_mass) < mass_bound)
					);
			},{"phi_four_momentum"})
			.Foreach(
			[&hist_both, &hist1, Pi] (const RVecLorCyl p) {	
				Float_t phi1 = std::fmod(p[0].Phi() - p[1].Phi() + 3*Pi, Pi);
				hist1->Fill(phi1,(p[0] + p[1]).M());
				hist_both->Fill(phi1,(p[0] + p[1]).M());
			}, {"phi_four_momentum"}
		);

		df2.Filter([phi_mass, mass_bound] (const RVecLorCyl p) {
			return 	( 
						(TMath::Abs(p[2].M() - phi_mass) < mass_bound) 
						|| (TMath::Abs(p[3].M() - phi_mass) < mass_bound)
					);
			},{"phi_four_momentum"})
			.Foreach(
			[&hist_both, &hist2, Pi] (const RVecLorCyl p) {	
				Float_t phi2 = std::fmod(p[2].Phi() - p[3].Phi() + 3*Pi, Pi);	
				hist2->Fill(phi2,(p[2] + p[3]).M());
				hist_both->Fill(phi2, (p[2] + p[3]).M());
			}, {"phi_four_momentum"}
		);		


		hist_both->SetTitle(TString::Format("#Delta#phi Between Possible #phi#phi Pairings;#Delta#phi (rad);Events  [%.3g]", hist_both->GetXaxis()->GetBinWidth(1)));
		hist1->SetTitle(TString::Format("#Delta#phi Between Possible #phi#phi Pairings;#Delta#phi (rad);Events  [%.3g]", hist1->GetXaxis()->GetBinWidth(1)));
		hist2->SetTitle(TString::Format("#Delta#phi Between Possible #phi#phi Pairings;#Delta#phi (rad);Events  [%.3g]", hist2->GetXaxis()->GetBinWidth(1)));
	
		RVecDraw dat1 = {{hist1, "HIST"}};
		RVecDraw dat2 = {{hist2, "HIST"}};
		RVecDraw dat3 = {{hist_both, "HIST"}};
		
		SaveCanvas(c1, dat1, TString::Format("media/root_files/phi_phi_reconstruction/phi_phi_plot/mass_delta_phi/"+topo+"/mass_cut=%.4gMeV.root", mass_bound*1e3), "RECREATE");
		SaveCanvas(c2, dat2, TString::Format("media/root_files/phi_phi_reconstruction/phi_phi_plot/mass_delta_phi/"+topo+"/mass_cut=%.4gMeV.root", mass_bound*1e3), "UPDATE");
		SaveCanvas(c3, dat3, TString::Format("media/root_files/phi_phi_reconstruction/phi_phi_plot/mass_delta_phi/"+topo+"/mass_cut=%.4gMeV.root", mass_bound*1e3), "UPDATE");
		
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




