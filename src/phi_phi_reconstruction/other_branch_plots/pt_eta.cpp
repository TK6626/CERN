#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"


#include "../../../lib/computations.h"
#include "../../../lib/admin_utils.h"
#include "../../../lib/plotting_params.h"


// !./bashing/run_file.sh phi_phi_reconstruction/other_branch_plots/pt_eta


int main() {

	SetPlotStyle();
	// set up transparant colours for use later
	Float_t kaon_mass = m_kaon_char / 1e3;
	Float_t phi_mass = m_phi / 1e3;
	Float_t mass_bound = 10 / 1e3; 

	RVecStr topology = {"20", "40"};
	ROOT::EnableImplicitMT();
	
		
	for (TString topo : topology) {

		TString file = TString::Format("data/phi_phi_reconstruction/uncut_SNR_"	+ topo + ".root");
		RDF df_df("tree", file);
		TCanvas* c1 = new TCanvas("c1");
		TCanvas* c2 = new TCanvas("c2");
		TCanvas* c3 = new TCanvas("c3");
		TH2F* h_2d = new TH2F("h_2d", "", 100, -15, 15, 100, 0, 8);
		TH1F* h_pt = new TH1F("pt", "", 100, 0, 8); 
		TH1F* h_eta = new TH1F("eta", "", 100, -15, 15); 


		RN df = df_df.Filter([phi_mass, mass_bound] (const RVecLorCyl p) {
			return 	( 
						(TMath::Abs(p[0].M() - phi_mass) < mass_bound) 
						|| (TMath::Abs(p[1].M() - phi_mass) < mass_bound)
					);
			},{"phi_four_momentum"});

		df.Foreach(
			[&h_2d, &h_pt, h_eta] (const RVecF pt, const RVecF p) {	
						h_2d->Fill(Sum(eta), Sum(pt));
						h_pt->Fill(Sum(pt));
						h_eta->Fill(Sum(eta));
			}, {"trk_pt", "kaon_four_momentum"}
		);

		h_2d->SetTitle(TString::Format("Energy loss against momentum;Total p_t (GeV/c) [%.3g]; Total #eta [%.3g]", h_2d->GetYaxis()->GetBinWidth(1),  h_2d->GetXaxis()->GetBinWidth(1)));
		h_pt->SetTitle(TString::Format("Total p_t Distribution ;Total p_t (GeV/c);Events [%.3g]", h_pt->GetBinWidth(1)));
		h_eta->SetTitle(TString::Format("Total #eta Distribution ;Total #eta (GeV/c);Events [%.3g]", h_eta->GetBinWidth(1)));
		
		SaveCanvas(c1, h_2d, "COL Z", TString::Format("media/root_files/phi_phi_reconstruction/other_branch_plots/pt_eta/"+topo+"/mass_cut=%.4g_one/compare.root", mass_bound*1e3), "RECREATE");
		SaveCanvas(c2, h_pt, "HIST", TString::Format("media/root_files/phi_phi_reconstruction/other_branch_plots/pt_eta/"+topo+"/mass_cut=%.4g_one/pt_distribution.root", mass_bound*1e3), "RECREATE");
		SaveCanvas(c3, h_eta, "HIST", TString::Format("media/root_files/phi_phi_reconstruction/other_branch_plots/pt_eta/"+topo+"/mass_cut=%.4g_one/eta_distribution.root", mass_bound*1e3), "RECREATE");
		
		delete c1;
		delete c2;
		delete c3;
		delete h_2d;
		delete h_pt;
		delete h_eta;
	}
	return 0;
}

