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
	RVecF mass = {10, 5e5}; 
	mass/=1e3;



	RVecStr topology = {"40"};
	ROOT::EnableImplicitMT();
	
	for(const Float_t mass_bound : mass) {	
	for (TString topo : topology) {

		TString file = TString::Format("data/phi_phi_reconstruction/uncut_SNR_"	+ topo + ".root");
		RDF df_df("tree", file);
		TCanvas* c1 = new TCanvas("c1");
		TCanvas* c2 = new TCanvas("c2");
		TCanvas* c3 = new TCanvas("c3");
		TH2F* h_2d = new TH2F("h_2d", "", 200, -5, 5, 200, 0, 2);
		TH1F* h_pt = new TH1F("pt", "", 200, 0, 2); 
		TH1F* h_eta = new TH1F("eta", "", 200, -5, 5); 

		c1->SetTopMargin(0.02);
		c2->SetRightMargin(0.12);
		h_2d->GetXaxis()->SetTitleOffset(1.2);
		gStyle->SetOptStat(1000);


		RN df1 = df_df.Filter([phi_mass, mass_bound] (const RVecLorCyl p) {
			return 	( 
						(TMath::Abs(p[0].M() - phi_mass) < mass_bound) 
						|| (TMath::Abs(p[1].M() - phi_mass) < mass_bound)
				);
			},{"phi_four_momentum"});

		df1.Foreach(
			[&h_2d, &h_pt, h_eta] (const RVecLorCyl p4) {	
				
				LorCyl P = p4[0] + p4[1];
				Float_t eta = P.Eta();
				Float_t pt = P.Pt();

				h_2d->Fill(eta, pt);
				h_pt->Fill(pt);
				h_eta->Fill(eta);
			}, {"phi_four_momentum"}
		);
		
		RN df2 = df_df.Filter([phi_mass, mass_bound] (const RVecLorCyl p) {
			return 	( 
						(TMath::Abs(p[2].M() - phi_mass) < mass_bound) 
						|| (TMath::Abs(p[3].M() - phi_mass) < mass_bound)
				);
			},{"phi_four_momentum"});

		df2.Foreach(
			[&h_2d, &h_pt, h_eta] (const RVecLorCyl p4) {	
				
				LorCyl P = p4[2] + p4[3];
				Float_t eta = P.Eta();
				Float_t pt = P.Pt();

				h_2d->Fill(eta, pt);
				h_pt->Fill(pt);
				h_eta->Fill(eta);
			}, {"phi_four_momentum"}
		);


		h_2d->SetTitle(TString::Format(";#phi p_t (GeV/c) [%.3g]; #phi #eta [%.3g]", h_2d->GetXaxis()->GetBinWidth(1),  h_2d->GetYaxis()->GetBinWidth(1)));
		h_pt->SetTitle(TString::Format("Total p_t Distribution ;Total p_t (GeV/c);Events [%.3g]", h_pt->GetBinWidth(1)));
		h_eta->SetTitle(TString::Format("Total #eta Distribution ;Total #eta (GeV/c);Events [%.3g]", h_eta->GetBinWidth(1)));
		
		SaveCanvas(c1, h_2d, "COLZ", TString::Format("media/root_files/phi_phi_reconstruction/other_branch_plots/pt_eta/"+topo+"/mass_cut=%.4g_one/compare.root", mass_bound*1e3), "RECREATE");
		SaveCanvas(c2, h_pt, "HIST", TString::Format("media/root_files/phi_phi_reconstruction/other_branch_plots/pt_eta/"+topo+"/mass_cut=%.4g_one/pt_distribution.root", mass_bound*1e3), "RECREATE");
		SaveCanvas(c3, h_eta, "HIST", TString::Format("media/root_files/phi_phi_reconstruction/other_branch_plots/pt_eta/"+topo+"/mass_cut=%.4g_one/eta_distribution.root", mass_bound*1e3), "RECREATE");
		
		delete c1;
		delete c2;
		delete c3;
		delete h_2d;
		delete h_pt;
		delete h_eta;
	}
	}
	return 0;
}

