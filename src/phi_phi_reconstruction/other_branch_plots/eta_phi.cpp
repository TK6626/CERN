#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"

#include "../../../lib/computations.h"
#include "../../../lib/admin_utils.h"
#include "../../../lib/plotting_params.h"


// !./bashing/run_file.sh phi_phi_reconstruction/other_branch_plots/eta_phi	

int main() {

	SetPlotStyle();
	// set up transparant colours for use later
	Float_t kaon_mass = m_kaon_char / 1e3;
	Float_t phi_mass = m_phi / 1e3;
	RVecF mass = {20, 50, 50000}; 
	mass /=1e3;

	gStyle->SetOptStat(0);

	RVecStr topology = {"20", "40"};
	ROOT::EnableImplicitMT();

	for (Float_t mass_bound : mass) {
		
	for (TString topo : topology) {

		TString file = TString::Format("data/phi_phi_reconstruction/uncut_SNR_"	+ topo + ".root");
		RDF df_df("tree", file);
		TCanvas* c1 = new TCanvas("c1", "c1", 600, 800);
		TCanvas* c2 = new TCanvas("c2", "c2", 600, 800);
		TCanvas* c3 = new TCanvas("c3", "c3", 600, 800);
		TH2F* hist_pair1 = new TH2F("hist_pair1", "", 200, -7, 7, 200, -3.2, 3.2);
		TH2F* hist_pair2 = new TH2F("hist_pair2", "", 200, -7, 7, 200, -3.2, 3.2);
		TH2F* hist_both = new TH2F("hist_both", "", 200, -7, 7, 200, -3.2, 3.2);


		RN df1 = df_df.Filter([phi_mass, mass_bound] (const RVecLorCyl p) {
			return 	( 
						(TMath::Abs(p[0].M() - phi_mass) < mass_bound) 
						|| (TMath::Abs(p[1].M() - phi_mass) < mass_bound)
					);
			},{"phi_four_momentum"});
		
		df1.Foreach(
			[&hist_pair1, &hist_both] (const RVecLorCyl p) {	
				
				
				Float_t eta0 = p[0].Eta();
				Float_t eta1 = p[1].Eta();
				Float_t phi0 = p[0].Phi(); 
				Float_t phi1 = p[1].Phi(); 
				
				hist_pair1->Fill(eta0, phi0);
				hist_pair1->Fill(eta1, phi1);
				hist_both->Fill(eta0, phi0);
				hist_both->Fill(eta1, phi1);
			}, {"phi_four_momentum"}
		);
	
		RN df2 = df_df.Filter([phi_mass, mass_bound] (const RVecLorCyl p) {
			return 	( 
						(TMath::Abs(p[2].M() - phi_mass) < mass_bound) 
						|| (TMath::Abs(p[3].M() - phi_mass) < mass_bound)
					);
			},{"phi_four_momentum"});
		
		df2.Foreach(
			[&hist_pair2, &hist_both] (const RVecLorCyl p) {	
				
				Float_t eta0 = p[2].Eta();
				Float_t eta1 = p[3].Eta();
				Float_t phi0 = p[2].Phi(); 
				Float_t phi1 = p[3].Phi(); 
				
				hist_pair2->Fill(eta0, phi0);
				hist_pair2->Fill(eta1, phi1);
				hist_both->Fill(eta0, phi0);
				hist_both->Fill(eta1, phi1);
			}, {"phi_four_momentum"}
		);
	
		hist_pair1->SetTitle(TString::Format("Total #phi-#eta Comparison;Total #eta [%.3g]; Total #phi (rad) [%.3g]", hist_pair1->GetXaxis()->GetBinWidth(1),  hist_pair1->GetYaxis()->GetBinWidth(1)));
		hist_pair2->SetTitle(TString::Format("Total #phi-#eta Comparison;Total #eta [%.3g]; Total #phi (rad) [%.3g]", hist_pair2->GetXaxis()->GetBinWidth(1),  hist_pair2->GetYaxis()->GetBinWidth(1)));
		hist_both->SetTitle(TString::Format("Total #phi-#eta Comparison;Total #eta [%.3g]; Total #phi (rad) [%.3g]", hist_both->GetXaxis()->GetBinWidth(1),  hist_both->GetYaxis()->GetBinWidth(1)));
	
		RVecDraw dat1 = {{hist_pair1, "COLZ"}};
		RVecDraw dat2 = {{hist_pair2, "COLZ"}};
		RVecDraw dat3 = {{hist_both, "COLZ"}};
		
		SaveCanvas(c1, dat1, TString::Format("media/root_files/phi_phi_reconstruction/other_branch_plots/eta_phi/"+topo+"/mass_cut=%.4gMeV.root", mass_bound*1e3), "RECREATE");
		SaveCanvas(c2, dat2, TString::Format("media/root_files/phi_phi_reconstruction/other_branch_plots/eta_phi/"+topo+"/mass_cut=%.4gMeV.root", mass_bound*1e3), "UPDATE");
		SaveCanvas(c3, dat3, TString::Format("media/root_files/phi_phi_reconstruction/other_branch_plots/eta_phi/"+topo+"/mass_cut=%.4gMeV.root", mass_bound*1e3), "UPDATE");
		
		delete c1;
		delete c2;
		delete c3;
		delete hist_pair1;
		delete hist_pair2;
		delete hist_both;
	}
	}
	return 0;
}



