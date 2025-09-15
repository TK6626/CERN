#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"

#include "../../../lib/plotting_params.h"
#include "../../../lib/admin_utils.h"



// !./bashing/run_file.sh phi_phi_reconstruction/other_branch_plots/eta_phi_single

int main() {

	SetPlotStyle();
	// set up transparant colours for use later
	Float_t kaon_mass = m_kaon_char / 1e3;
	Float_t phi_mass = m_phi / 1e3;
	RVecF mass = {20, 50, 5000}; 
	mass /=1e3;





	RVecStr topology = {"4"};
	ROOT::EnableImplicitMT();

	for (Float_t mass_bound : mass) {
	for (TString topo : topology) {
		


		TString file = TString::Format("data/phi_phi_reconstruction/uncut_SNR_"	+ topo + ".root");
		RDF df_df("tree", file);

		TCanvas* c1 = new TCanvas("c1", "c1", 600, 800);
		TH2F* hist_pair1 = new TH2F("hist_pair1", "", 200, -3.2, 3.2, 200, -3.2, 3.2);

		c1->SetTopMargin(0.02);
		c1->SetRightMargin(0.1);
		hist_pair1->GetXaxis()->SetTitleOffset(1.2);
		gStyle->SetOptStat(1110); 

		RN df1 = df_df.Filter([phi_mass, mass_bound] (const RVecLorCyl p) {
			return 	( 
						(TMath::Abs(p[0].M() - phi_mass) < mass_bound) 
						|| (TMath::Abs(p[1].M() - phi_mass) < mass_bound)
					);
			},{"phi_four_momentum"});
		
		df1.Foreach(
			[&hist_pair1] (const RVecLorCyl p) {	
				for (Int_t i; i < 4; ++i) {
				hist_pair1->Fill(p[i].Eta(), p[i].Phi());
				}
			}, {"kaon_four_momentum"}
		);

		hist_pair1->SetTitle(TString::Format(";#eta [%.3g]; #phi (rad) [%.3g]", hist_pair1->GetXaxis()->GetBinWidth(1),  hist_pair1->GetYaxis()->GetBinWidth(1)));
		RVecDraw dat1 = {{hist_pair1, "COLZ"}};
		
		SaveCanvas(c1, dat1, TString::Format("media/root_files/phi_phi_reconstruction/kaons/eta_phi/"+topo+"/mass_cut=%.4gMeV.root", mass_bound*1e3), "RECREATE");
		delete c1;
		delete hist_pair1;
	}
	}
	return 0;
}



