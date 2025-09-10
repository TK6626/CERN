#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"


#include "../../../lib/computations.h"
#include "../../../lib/admin_utils.h"
#include "../../../lib/plotting_params.h"


// !./bashing/run_file.sh phi_phi_reconstruction/protons/delta_phi

int main() {

	SetPlotStyle();
	// set up transparant colours for use later
	Float_t kaon_mass = m_kaon_char / 1e3;
	Float_t phi_mass = m_phi / 1e3;
	Float_t Pi = pi;
	RVecF mass = {10, 20, 30, 50, 50000};
	mass /= 1e3; 

	RVecStr topology = {"20", "40"};
	for (Float_t mass_bound : mass) {
	ROOT::EnableImplicitMT();
	
		
	for (TString topo : topology) {

		TString file = TString::Format("data/phi_phi_reconstruction/uncut_SNR_"	+ topo + ".root");
		RDF df_df("tree", file);
		TCanvas* c1 = new TCanvas("c1", "c1", 600, 800);
		TH1F* hist_both = new TH1F("hist_both", "", 200, -Pi, Pi);
		gStyle->SetOptStat(111111); // include overflow bin counter to verify info not lost

		RN df = df_df.Filter([phi_mass, mass_bound] (const RVecLorCyl p) {
			return 	( 
					((TMath::Abs(p[0].M() - phi_mass) < mass_bound) 
					&& (TMath::Abs(p[1].M() - phi_mass) < mass_bound))

					||
					
					((TMath::Abs(p[2].M() - phi_mass) < mass_bound) 
					&& (TMath::Abs(p[3].M() - phi_mass) < mass_bound))
				);
			},{"phi_four_momentum"});
		
		df.Foreach(
			[&hist_both, Pi] (const Double_t px_a, const Double_t py_a, const Double_t px_b, const Double_t py_b) {	
				Double_t phi_a = TMath::ATan2(py_a, px_a);
				Double_t phi_b = TMath::ATan2(py_b, px_b);
				Double_t d_phi = std::fmod(phi_b - phi_a + 3*Pi, 2*Pi) - Pi;
				
				hist_both->Fill(d_phi);


			}, {"pr_px_a", "pr_py_a", "pr_px_b", "pr_py_b"}
		);
	
		hist_both->SetTitle(TString::Format("Proton #Delta #phi;#Delta #phi (rad); Events [%.3g]", hist_both->GetXaxis()->GetBinWidth(1)));
		RVecDraw dat = {{hist_both, "HIST"}};
		SaveCanvas(c1, dat, TString::Format("media/root_files/phi_phi_reconstruction/protons/delta_phi/"+topo+"/mass_cut=%.4gMeV.root", mass_bound*1e3), "UPDATE");
		
		delete c1;
		delete hist_both;
	}
	}
	return 0;
}
