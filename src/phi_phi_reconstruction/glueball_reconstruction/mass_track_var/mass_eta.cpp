#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"


#include "../../../../lib/computations.h"
#include "../../../../lib/admin_utils.h"
#include "../../../../lib/plotting_params.h"


// !./bashing/run_file.sh phi_phi_reconstruction/glueball_reconstruction/mass_track_var/mass_eta


int main() {

	SetPlotStyle();
	// set up transparant colours for use later
	Float_t kaon_mass = m_kaon_char / 1e3;
	Float_t phi_mass = m_phi / 1e3;
	Float_t mass_bound = 35 / 1e3; 

	RVecStr topology = {"20", "40"};
	ROOT::EnableImplicitMT();
	
		
	for (TString topo : topology) {

		TString file = TString::Format("data/phi_phi_reconstruction/uncut_SNR_"	+ topo + ".root");
		RDF df_df("tree", file);
		TCanvas* c1 = new TCanvas("c1", "c1", 600, 800);
		TCanvas* c2 = new TCanvas("c2", "c2", 600, 800);
		TCanvas* c3 = new TCanvas("c3", "c3", 600, 800);
		TH2F* hist_pair1 = new TH2F("hist_pair1", "", 200, -7, 7, 200, 2000, 3500);
		TH2F* hist_pair2 = new TH2F("hist_pair2", "", 200, -7, 7, 200, 2000, 3500);
		TH2F* hist_both = new TH2F("hist_both", "", 200, -7, 7, 200, 2000, 3500);


		RN df1 = df_df.Filter([phi_mass, mass_bound] (const RVecLorCyl p) {
			return 	( 
						(TMath::Abs(p[0].M() - phi_mass) < mass_bound) 
						|| (TMath::Abs(p[1].M() - phi_mass) < mass_bound)
					);
			},{"phi_four_momentum"});
		
		df1.Foreach(
			[&hist_pair1, &hist_both] (const RVecLorCyl p4, const RVecLorCyl p) {	
				
				
				LorCyl P = p4[0] + p4[1] + p4[2] + p4[3];
				Float_t eta = P.Eta();
				Float_t mass = (p[0] + p[1]).M()*1e3;
				hist_pair1->Fill(eta, mass);
				hist_both->Fill(eta, mass);
			}, {"kaon_four_momentum", "phi_four_momentum"}
		);
	
		RN df2 = df_df.Filter([phi_mass, mass_bound] (const RVecLorCyl p) {
			return 	( 
						(TMath::Abs(p[2].M() - phi_mass) < mass_bound) 
						|| (TMath::Abs(p[3].M() - phi_mass) < mass_bound)
					);
			},{"phi_four_momentum"});
		
		df2.Foreach(
			[&hist_pair2, &hist_both] (const RVecLorCyl p4, const RVecLorCyl p) {	
				
				LorCyl P = p4[0] + p4[1] + p4[2] + p4[3];
				Float_t eta = P.Eta();
				Float_t mass = (p[2] + p[3]).M() * 1e3;
				hist_pair2->Fill(eta, mass);
				hist_both->Fill(eta, mass);
			}, {"kaon_four_momentum", "phi_four_momentum"}
		);

		hist_pair1->SetTitle(TString::Format("Mass #eta Comparison;Total #eta [%.3g]; X Mass (MeV/c^{2}) [%.3g]", hist_pair1->GetXaxis()->GetBinWidth(1),  hist_pair1->GetYaxis()->GetBinWidth(1)));
		hist_pair2->SetTitle(TString::Format("Mass #eta Comparison;Total #eta [%.3g]; X Mass (MeV/c^{2}) [%.3g]", hist_pair2->GetXaxis()->GetBinWidth(1),  hist_pair2->GetYaxis()->GetBinWidth(1)));
		hist_both->SetTitle(TString::Format("Mass #eta Comparison;Total #eta [%.3g]; X Mass (MeV/c^{2}) [%.3g]", hist_both->GetXaxis()->GetBinWidth(1),  hist_both->GetYaxis()->GetBinWidth(1)));
	
		RVecDraw dat1 = {{hist_pair1, "COLZ"}};
		RVecDraw dat2 = {{hist_pair2, "COLZ"}};
		RVecDraw dat3 = {{hist_both, "COLZ"}};
		
		SaveCanvas(c1, dat1, TString::Format("media/root_files/phi_phi_reconstruction/glueball_reconstruction/mass_track_var/eta/"+topo+"/mass_cut=%.4g.root", mass_bound*1e3), "RECREATE");
		SaveCanvas(c2, dat2, TString::Format("media/root_files/phi_phi_reconstruction/glueball_reconstruction/mass_track_var/eta/"+topo+"/mass_cut=%.4g.root", mass_bound*1e3), "UPDATE");
		SaveCanvas(c3, dat3, TString::Format("media/root_files/phi_phi_reconstruction/glueball_reconstruction/mass_track_var/eta/"+topo+"/mass_cut=%.4g.root", mass_bound*1e3), "UPDATE");
		
		delete c1;
		delete c2;
		delete hist_pair1;
		delete hist_pair2;
	}
	return 0;
}


