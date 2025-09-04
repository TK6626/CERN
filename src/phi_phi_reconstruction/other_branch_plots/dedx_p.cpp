#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"


#include "../../../../lib/computations.h"
#include "../../../../lib/admin_utils.h"
#include "../../../../lib/plotting_params.h"


// !./bashing/run_file.sh phi_phi_reconstruction/simple_analysis/other_branch_plots/dedx_p


int main() {

	SetPlotStyle();
	
	// set up transparant colours for use later
	
	Float_t kaon_mass = m_kaon_char;
	Float_t phi_mass = m_phi;
	Float_t mass_bound = 100; 


	RVecStr topology = {"20", "40"};
	ROOT::EnableImplicitMT();
	
		
	for (TString topo : topology) {

		TString file = TString::Format("data/phi_phi_reconstruction/uncut_SNR_"	+ topo + ".root");
		RDF df_df("tree", file);
		RN df = df_df;
		TCanvas* c = new TCanvas("c");
		TH2F* h_2d = new TH2F("h_2d", "", 100, 0, 11, 100, -46, 40);
		TH1F* h_dedx = new TH1F("dedx", "", 100, -46, 40); 

		df.Foreach(
			[&h_2d, &h_dedx] (const RVecF p, const RVecF dedx) {
				for (Int_t i=0; i < 4; ++i) {
					Float_t ln_dedx = TMath::Log10(dedx[i]);
					
					if (TMath::Finite(ln_dedx)) {
						h_2d->Fill(p[i], ln_dedx);
						h_dedx->Fill(ln_dedx);	
					}
				}
			}, {"trk_p", "trk_dedx"}
		);

		h_dedx->SetTitle(TString::Format("Energy loss; log_{10}(dEdx) (GeV/c dm) ;Events [%.3g GeV/c cm]", h_dedx->GetXaxis()->GetBinWidth(1)));
		h_2d->SetTitle(TString::Format("Energy loss against momentum;p (GeV/c)[%.3g]; log_{10}(dE/dx) (GeV/c dm) [%.3g GeV/c]", h_2d->GetYaxis()->GetBinWidth(1),  h_2d->GetXaxis()->GetBinWidth(1)));
		
		SaveCanvas(c, h_2d, "COLZ", TString::Format("media/root_files/phi_phi_reconstruction/simple_analysis/other_branch_plots/"+topo+"/ln_dedx_p.root"), "RECREATE");
		SaveCanvas(c, h_dedx, "E1P", TString::Format("media/root_files/phi_phi_reconstruction/simple_analysis/other_branch_plots/"+topo+"/ln_dedx.root"), "RECREATE");
		
		delete c;
		delete h_2d;
		delete h_dedx;
	}
	return 0;
}
