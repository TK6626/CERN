#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"
#include "TLorentzVector.h"

#include "../../../lib/cut_branch.h"
#include "../../../lib/computations.h"
#include "../../../lib/apply_cuts.cpp"
#include "../../../lib/plotting_params.h"
#include "../../../lib/admin_utils.h"

// !./bashing/run_file.sh phi_phi_reconstruction/glueball_reconstruction/glueball_trk_dxy


/** make cuts on data to 
 * try and remove extraneous kaon backroudn
 * due arising from secodary interactions
 */
int main() {

	
	// constants 
	Float_t kaon_mass = m_kaon_char / 1e3;
	Float_t phi_mass = m_phi /1e3;
	RVecF mass= {50, 40, 30, 20, 10, 5};
	mass /= 1e3;
	
	RVecF track_bound = {1, 0.7, 0.5, 0.3, 0.1, 0.05};
	const char* branch = "trk_dxy";




	for (const Float_t bound : track_bound) {
	for (const Float_t mass_bound : mass) {

	Float_t lower_bound = -bound;
	Float_t upper_bound = bound; 
	RVecStr topology = {"20", "40"};
	ROOT::EnableImplicitMT();
	for (TString topo : topology) {
		
		TString file = TString::Format("data/phi_phi_reconstruction/uncut_SNR_"	+ topo + ".root");
		RDF df_df("tree", file);
	
		TCanvas* c = new TCanvas("c");
		TH1F* h = new TH1F("h", "h", 200, 2000, 3000);
	
		RN df = df_df;
		df.Filter( 
			[lower_bound, upper_bound](const RVecF& q) {
        	for (const Float_t val : q) {
            	if (!(lower_bound < val && val < upper_bound)) {
                	return false;
            	}
        	}
        	return true;
		}, {branch}
		)
		.Foreach(
			[&h, mass_bound, phi_mass](const RVecLorCyl p) {
			
			if (
				(TMath::Abs((p[0].M() - phi_mass)) < mass_bound)  
				&& 
				(TMath::Abs((p[1].M() - phi_mass)) < mass_bound)
			) {h->Fill((p[0] + p[1]).M() * 1e3);}
			
			if (
				(TMath::Abs((p[2].M() - phi_mass)) < mass_bound)  
				&& 
				(TMath::Abs((p[3].M() - phi_mass)) < mass_bound)
			) {h->Fill((p[2] + p[3]).M() * 1e3);}
			}, {"phi_four_momentum"}
		);
		h->SetTitle(TString::Format("Invariant X Mass;M (MeV/c^{2});Events [%.2g MeV]", h->GetXaxis()->GetBinWidth(1)));
		TLine* l1 = DrawLine(2220, 0, h, kVertical, kRed, 3, kDashed);
		RVecDraw dat = {{h, "HIST"}, {l1, ""}};
		SaveCanvas(c, dat, TString::Format("media/root_files/phi_phi_reconstruction/glueball_reconstruction/%s/"+ topo +"/<%.3gcm/mass_bound=%.3g_MeV.root",branch, bound,  mass_bound*1e3), "RECREATE");
		
		delete h;
		delete c;
	}
	}
	}
	return 0;
}



