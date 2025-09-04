#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"

#include "../../../../lib/cut_branch.h"
#include "../../../../lib/computations.h"
#include "../../../../lib/apply_cuts.cpp"
#include "../../../../lib/plotting_params.h"
#include "../../../../lib/admin_utils.h"

// !./bashing/run_file.sh phi_phi_reconstruction/phi_phi_plot/rm_bg/p
                

/** make cuts on data to 
 * try and remove extraneous kaon backroudn
 * due arising from secodary interactions
 */
int main() {

	
	// constants 
	Float_t kaon_mass = m_kaon_char / 1e3;
	Float_t phi_mass = m_phi /1e3;
	
	TString branch = "p";
	RVecF cut = {10, 8, 5, 4, 2, 1.5,  1, 0.7, 0.5, 0.3, 0.1};

	RVecStr topology = {"20", "40"};
	ROOT::EnableImplicitMT();
	for (TString topo : topology) {
		for (Float_t val : cut) {
		
		TString file = TString::Format("data/phi_phi_reconstruction/uncut_SNR_"	+ topo + ".root");
		RDF df_df("tree", file);
	
		TCanvas* c = new TCanvas("c");
		TH2F* h = new TH2F("h", "h", 300, 980, 1400, 300, 980, 1400);
	
		RN df = CutApplier(df_df)
			.apply_p_cut(val)
			.result();

		df.Foreach(
			[&h](const RVecLorCyl p) {
			
			h->Fill(p[0].M()*1e3, p[1].M() * 1e3);
			h->Fill(p[2].M()*1e3, p[3].M() * 1e3);
			
			}, {"phi_four_momentum"}
		);
		h->SetTitle(TString::Format("Invariant X Mass;M (MeV/c^{2});Events [%.2g MeV]", h->GetXaxis()->GetBinWidth(1)));
		h->Draw();
		RVecDraw dat = {{h, "COLZ"}};
		SaveCanvas(c, dat, TString::Format("media/root_files/phi_phi_reconstruction/phi_phi_plot/"+branch+"/"+topo+"/%.2f_cut.root", val), "RECREATE");
		
		delete h;
		delete c;
	}
	}
	return 0;
}


