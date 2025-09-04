#include <stdio.h>

#include <iostream>

#include "TH2.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLorentzVector.h"

#include "../../../lib/custom_definitions.h"
#include "../../../lib/cut_branch.h"
#include "../../../lib/computations.h"
#include "../../../lib/apply_cuts.cpp"
#include "../../../lib/plotting_params.h"
#include "../../../lib/admin_utils.h"

// !./bashing/run_file.sh phi_phi_reconstruction/mass_cut_invariant_phi_mass


/** Create simple plot of constructed invariant kaon masses
 * subject to simple mass bound cut
 */

int main(){

	// constants 
	Float_t kaon_mass = m_kaon_char / 1e3;
	Float_t phi_mass = m_phi /1e3;

	RVecStr topology = {"20", "40"};
	ROOT::EnableImplicitMT();
	for (TString topo : topology) {
		
		TString file = TString::Format("data/phi_phi_reconstruction/uncut_SNR_"	+ topo + ".root");
		RDF df_df("tree", file);
		TCanvas* c = new TCanvas("c", "c", 600, 800);
		TH2F* h = new TH2F("h", "Invariant #phi#phi Mass;M_{1} (MeV/c^{2});M_{2} (MeV/c^{2})", 400, 985, 1450, 400, 985, 1450);
	
		RN df = df_df;
		df.Foreach(
			[&h] (const RVecLor p4) {
			h->Fill(p4[0].M()*1e3, p4[1].M()*1e3);
			h->Fill(p4[2].M()*1e3, p4[3].M()*1e3);
		}, {"phi_four_momentum"}
		);
		
		h->Draw("COLZ");
		DrawLine(phi_mass * 1e3, 0, h, kVertical, kRed, 2, kDashed)->Draw("SAME");
		DrawLine(phi_mass * 1e3, 0, h, kHorizontal, kRed, 2, kDashed)->Draw("SAME");
		c->Update();
		SaveOnlyCanvas(c, TString::Format("media/root_files/phi_phi_reconstruction/phi_phi_cuts/general_cuts/" + topo + "/mass_cut_invariant_phi_mass.root"), "RECREATE");
	delete h;
	delete c;
	
	}
	return 0;
}

