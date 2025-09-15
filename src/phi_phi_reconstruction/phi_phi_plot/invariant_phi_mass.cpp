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

// !./bashing/run_file.sh phi_phi_reconstruction/phi_phi_plot/invariant_phi_mass


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
		RN df = df_df;
		TCanvas* c = new TCanvas("c", "c", 600, 800);
		TH2F* h_invar = new TH2F("h_invar", ";M_{1} (MeV/c^{2});M_{2} (MeV/c^{2})", 100, 985, 1300, 100, 985,1300);
		
		c->SetTopMargin(0.02);
		c->SetRightMargin(0.1);
		h_invar->GetXaxis()->SetTitleOffset(1.2);
		gStyle->SetOptStat(1110); 

	
		df.Foreach(
			[&h_invar] (const RVecLorCyl p4) {

			h_invar->Fill(p4[0].M()*1e3, p4[1].M()*1e3);
			h_invar->Fill(p4[2].M()*1e3, p4[3].M()*1e3);
		}, {"phi_four_momentum"}
		);
		
		TLine* l1 = DrawLine(phi_mass * 1e3, 0, h_invar, kVertical, kRed, 2, kDashed);
		TLine* l2 = DrawLine(phi_mass * 1e3, 0, h_invar, kHorizontal, kRed, 2, kDashed);
		c->Update();
		RVecDraw dat = {
			{h_invar, "COLZ"}
	//		, {l1, ""}, {l2, ""}
		};
		SaveCanvas(c, dat, TString::Format("media/root_files/phi_phi_reconstruction/phi_phi_plot/initial_reconstruction/" + topo + "/Invariant_KK_mass.root"), "RECREATE");
	delete h_invar;
	delete c;
	delete l1;
	delete l2;
	
	}
	return 0;
}

