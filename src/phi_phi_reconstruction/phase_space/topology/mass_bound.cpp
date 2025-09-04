#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "THStack.h"

#include "../../../../lib/cut_branch.h"
#include "../../../../lib/computations.h"
#include "../../../../lib/apply_cuts.cpp"
#include "../../../../lib/plotting_params.h"
#include "../../../../lib/admin_utils.h"

// !./bashing/run_file.sh phi_phi_reconstruction/phase_space/topology/mass_bound


/** make cuts on data to 
 */
int main() {

	
	// constants 
	Float_t kaon_mass = m_kaon_char / 1e3; //GeV
	Float_t phi_mass = m_phi /1e3;
	
	const char* branch = "trk_pt"; // what we are looking at how it varies due to mass cuts 

	// Load in aggregated uncut data
	ROOT::EnableImplicitMT();
	
	RVecStr topology = {"20", "40"};
	RVecI mass_cuts = {30, 25, 20, 15, 10};	

	TH1* h = new TH1F("h", "p_{t} Distribution of 4-track;p_{t} MeV/c;Events", 100, 0, 2);
	TCanvas* c = new TCanvas("c");
	TLegend* leg = new TLegend(0.7,0.6, 0.9,0.8, "Select mass cuts");
	THStack* stack = new THStack("stack", "");

	for (TString topo : topology) {
		TString file = "data/phi_phi_reconstruction/uncut_SNR_" + topo+ ".root";
		RDF df_df("tree", file);
		RN df = df_df;

		Int_t i = 0;
		for (Int_t val : mass_cuts) {
			df.Filter( [&val, &phi_mass] (const RVecLorCyl p) { // filter for mass cuts around the resonance
				return (TMath::Abs(p[0].M() - phi_mass) < val && TMath::Abs(p[1].M() < val));
			}, {"phi_four_momentum"}
			)
			.Foreach( [&h] (const RVecF p4) {
				for (Float_t p : p4){
					h->Fill(p);
				}
			}, {branch}
			);
		
		h->SetLineColor(kBlack);
    	h->SetFillColor(color_scheme[i % color_scheme.size()]);
    	leg->AddEntry(h, TString::Format("%s  < %.3g", branch, val), "f");
    	stack->Add(h);
		++i;	
		
		h->Reset();
		}

		
	stack->Draw("HIST F nostack");
	stack->SetTitle("Invariant Mass;M (MeV/c^{2});" + TString::Format("Events [%.2g MeV/c^{2}]", stack->GetHistogram()->GetXaxis()->GetBinWidth(1)));
	RVecDraw dat = {{stack, "HIST F nostack"}, {leg, ""}};
	SaveCanvas(c, dat, TString::Format("media/root_files/phi_phi_reconstruction/phase_space/mass_cuts/%s/topology=" + topo + ".root", branch), "RECREATE");



	}
	return 0;
}

