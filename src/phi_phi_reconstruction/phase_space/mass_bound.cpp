#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TLegend.h"
#include "THStack.h"

#include "../../../lib/cut_branch.h"
#include "../../../lib/computations.h"
#include "../../../lib/apply_cuts.cpp"
#include "../../../lib/plotting_params.h"
#include "../../../lib/admin_utils.h"

// !./bashing/run_file.sh phi_phi_reconstruction/phase_space/mass_bound


/** make cuts on data to 
 */
int main() {

	
	// constants 
	Float_t kaon_mass = m_kaon_char / 1e3; //GeV
	Float_t phi_mass = m_phi /1e3;
	RVecF mass_cuts = {50, 30, 25, 20, 15, 10};
	mass_cuts *=1e-3;
	Int_t n = mass_cuts.size() + 1;	
	

	const char* branch = "trk_pt"; // what we are looking at how it varies due to mass cuts 

	ROOT::EnableImplicitMT();


	TString topo = "20";
		TString file = "data/phi_phi_reconstruction/uncut_SNR_" + topo + ".root";                         
		RDF df_df("tree", file);                                                                         
		RN df = df_df;                                                                                   

		TCanvas* c = new TCanvas(TString::Format("c%s", topo.Data()), "3D Histograms", 1200, 600);
		TLegend* leg = new TLegend(0.7,0.45,0.9,0.7,"Select mass cuts");
		RVecDraw dat(n);
		
		Int_t i = 0;
		for (Float_t val : mass_cuts) {                                                                    
    		TString hname = TString::Format("h_%d_" + topo, i);
    		TH1F* h = new TH1F(hname, "p_{t} Distribution of 4-track;p_{t} MeV/c;Events", 100, 0, 1000);

    		df.Filter([val, phi_mass](const RVecLorCyl p) { // filter for mass cuts around the resonance
				return ((TMath::Abs(p[0].M() - phi_mass) < val) && (TMath::Abs(p[1].M() - phi_mass) < val));
        	}, {"phi_four_momentum"}
    		)
			.Foreach([&h](const RVecF& p4) {                                                            
				for (Float_t p : p4){                                                                    
                	h->Fill(p *1e3);                                                                          
            	}                                                                                        
        	}, {branch}
    		);                                                                                           

		h->SetLineColor(kBlack);
    	h->SetFillColor(color_scheme[i % color_scheme.size()]);
    	h->SetFillStyle(1001);
    	h->SetBarWidth(0.8);        // bar width in X
    	h->SetBarOffset(i*1.0);    // offset along Y to separate histograms
    	leg->AddEntry(h, TString::Format("#Delta#M < %.3g MeV/c", val * 1e3), "f");
    	Drawable data {h, "hist"};
		dat[i] = data;
		++i;                                                                                             
		}                                                                                                

		Drawable leg_draw = {leg, ""};
		dat[n-1] = leg_draw;
		SaveCanvas(c, dat, TString::Format("media/root_files/phi_phi_reconstruction/phase_space/mass_cuts/%s/topology=" + topo + ".root", branch), "RECREATE");
	
	return 0;
}

