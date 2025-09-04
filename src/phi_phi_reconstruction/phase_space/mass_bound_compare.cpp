#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "TCanvas.h"
#include "TH1.h"
#include "TColor.h"
#include "TLegend.h"
#include "THStack.h"

#include "../../../lib/cut_branch.h"
#include "../../../lib/computations.h"
#include "../../../lib/apply_cuts.cpp"
#include "../../../lib/plotting_params.h"
#include "../../../lib/admin_utils.h"

// !./bashing/run_file.sh phi_phi_reconstruction/phase_space/mass_bound_compare


/** make cuts on data to 
 */
int main() {

	
// constants 
	Float_t kaon_mass = m_kaon_char / 1e3; //GeV
	Float_t phi_mass = m_phi /1e3;
	RVecF mass_cuts = {50, 30, 25, 20, 15, 10, 5};
	mass_cuts *=1e-3;
	Int_t n = mass_cuts.size() + 1;	
	

	const char* branch = "trk_pt"; // what we are looking at how it varies due to mass cuts 

	ROOT::EnableImplicitMT();
	
	TString topo_40 = "40";
	TString topo_20 = "20";
		
	TString file_20 = "data/phi_phi_reconstruction/uncut_SNR_" + topo_20 + ".root";                         
	TString file_40 = "data/phi_phi_reconstruction/uncut_SNR_" + topo_40 + ".root";                         
		
	RDF df_df_40("tree", file_40);                                                                         
		RN df_40 = df_df_40;                                                                                   
	RDF df_df_20("tree", file_20);                                                                         
		RN df_20 = df_df_20;                                                                                   

		
		for (Float_t val : mass_cuts) {                                                                    
			TCanvas* c = new TCanvas(TString::Format("c_mass_bound=%.3g", val * 1e3) , "3D Histograms", 1200, 600);
    		TLegend* leg = new TLegend(0.7,0.45, 0.9,0.7);	
			TString hname_20 = TString::Format("diagonal %.3gMeV", val*1e3);
    		TH1F* h_20 = new TH1F(hname_20, "", 150, 0, 1000);
			h_20->SetTitle(TString::Format("Comparison of Topologies;p_{t} (MeV/c);Events %.3g [MeV/c]", h_20->GetBinWidth(1)));
			h_20->SetLineColor(kRed);
			h_20->SetLineWidth(1);
    		
			Int_t trans_red = TColor::GetColorTransparent(kRed, 0.3);
			TString hname_40 = TString::Format("parralel %.3gMeV", val*1e3);
    		TH1F* h_40 = new TH1F(hname_40, "", 150, 0, 1000);
			h_40->SetTitle(TString::Format("Comparison of Topologies;p_{t} (MeV/c);Events %.3g [MeV/c]", h_40->GetBinWidth(1)));
			h_40->SetLineColor(kBlack);
			h_40->SetLineWidth(1);

			df_20.Filter([val, phi_mass](const RVecLorCyl p) { // filter for mass cuts around the resonance
				return ((TMath::Abs(p[0].M() - phi_mass) < val) && (TMath::Abs(p[1].M() - phi_mass) < val));
        	}, {"phi_four_momentum"}
    		)
			.Foreach([&h_20](const RVecF& p4) {                                                            
				for (Float_t p : p4){                                                                    
                	h_20->Fill(p *1e3);                                                                          
            	}                                                                                        
        	}, {branch}
    		);

    		df_40.Filter([val, phi_mass](const RVecLorCyl p) { // filter for mass cuts around the resonance
				return ((TMath::Abs(p[0].M() - phi_mass) < val) && (TMath::Abs(p[1].M() - phi_mass) < val));
        	}, {"phi_four_momentum"}
    		)
			.Foreach([&h_40](const RVecF& p4) {                                                            
				for (Float_t p : p4){                                                                    
                	h_40->Fill(p *1e3); 
            	}                                                                                        
        	}, {branch}
    		);
		leg->AddEntry(h_20, "Diagaonal", "L");
		leg->AddEntry(h_40, "Parralel", "L");
		RVecDraw dat = {{h_20, "C"},  {h_40, "C"}, {leg, ""}};
		SaveCanvas(c, dat, TString::Format("media/root_files/phi_phi_reconstruction/phase_space/mass_cuts/%s/compare_topologies/mass_bound = %.3g MeV.root", branch, val*1e3), "RECREATE");
    	
		delete c;
		delete h_20;
		delete h_40;
		};
	
		//delete c;
		//delete leg;
	return 0;
}


