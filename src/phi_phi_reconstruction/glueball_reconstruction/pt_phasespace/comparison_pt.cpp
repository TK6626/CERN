#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TColor.h"
#include "TFitResult.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"

#include "../../../../lib/cut_branch.h"
#include "../../../../lib/computations.h"
#include "../../../../lib/apply_cuts.cpp"
#include "../../../../lib/plotting_params.h"
#include "../../../../lib/admin_utils.h"
#include "../../../../lib/roofithelper.h"

// !./bashing/run_file.sh phi_phi_reconstruction/glueball_reconstruction/pt_phasespace/comparison_pt


/** make cuts on data to 
 * try and remove extraneous kaon backroudn
 * due arising from secodary interactions
 */

Bool_t massWindowCheck(const RVecLorCyl& p, Float_t phi_mass, Float_t mass_bound, Int_t idx1, Int_t idx2);

int main() {

	
	// constants 
	Float_t kaon_mass = m_kaon_char / 1e3;
	Float_t phi_mass = (m_phi) /1e3;
	RVecF mass= {25, 60, 50000};
	mass /= 1e3;
	Int_t nbins = 700;

	RVecStr pairs = {"first", "second", "both"};

	for (const Float_t mass_bound : mass) {
	ROOT::EnableImplicitMT();
		
	for (const TString pairing : pairs) {
		TString file2 = TString::Format("data/phi_phi_reconstruction/uncut_SNR_20.root");
		TString file4 = TString::Format("data/phi_phi_reconstruction/uncut_SNR_40.root");
		RDF df_df_2("tree", file2);
		RDF df_df_4("tree", file4);
	
		TH1F* hist_both = new TH1F("hist", "hist", nbins, 0, 1.7); 
		hist_both->SetTitle(TString::Format("Total p_{t} comparison across Topologies;p_{t} (GeV/c);Events [%.3g GeV/c]", hist_both->GetXaxis()->GetBinWidth(1))); 
		hist_both->SetLineWidth(2);
		
		TH1F* hist_p2 = (TH1F*)hist_both->Clone(TString::Format("#phi Diagonal " + pairing)); hist_p2->SetLineColor(kGreen);
		TH1F* hist_p4 = (TH1F*)hist_both->Clone(TString::Format("#phi Parallel " + pairing)); hist_p4->SetLineColor(kRed);
		TH1F* hist_kaon_2 = (TH1F*)hist_both->Clone(TString::Format("K Diagonal " + pairing));  hist_kaon_2->SetLineColor(kCyan); 
		TH1F* hist_kaon_4 = (TH1F*)hist_both->Clone(TString::Format("K Parallel " + pairing));  hist_kaon_4->SetLineColor(kBlue); 
		gStyle->SetOptStat(111111);


		df_df_2.Filter(
			[phi_mass, mass_bound, pairing](const RVecLorCyl p) {
				if (pairing =="first") {
					return massWindowCheck(p, phi_mass, mass_bound, 0, 1);
				} else if  (pairing == "second") {
					return massWindowCheck(p, phi_mass, mass_bound, 2, 3);
				} else if (pairing == "both") { 
					return (massWindowCheck(p, phi_mass, mass_bound, 0, 1)
					&&
					massWindowCheck(p, phi_mass, mass_bound, 2, 3));
				} else std::cout << "something went wrong" << std::endl; return false;
		}, {"phi_four_momentum"})
			.Foreach(
    		[hist_both, hist_p2, hist_kaon_2, pairing](const RVecLorCyl &pk, const RVecLorCyl &pp) {
        	Float_t P_k = (pk[0] + pk[1] + pk[2] + pk[3]).Pt();
			hist_kaon_2->Fill(P_k);
        	
			if (pairing == "first") {
            	hist_p2->Fill((pp[0] + pp[1]).Pt());
        	} else if (pairing == "second") {
            	hist_p2->Fill((pp[2] + pp[3]).Pt());
        	} else if (pairing == "both") {
            	hist_p2->Fill((pp[0] + pp[1]).Pt());
            	hist_p2->Fill((pp[2] + pp[3]).Pt());
        		hist_kaon_2->Fill(P_k); // account for double filling of phis
        	}
   			}, {"kaon_four_momentum", "phi_four_momentum"}
		);

		df_df_4.Filter(
			[phi_mass, mass_bound, pairing](const RVecLorCyl p) {
				if (pairing =="first") {
					return massWindowCheck(p, phi_mass, mass_bound, 0, 1);
				} else if  (pairing == "second") {
					return massWindowCheck(p, phi_mass, mass_bound, 2, 3);
				} else if (pairing == "both") { 
					return (massWindowCheck(p, phi_mass, mass_bound, 0, 1)
					&&
					massWindowCheck(p, phi_mass, mass_bound, 2, 3));
				} else std::cout << "something went wrong" << std::endl; return false;
		}, {"phi_four_momentum"})

			.Foreach(
			[hist_both, hist_p4, hist_kaon_4, pairing](const RVecLorCyl pk, const RVecLorCyl pp) {
        	
			Float_t P_k = (pk[0] + pk[1] + pk[2] + pk[3]).Pt();
        	hist_kaon_4->Fill(P_k);

        	if (pairing == "first") {
            	hist_p4->Fill((pp[0] + pp[1]).Pt());
        	} else if (pairing == "second") {
            	hist_p4->Fill((pp[2] + pp[3]).Pt());
        	} else if (pairing == "both") {
            	hist_p4->Fill((pp[0] + pp[1]).Pt());
            	hist_p4->Fill((pp[2] + pp[3]).Pt());
        		hist_kaon_4->Fill(P_k); // account for double filling of phis
			}	
			}, {"kaon_four_momentum", "phi_four_momentum"}
		);
		
		
		TCanvas* c = new TCanvas(TString::Format("c " + pairing), TString::Format("c " + pairing));

		TLegend* leg = new TLegend(0.55,0.4, 0.9,0.65);
		//leg1->AddEntry(hist_both, "Diagonal + Parallel");
		leg->AddEntry(hist_p2, "#phi Diagonal");
		leg->AddEntry(hist_p4, "#phi Parallel");
		leg->AddEntry(hist_kaon_2, "K Diagonal");
		leg->AddEntry(hist_kaon_4, "K Parallel");
		RVecDraw dat = {
			{hist_kaon_2, "HIST"},
			{hist_p2, "HIST"},
			{hist_kaon_4, "HIST"},
			{hist_p4, "HIST"},
			{leg, ""}
		};

		if (pairing == "first"){
			SaveCanvas(c, dat, TString::Format("media/root_files/phi_phi_reconstruction/glueball_reconstruction/pt_phasespace/topo_comparison/mass_bound=%.3g_MeV.root", mass_bound*1e3), "RECREATE");
		} else {
			SaveCanvas(c, dat, TString::Format("media/root_files/phi_phi_reconstruction/glueball_reconstruction/pt_phasespace/topo_comparison/mass_bound=%.3g_MeV.root", mass_bound*1e3), "UPDATE");
		}

		delete c;
		delete leg;
		delete hist_both;
		delete hist_p2;
		delete hist_p4;
		delete hist_kaon_2;
		delete hist_kaon_4;
	}	
	}	
	return 0;
}

Bool_t massWindowCheck(const RVecLorCyl& p, Float_t phi_mass, Float_t mass_bound, Int_t idx1, Int_t idx2) {
   		 	return (
				(TMath::Abs(p[idx1].M() - phi_mass) < mass_bound) 
				&&
           		(TMath::Abs(p[idx2].M() - phi_mass) < mass_bound)
				);
		}
