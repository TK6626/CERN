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

#include "../../../lib/cut_branch.h"
#include "../../../lib/computations.h"
#include "../../../lib/apply_cuts.cpp"
#include "../../../lib/plotting_params.h"
#include "../../../lib/admin_utils.h"

// !./bashing/run_file.sh phi_phi_reconstruction/protons/pt_topo_comparison


/** make cuts on data to 
 * try and remove extraneous kaon backroudn
 * due arising from secodary interactions
 */

Bool_t massWindowCheck(const RVecLorCyl& p, Float_t phi_mass, Float_t mass_bound, Int_t idx1, Int_t idx2);

int main() {

	
	// constants 
	Float_t kaon_mass = m_kaon_char / 1e3;
	Float_t phi_mass = (m_phi + 50) /1e3;
	RVecF mass= {15, 20, 25};
	mass /= 1e3;
	Int_t nbins = 150;
	
	Float_t pt_low = 0.05;
	Float_t pt_high = 0.15;

	RVecStr pairs = {"first", "second", "both"};
	for (const Float_t mass_bound : mass) {
	ROOT::EnableImplicitMT();
	gStyle->SetOptStat(0);
		
	for (const TString pairing : pairs) {
		TString file2 = TString::Format("data/phi_phi_reconstruction/uncut_SNR_20.root");
		TString file4 = TString::Format("data/phi_phi_reconstruction/uncut_SNR_40.root");
		RDF df_df_2("tree", file2);
		RDF df_df_4("tree", file4);
	
		TH1F* hist_both = new TH1F("hist 1", "hist 1", nbins, 0, 1.7); 
		hist_both->SetTitle(TString::Format("Total p_t Comparison Across Topologies;p_{t} (GeV/c);Events [%.3g GeV/c]", hist_both->GetXaxis()->GetBinWidth(1))); 
		hist_both->SetLineWidth(3);
		
		TH1F* hist_2 = (TH1F*)hist_both->Clone(TString::Format("Proton Diagonal "+pairing)); hist_2->SetLineColor(kRed);
		TH1F* hist_4 = (TH1F*)hist_both->Clone(TString::Format("Proton Parallel "+pairing)); hist_4->SetLineColor(kBlue);
		
		Int_t n = 220;
		TH1F* hist_phi = new TH1F("hist mass", "hist mass", n, 1.8, 5); 
		hist_phi->SetTitle(TString::Format("X mass %s;M (GeV/c);Events [%.3g GeV/c]", pairing.Data(), hist_both->GetXaxis()->GetBinWidth(1))); 
		hist_phi->SetLineWidth(3);
		
		TH1F* hist_phi_2 = (TH1F*)hist_phi->Clone(TString::Format("X mass Diagonal "+pairing)); hist_phi_2->SetLineColor(kRed);
		TH1F* hist_phi_4 = (TH1F*)hist_phi->Clone(TString::Format("X mass Parallel "+pairing)); hist_phi_4->SetLineColor(kBlue);
		//gStyle->SetOptStat(111111);
		gStyle->SetOptStat(0);

		df_df_2.Filter(
			[phi_mass, mass_bound, pairing](const RVecLorCyl p) {
				if (pairing =="first") {
					return massWindowCheck(p, phi_mass, mass_bound, 0, 1);
				} else if  (pairing == "second") {
					return massWindowCheck(p, phi_mass, mass_bound, 2, 3);
				} else if (pairing == "both") { 
					return (massWindowCheck(p, phi_mass, mass_bound, 0, 1)
					||
					massWindowCheck(p, phi_mass, mass_bound, 2, 3));
				} else std::cout << "something went wrong" << std::endl; return false;
		}, {"phi_four_momentum"})
			.Foreach(
    		[hist_2, hist_phi_2, pairing, pt_low, pt_high](const Double_t px_a, const Double_t py_a, const Double_t px_b, const Double_t py_b, const RVecLorCyl p) {	
				Double_t px = px_a + px_b;
				Double_t py = py_a + py_b;
				Double_t pt = TMath::Sqrt(px*px + py*py);
				
				hist_2->Fill(pt);
			
				if ((pt_low < pt) && (pt < pt_high)) { // histogram of glueball mass depending what the associated proton transvers momentum was
					if (pairing =="first") {
					hist_phi_2->Fill((p[0] + p[1]).M());
					} else if  (pairing == "second") {
					hist_phi_2->Fill((p[2] + p[3]).M());
					} else if (pairing == "both") { 
					hist_phi_2->Fill((p[0] + p[1]).M());
					hist_phi_2->Fill((p[2] + p[3]).M());
				}
				}
			}, {"pr_px_a", "pr_py_a", "pr_px_b", "pr_py_b", "phi_four_momentum"}
		);


		df_df_4.Filter(
			[phi_mass, mass_bound, pairing](const RVecLorCyl p) {
				if (pairing =="first") {
					return massWindowCheck(p, phi_mass, mass_bound, 0, 1);
				} else if  (pairing == "second") {
					return massWindowCheck(p, phi_mass, mass_bound, 2, 3);
				} else if (pairing == "both") { 
					return (massWindowCheck(p, phi_mass, mass_bound, 0, 1)
					||
					massWindowCheck(p, phi_mass, mass_bound, 2, 3));
				} else std::cout << "something went wrong" << std::endl; return false;
		}, {"phi_four_momentum"})
		.Foreach(
    		[hist_4, hist_phi_4, pairing, pt_low, pt_high](const Double_t px_a, const Double_t py_a, const Double_t px_b, const Double_t py_b, const RVecLorCyl p) {	
				Double_t px = px_a + px_b;
				Double_t py = py_a + py_b;
				Double_t pt = TMath::Sqrt(px*px + py*py);
				hist_4->Fill(pt);


				if ((pt_low < pt) && (pt < pt_high)) {
					if (pairing =="first") {
					hist_phi_4->Fill((p[0] + p[1]).M());
					} else if  (pairing == "second") {
					hist_phi_4->Fill((p[2] + p[3]).M());
					} else if (pairing == "both") { 
					hist_phi_4->Fill((p[0] + p[1]).M());
					hist_phi_4->Fill((p[2] + p[3]).M());
					}
				}
			}, {"pr_px_a", "pr_py_a", "pr_px_b", "pr_py_b", "phi_four_momentum"}
		);
		
		
		TCanvas* c = new TCanvas(TString::Format("c " + pairing), TString::Format("c " + pairing));
		TCanvas* c_phi = new TCanvas(TString::Format("c mass " + pairing), TString::Format("c " + pairing));

		TLegend* leg = new TLegend(0.7,0.7, 0.9,0.9);
		//leg1->AddEntry(hist_both, "Diagonal + Parallel");
		leg->AddEntry(hist_2, "Diagonal");
		leg->AddEntry(hist_4, "Parallel");
		RVecDraw dat = {
			{hist_2, "HIST"},
			{hist_4, "HIST"},
			{leg, ""}
		};
		
		TLine* l1 = DrawLine(2*phi_mass, 0, hist_phi_2, kVertical, kBlack, 3, kDashed);
		TLegend* leg2 = new TLegend(0.7,0.7, 0.9,0.9);
		leg2->AddEntry(hist_phi_2, "Diagonal");
		leg2->AddEntry(hist_phi_4, "Parallel");
		leg2->AddEntry(l1, TString::Format("M=%.2fGeV/c",2*phi_mass));
		RVecDraw dat2 = {
			{hist_phi_2, "HIST"},
			{hist_phi_4, "HIST"},
			{l1, ""},
			{leg2, ""}
		};

		if (pairing == "first"){
			SaveCanvas(c, dat, TString::Format("media/root_files/phi_phi_reconstruction/protons/pt/topo_comparison/pt_bounds=%.3g-%.3g/mass_bound=%.3g_MeV.root",pt_low, pt_high, mass_bound*1e3), "RECREATE");
			SaveCanvas(c_phi, dat2, TString::Format("media/root_files/phi_phi_reconstruction/protons/pt/topo_comparison/pt_bounds=%.3g-%.3g/mass_bound=%.3g_MeV.root",pt_low, pt_high, mass_bound*1e3), "UPDATE");
		} else {
			SaveCanvas(c, dat, TString::Format("media/root_files/phi_phi_reconstruction/protons/pt/topo_comparison/pt_bounds=%.3g-%.3g/mass_bound=%.3g_MeV.root",pt_low, pt_high, mass_bound*1e3), "UPDATE");
			SaveCanvas(c_phi, dat2, TString::Format("media/root_files/phi_phi_reconstruction/protons/pt/topo_comparison/pt_bounds=%.3g-%.3g/mass_bound=%.3g_MeV.root",pt_low, pt_high, mass_bound*1e3), "UPDATE");
		}
		

		delete c;
		delete leg;
		delete c_phi;
		delete leg2;
		delete hist_both;
		delete hist_2;
		delete hist_4;
		delete hist_phi;
		delete hist_phi_2;
		delete hist_phi_4;
	}	
	}	
	return 0;
}

Bool_t massWindowCheck(const RVecLorCyl& p, Float_t phi_mass, Float_t m, Int_t idx1, Int_t idx2) {
   		 	return (
				(TMath::Abs(p[idx1].M() - phi_mass) < m) 
				&&
           		(TMath::Abs(p[idx2].M() - phi_mass) < m)
				);
		}

