#include <stdio.h>
#include <stdlib.h>

#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TColor.h"
#include "TString.h"
#include "../../../../lib/custom_definitions.h"
#include "../../../../lib/computations.h"
#include "../../../../lib/plotting_params.h"
#include "../../../../lib/admin_utils.h"

//  !./bashing/run_file.sh reconstruct_all_data/rho_reconstruction/elastic_background/x_y_momentum


/** make cuts on data to 
 * try and remove extraneous kaon backroudn
 * due arising from secodary interactions
 */
int main() {

	RVecChar topologies = {"20", "40"};
	for (const TString& topology: topologies) {	
			
		ROOT::EnableImplicitMT();
		TString file = "data/glueball_mass_reconstruction/data_" + topology + "_SNR.root";
		RDF df_df("tree", file);
		RN df = df_df;

		RVecChar pots = {"a", "b"};
		for (const char* pot: pots) { 
			
			TString px_name = TString::Format("pr_px_%s", pot);
			TString py_name = TString::Format("pr_py_%s", pot);


			// look at distribution of transverse momentum of each of the protons at the roman pots
			TCanvas* c = new TCanvas("c", "c");
			TH1F* h_pt = new TH1F("h_pt", TString::Format("p_{t} pot %s; p_{t} (MeV/c)", pot), 500, 0, 1500);
			TH2F* h_px_py = new TH2F("h_px_py", TString::Format("p_{x} p_{y} Comparison, pot %s", pot), 500, -1200, -1000, 500, 1200, 1000); 

			df.Foreach(
				[&h_pt, &h_px_py](Double_t px, Double_t py) {
				Double_t pt = TMath::Sqrt(px * px + py * py);
		
			h_pt->Fill(pt * 1e3);
			h_px_py->Fill(px * 1e3, py * 1e3);
			}, {px_name.Data(), py_name.Data()}
			);

			h_pt->GetYaxis()->SetTitle(Form("Events [%.2g MeV/c]", h_pt->GetBinWidth(1)));
			h_px_py->GetYaxis()->SetTitle(Form(" p_{y} (MeV) [%.2g MeV/c]", h_px_py->GetYaxis()->GetBinWidth(1)));
			h_px_py->GetXaxis()->SetTitle(Form(" p_{x} (MeV) [%.2g MeV/c]", h_px_py->GetXaxis()->GetBinWidth(1)));
			

			SaveCanvas(c, h_pt, "E1", TString::Format("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/elastic_background/px_and_py/pt_"+ topology + "_pot_%s.root", pot), "RECREATE");
		SaveCanvas(c, h_px_py, "COLZ", TString::Format("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/elastic_background/px_and_py/px_py_" + topology + "_pot_%s.root", pot), "RECREATE");
		h_pt->Delete();
		h_px_py->Delete();

		}
	}
	return 0;
}

