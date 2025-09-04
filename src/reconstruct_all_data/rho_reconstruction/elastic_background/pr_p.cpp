#include <stdio.h>
#include <stdlib.h>

#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TColor.h"

#include "../../../../lib/custom_definitions.h"
#include "../../../../lib/computations.h"
#include "../../../../lib/plotting_params.h"
#include "../../../../lib/admin_utils.h"

//  !./bashing/run_file.sh reconstruct_all_data/rho_reconstruction/elastic_background/pr_p


/** make cuts on data to 
 * try and remove extraneous kaon backroudn
 * due arising from secodary interactions
 */
int main() {

	
	RVecStr topologies = {"20", "40"};
	for (const TString& topology: topologies) {	
		ROOT::EnableImplicitMT();
		TString file = "data/glueball_mass_reconstruction/data_" + topology + "_SNR.root";
		RDF df_df("tree", file);
		std::cout << df_df.Describe() << std::endl;
		RN df = df_df;
		
		// look at distribution of transverse momentum of each of the protons at the roman pots
		TCanvas* c = new TCanvas("c", "Title");
		TH1F* p_sum = new TH1F("p_sum", "p Sum Across Pots;Events;p_{t} (GeV)", 500, 0.3, 2);
		TH1F* p_diff = new TH1F("p_diff", "p Difference Across Pots;Events;p_{t} (GeV)", 500, -1.5, 1.5);
		TH2F* p_comp = new TH2F("p_comp", "p Comparison", 500, 0.1, 1.2, 500, 0.1, 1.2); 

		df.Foreach(
			[&p_sum, &p_diff, & p_comp](Double_t px_a, Double_t py_a, Double_t pz_a, Double_t px_b, Double_t py_b, Double_t pz_b) {
			pz_a -= 6500;
			pz_b += 6500;
			Double_t p_a = TMath::Sqrt(px_a * px_a + py_a * py_a + pz_a * pz_a);
			Double_t p_b = TMath::Sqrt(px_b * px_b + py_b * py_b + pz_b * pz_b);

			p_sum->Fill(p_a + p_b);
			p_diff->Fill(p_a - p_b);
			p_comp->Fill(p_a, p_b);

			}, {"pr_px_a", "pr_py_a", "pr_pz_a", "pr_px_b", "pr_py_b", "pr_pz_b"});

		p_sum->GetYaxis()->SetTitle(Form("Events [%.2g GeV/c^{2}]", p_diff->GetBinWidth(1)));
		p_diff->GetYaxis()->SetTitle(Form("Events [%.2g GeV/c^{2}]", p_sum->GetBinWidth(1)));
		p_comp->GetYaxis()->SetTitle(Form(" p_{t} Pot A (GeV) [%.2g GeV/c^{2}]", p_comp->GetYaxis()->GetBinWidth(1)));
		p_comp->GetXaxis()->SetTitle(Form(" p_{t} Pot B (GeV) [%.2g GeV/c^{2}]", p_comp->GetXaxis()->GetBinWidth(1)));


		SaveCanvas(c, p_sum, "E1", TString::Format("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/elastic_background/transverse_momentum/p_" + topology+ "_sum.root"), "RECREATE");
		SaveCanvas(c, p_diff, "E1", TString::Format("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/elastic_background/transverse_momentum/p_" + topology+ "_difference.root"), "RECREATE");
		SaveCanvas(c, p_comp, "COLZ", TString::Format("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/elastic_background/transverse_momentum/p_" + topology+ "_comparison.root"), "RECREATE");
		}
	return 0;
}

