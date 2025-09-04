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

//  !./bashing/run_file.sh reconstruct_all_data/rho_reconstruction/elastic_background/pr_pz


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
		TH1F* pr_z_sum = new TH1F("pr_z_sum", "pz Sum Across Pots;Events;pr_z (GeV)", 500, 0.3, 2);
		TH1F* pr_z_diff = new TH1F("pr_z_diff", "pz Difference Across Pots;Events;pr_z (GeV)", 500, -1.5, 1.5);
		TH2F* pr_z_comp = new TH2F("pr_z_comp", "pz Comparison", 500, 0.1, 1.2, 500, 0.1, 1.2); 

		df.Foreach(
			[&pr_z_sum, &pr_z_diff, & pr_z_comp](Double_t pz_a, Double_t pz_b) {
			pz_a -= 6500;
			pz_b += 6500;
			Double_t pr_z_a = pz_b;
			Double_t pr_z_b = pz_b;

			pr_z_sum->Fill(pr_z_a + pr_z_b);
			pr_z_diff->Fill(pr_z_a - pr_z_b);
			pr_z_comp->Fill(pr_z_a, pr_z_b);

			}, {"pr_pz_a", "pr_pz_b"});

		pr_z_sum->GetYaxis()->SetTitle(Form("Events [%.2g GeV/c^{2}]", pr_z_diff->GetBinWidth(1)));
		pr_z_diff->GetYaxis()->SetTitle(Form("Events [%.2g GeV/c^{2}]", pr_z_sum->GetBinWidth(1)));
		pr_z_comp->GetYaxis()->SetTitle(Form(" pr_z_{t} Pot A (GeV) [%.2g GeV/c^{2}]", pr_z_comp->GetYaxis()->GetBinWidth(1)));
		pr_z_comp->GetXaxis()->SetTitle(Form(" pr_z_{t} Pot B (GeV) [%.2g GeV/c^{2}]", pr_z_comp->GetXaxis()->GetBinWidth(1)));


		SaveCanvas(c, pr_z_sum, "E1", TString::Format("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/elastic_background/transverse_momentum/pr_z_" + topology+ "_sum.root"), "RECREATE");
		SaveCanvas(c, pr_z_diff, "E1", TString::Format("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/elastic_background/transverse_momentum/pr_z_" + topology+ "_difference.root"), "RECREATE");
		SaveCanvas(c, pr_z_comp, "COLZ", TString::Format("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/elastic_background/transverse_momentum/pr_z_" + topology+ "_comparison.root"), "RECREATE");
		}
	return 0;
}

