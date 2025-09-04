#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "THStack.h"

#include "../../../../lib/custom_definitions.h"
#include "../../../../lib/cut_branch.h"
#include "../../../../lib/computations.h"
#include "../../../../lib/apply_cuts.cpp"
#include "../../../../lib/plotting_params.h"
#include "../../../../lib/admin_utils.h"

//  !./bashing/run_file.sh reconstruct_all_data/rho_reconstruction/vary_cuts_glueball/dxy_glueball


/** make cuts on data to 
 * see how varying dxy effects the (as so far of yet nonexistant peak at 2220Mev)
 */
int main() {

    // relevent constnats
    const Float_t pion_mass = 1e-3 * m_pi_char; // change to GeV)   

    // Load in aggregated uncut data
    ROOT::EnableImplicitMT();
    TString topology = "40";
    TString file = "data/glueball_mass_reconstruction/data_" + topology + "_SNR.root";
    RDF df_df("tree", file);



    RN df = df_df;

	return 0;
}
