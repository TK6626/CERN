#include <iostream>

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"

#include "../../../lib/custom_definitions.h"
#include "../../../lib/admin_utils.h"

// !bashing/run_file.sh phi_phi_reconstruction/phase_space/mass_bound_compare_fit
using namespace RooFit;

int main() {
   
	const char* branch = "trk_pt";
	//RVecF mass_cuts = {50, 30, 25, 20, 15, 10, 5};
	Float_t val = 5; 

	TFile *f = TFile::Open("media/root_files/phi_phi_reconstruction/phase_space/mass_cuts/trk_pt/compare_topologies/mass_bound = 5 MeV.root");
    TCanvas* can = (TCanvas*) f->Get("c_mass_bound_20");  // change name to your histogram
 	canls();
	
return 0;
}
