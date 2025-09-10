#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"


#include "../../../lib/roofithelper.h"
#include "../../../lib/computations.h"
#include "../../../lib/admin_utils.h"
#include "../../../lib/plotting_params.h"


// !./bashing/run_file.sh phi_phi_reconstruction/other_branch_plots/eta
using namespace RooFit;

int main() {

	SetPlotStyle();
	// set up transparant colours for use later
	Float_t kaon_mass = m_kaon_char / 1e3;
	Float_t phi_mass = m_phi / 1e3;
	Float_t mass_bound = 10 / 1e3; 

	Int_t nbins = 100;
	RVecStr topology = {"20", "40"};
	ROOT::EnableImplicitMT();
	
		
	for (TString topo : topology) {

		TString file = TString::Format("data/phi_phi_reconstruction/uncut_SNR_"	+ topo + ".root");
		RDF df_df("tree", file);
		TH1F* hist = new TH1F("eta", "", nbins, -7, 7); 

		RN df1 = df_df;
		RN df2 = df_df;
		
		df.Foreach(
			[&hist] (const RVecLorCyl p4) {	
				
				LorCyl P = p4[0] + p4[1] + p4[2] + p4[3];
				Float_t eta = P.Eta();
				hist->Fill(eta);
			}, {"phi_four_momentum"}
		);

		TCanvas* can1 = new TCanvas("can1", "c", 600, 800);
		RVecDraw dat1 = {
			{hist, "EP"}, 
			{leg, ""}};
	

		hist->SetTitle(TString::Format("Total #eta Distribution ;Total #eta (GeV/c);Events [%.3g]", hist->GetBinWidth(1)));
		hist->SetMarkerStyle(20);
		
		SaveCanvas(can, dat, TString::Format("media/root_files/phi_phi_reconstruction/other_branch_plots/eta/"+topo+"/mass_cut=%.4g.root", mass_bound*1e3), "RECREATE");
		
		delete can;
		delete hist;
		delete leg;
	}
	return 0;
}



