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

		RN df = df_df;
		df.Foreach(
			[&hist] (const RVecLorCyl p4) {	
				
				LorCyl P = p4[0] + p4[1] + p4[2] + p4[3];
				Float_t eta = P.Eta();
				hist->Fill(eta);
			}, {"kaon_four_momentum"}
		);

		double x_min = hist->GetXaxis()->GetXmin();
		double x_max = hist->GetXaxis()->GetXmax();
		RooRealVar x_var("x_var", "eta", x_min, x_max);
	
		RooDataHist data("data", "Dataset from histogram", RooArgList(x_var), hist);
		
		RooRealVar mean_1("mean_1", "mean of gauss 1", -2.5, -3, -2);
		RooRealVar sig1("sigma_1", "width of gauss 1", 0.5, 0.1, 3);
		
		RooRealVar mean_2("mean_2", "mean of gauss 2", 2.5, 2, 3);
		RooRealVar sigma_2("sigma_2", "width of gauss 2", 0.5, 0.1, 3);
		RooRealVar frac("frac", "ratio", 0.5, 0.3, 0.7);
		
		RooAbsPdf* bigauss = BiGaussianPdf("model", "signal gaussian", x_var, mean_1, sigma_1, mean_2, sigma_2, frac);
	
		RooRealVar nsig("nsig", "signal yield", hist->Integral() * 1, 0.9999 * hist->Integral(), hist->Integral());
   		RooAddPdf* model = new RooAddPdf("model", "sig", RooArgList(*bigauss), RooArgList(nsig));
		
		model->fitTo(data);
		calc_chi2(x_var, data, model, nbins);
	
		double binWidth = hist->GetXaxis()->GetBinWidth(1);
		TH1* hModel = model->createHistogram("hModel", x_var, Binning(nbins));
		/*
		TH1* hSignal = createHistogram("hSignal", x_var, Binning(nbins));
		hSignal->SetFillColor(kCyan);
		hSignal->Scale(nsig.getVal());
*/	
		TLegend* leg = new TLegend(0.7,0.45, 0.9, 0.7);
		leg->AddEntry(hist, "Data", "lep");
		leg->AddEntry(hModel, "Model Fit", "l");
//		leg->AddEntry(hSignal, "Gaussian Signal", "f");

		TCanvas* can = new TCanvas("can", "c", 600, 800);
		RVecDraw dat = {
			{hModel, "C"},
			{hist, "EP"}, 
//			{hSignal, "HIST"}, 
//			{hist, "EP"},
			{leg, ""}};
	

		hist->SetTitle(TString::Format("Total #eta Distribution ;Total #eta (GeV/c);Events [%.3g]", hist->GetBinWidth(1)));
		hist->SetMarkerStyle(20);
		
		SaveCanvas(can, dat, TString::Format("media/root_files/phi_phi_reconstruction/other_branch_plots/eta/"+topo+"/mass_cut=%.4g.root", mass_bound*1e3), "RECREATE");
		
		delete can;
		delete hist;
		//delete hModel;
		delete leg;
//		delete model;
//		delete hSignal;
	}
	return 0;
}


