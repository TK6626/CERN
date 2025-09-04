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

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooArgList.h"
#include "RooConstVar.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"

#include "../../../../lib/cut_branch.h"
#include "../../../../lib/computations.h"
#include "../../../../lib/apply_cuts.cpp"
#include "../../../../lib/plotting_params.h"
#include "../../../../lib/admin_utils.h"
#include "../../../../lib/roofithelper.h"

// !./bashing/run_file.sh phi_phi_reconstruction/glueball_reconstruction/full_glueball/vary_pt_window


/** make cuts on data to 
 * try and remove extraneous kaon backroudn
 * due arising from secodary interactions
 */

using namespace RooFit;
int main() {

	
	// constants 
	Float_t kaon_mass = m_kaon_char / 1e3;
	Float_t phi_mass = m_phi /1e3;
	Float_t mass_bound= 20;
	Int_t nbins = 70;
	RVecF pt_centre = {0. 83, 0.85, 0.87, 0.89, 0.91, 0.93};

	Float_t sf = 1e-3; // introduce a scale factor to make fitting easier convert from MeV to whateve
	for (const Float_t pt : pt_centre) {

	RVecStr topology = {"20", "40"};
	ROOT::EnableImplicitMT();
	for (TString topo : topology) {
		
		TString file = TString::Format("data/phi_phi_reconstruction/uncut_SNR_"	+ topo + ".root");
		RDF df_df("tree", file);
	
		TH1F* hist = new TH1F("hist", "hist", nbins, 2000*sf, 3000*sf); // make fitting easier by converting to TeV 

		Float_t low_pt = pt - 0.02; // GeV
		Float_t high_pt = pt + 0.01; // GeV
		Float_t low_eta = 3; // GeV
		Float_t high_eta = 9.0; // GeV

		RN df = df_df.Filter(
				[low_pt, high_pt](const RVecF pt) {
					return ((low_pt < Sum(pt)) && (Sum(pt) < high_pt));
				}, {"trk_pt"})
		.Filter([low_eta, high_eta] (const RVecF eta) {
					return ((low_eta < Sum(eta)) && (Sum(eta) < high_eta));
				},{"trk_eta"})
		.Filter([&hist, mass_bound, phi_mass](const RVecLorCyl p) {
			return (
				(TMath::Abs((p[0].M() - phi_mass)) < mass_bound)  
				&& 
				(TMath::Abs((p[1].M() - phi_mass)) < mass_bound)
			);}, {"phi_four_momentum"}
		);
	
		df.Foreach([&hist, sf, phi_mass](const RVecLorCyl p) {
			{
			hist->Fill((p[0] + p[1]).M() *sf*1e3);
			} // convert to TeV FROM GeV!!!!
		}, {"phi_four_momentum"});
	
		hist->SetTitle(TString::Format("Invariant X Mass;M (MeV/c^{2});Events [%.2g MeV]", hist->GetXaxis()->GetBinWidth(1)));
		//hist->Scale(1.0 / hist->Integral()); // to make fitting easier for program
		Double_t x_min = hist->GetXaxis()->GetXmin();
		Double_t x_max = hist->GetXaxis()->GetXmax();
		RooRealVar x_var("x_var", "mass", x_min, x_max);
	
		RooDataHist data("data", "Dataset from histogram", RooArgList(x_var), hist);
		
		RooRealVar mean("mean", "mean of gauss", 2220 * sf, 2210 *sf, 2240*sf); //Mev -> TeV
		RooRealVar sigma("sigma", "width of gauss", 10*sf, sf, 20*sf);
		RooGaussian gauss("gauss", "signal gaussian", x_var, mean, sigma);
	
		RooRealVar x_var0("x_var0", "endpoint", x_max, x_min, x_max*2);
		RooRealVar p("p", "power", 1, 1e-5, 20);
		RooRealVar c("c", "curvature", 5, 0.001, 40);
		RooRealVar shift("shift", "pt offset", 2000, 1800, 2100);
		RooAbsPdf* argus = ShiftedArgusPdf("argus", "shifted ARGUS", x_var, x_var0, c, p, shift);

		RooRealVar x_0("x_0", "threshold", x_min, x_min*0.5, x_min*1.2);
		RooRealVar p0("p0", "power", 17, 8, 25);
		RooRealVar a0("a0", "a0", 0, -20, 20);
		RooRealVar b0("b0", "b0", 0, -20, 40);
		RooRealVar c0("c0", "c0", -80, -80, 1);

		RooAbsPdf* Background = ThresholdBackgroundPdf("threshbkg", "threshbkg", x_var, x_0, p0, a0, b0, c0);
	
	// Combine signal and background
		RooRealVar nsig("nsig", "signal yield", hist->Integral() * 0.5, 0, hist->Integral());
		RooRealVar nbkg("nbkg", "background yield", hist->Integral() * 0.5, 0, hist->Integral());
   
		//RooAbsPdf* Background = &argus;
		RooAbsPdf* Signal = &gauss;
		//RooAddPdf* model = new RooAddPdf("model", "sig+bg", RooArgList(*Background), RooArgList(nbkg));
		RooAddPdf* model = new RooAddPdf("model", "sig+bg", RooArgList(*Signal, *Background), RooArgList(nsig,nbkg));
	
		model->fitTo(data);
		calc_chi2(x_var, data, model, nbins);
	
		double binWidth = hist->GetXaxis()->GetBinWidth(1);
		TH1* hModel = model->createHistogram("hModel", x_var, Binning(nbins));
		hist->SetMarkerStyle(20);

		TH1* hSignal = Signal->createHistogram("hSignal", x_var, Binning(nbins));
		hSignal->SetFillColor(kCyan);
		hSignal->Scale(nsig.getVal());
	
		TH1* hBackground = Background->createHistogram("hBackground", x_var, Binning(nbins));
		hBackground->Scale(nbkg.getVal());
		hBackground->SetFillColor(kRed);

		TLegend* leg = new TLegend(0.7,0.45, 0.9, 0.7);
		leg->AddEntry(hist, "Data", "lep");
		leg->AddEntry(hModel, "Model Fit", "l");
		leg->AddEntry(hSignal, "Gaussian Signal", "f");
		leg->AddEntry(hBackground, "Background", "f");

		TLine* l1 = DrawLine(2220, 0, hist, kVertical); 
		TCanvas* can = new TCanvas("can", "c", 600, 800);
		RVecDraw dat = {{hist, "EP"}, {hModel, "C"}, {hBackground, "HIST"}, {hSignal, "HIST"}, {hist, "EP"}, {leg, ""}};
	 		
		
		SaveCanvas(can, dat, TString::Format("media/root_files/phi_phi_reconstruction/glueball_reconstruction/full_glueball/vary_window/vary_pt_window/"+ topo +"/mass_bound=20MeV_pt_centre=%.4g_MeV.root", pt*1e3), "RECREATE");
	
		delete can;
		delete l1;
		delete hist;
		delete hModel;
		delete hSignal;
		delete hBackground;	
		
	}	
	}
	return 0;
}


