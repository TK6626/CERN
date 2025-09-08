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

// !./bashing/run_file.sh phi_phi_reconstruction/glueball_reconstruction/eta_phasespace/glueball


/** make cuts on data to 
 * try and remove extraneous kaon backroudn
 * due arising from secodary interactions
 */

using namespace RooFit;
int main() {

	
	// constants 
	Float_t kaon_mass = m_kaon_char / 1e3;
	Float_t phi_mass = m_phi /1e3;
	RVecF mass= {25, 20, 15, 10, 5};


	Float_t high_eta = 10; // GeV
	RVecF eta_low_bounds = {0, 0.5, 1, 1.5, 2, 2.5, 3};
	mass /= 1e3;
	Int_t nbins = 150;

	for (const Float_t low_eta : eta_low_bounds){
	for (const Float_t mass_bound : mass) {

	RVecStr topology = {"20", "40"};
	ROOT::EnableImplicitMT();
	for (TString topo : topology) {
		
		TString file = TString::Format("data/phi_phi_reconstruction/uncut_SNR_"	+ topo + ".root");
		RDF df_df("tree", file);
	
		TH1F* hist = new TH1F("hist", "hist", nbins, 2000, 3000);


		RN df = df_df.Filter(
				[&low_eta, &high_eta](const RVecLorCyl p) {
					LorCyl P = p[0] + p[1] + p[2] + p[3]; 
					Float_t eta = P.Eta();
				
				return (
						(low_eta < eta) && (eta < high_eta)
						||
						(-low_eta > eta) && (eta > -high_eta));
				}, {"kaon_four_momentum"});

		df.Foreach(
			[&hist, mass_bound, phi_mass](const RVecLorCyl p) {
			if (
				(TMath::Abs((p[0].M() - phi_mass)) < mass_bound)  
				&& 
				(TMath::Abs((p[1].M() - phi_mass)) < mass_bound)
			) {hist->Fill((p[0] + p[1]).M() * 1e3);}
			
			}, {"phi_four_momentum"}
		);
		hist->SetTitle(TString::Format("Invariant X Mass;M (MeV/c^{2});Events [%.2g MeV]", hist->GetXaxis()->GetBinWidth(1)));
	
		double x_min = hist->GetXaxis()->GetXmin();
		double x_max = hist->GetXaxis()->GetXmax();
		RooRealVar x_var("x_var", "eta", x_min, x_max);
	

		RooDataHist data("data", "Dataset from histogram", RooArgList(x_var), hist);
		
		RooRealVar mean("mean", "mean of gauss", 2220, 2210, 2230);
		RooRealVar sigma("sigma", "width of gauss", 10, 0.1, 50);
		RooGaussian gauss("gauss", "signal gaussian", x_var, mean, sigma);
	
		RooRealVar x_var0("x_var0", "endpoint", x_max, x_min, x_max*2);
		RooRealVar p("p", "power", 1, 1e-5, 20);
		RooRealVar c("c", "curvature", 5, 0.001, 40);
		RooRealVar shift("shift", "eta offset", 2040, 2000, 2100);
		RooAbsPdf* argus = ShiftedArgusPdf("argus", "shifted ARGUS", x_var, x_var0, c, p, shift);
	
	// Combine signal and background
		RooRealVar nsig("nsig", "signal yield", hist->Integral() * 0.5, 0, hist->Integral());
		RooRealVar nbkg("nbkg", "background yield", hist->Integral() * 0.5, 0, hist->Integral());
   
		RooAbsPdf* Background = argus;
		RooAbsPdf* Signal = &gauss;
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
		RVecDraw dat = {{hist, "EP"}, {hModel, "C"}, {hSignal, "HIST"}, {hBackground, "HIST"}, {hist, "EP"}, {leg, ""}};
		SaveCanvas(can, dat, TString::Format("media/root_files/phi_phi_reconstruction/glueball_reconstruction/eta_phasespace/"+ topo +"/%.3g/mass_bound=%.3g_MeV.root", low_eta, mass_bound*1e3), "RECREATE");
		
		delete can;
		delete l1;
		delete hist;
		delete hModel;
		delete hSignal;
		delete hBackground;	
	}	
	}
	}
	return 0;
}

