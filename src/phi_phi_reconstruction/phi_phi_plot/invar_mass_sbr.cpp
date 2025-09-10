#include <stdio.h>
#include <iostream>

#include "TH2.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLorentzVector.h"

#include "../../../lib/custom_definitions.h"
#include "../../../lib/cut_branch.h"
#include "../../../lib/computations.h"
#include "../../../lib/apply_cuts.cpp"
#include "../../../lib/plotting_params.h"
#include "../../../lib/admin_utils.h"
#include "../../../lib/roofithelper.h"

#include "RooProdPdf.h"
#include "RooUniform.h"

using namespace RooFit;

// !./bashing/run_file.sh phi_phi_reconstruction/phi_phi_plot/invar_mass_sbr


/** Create simple plot of constructed invariant kaon masses
 * subject to simple mass bound cut
 */

int main(){

	// constants 
	Float_t kaon_mass = m_kaon_char / 1e3;
	Float_t phi_mass = m_phi /1e3;
	Int_t nbins = 100;


	RVecStr topology = {"20", "40"};
	ROOT::EnableImplicitMT();
	for (TString topo : topology) {
		
		TString file = TString::Format("data/phi_phi_reconstruction/uncut_SNR_"	+ topo + ".root");
		RDF df_df("tree", file);
		RN df = df_df;
		TH2F* hist = new TH2F("hist", "Invariant K^{-}K^{+} Mass;M_{1} (MeV/c^{2});M_{2} (MeV/c^{2})", nbins, 985, 1200, nbins, 985,1200);
		
		Float_t dxy_bound = 0.5;
		df.Filter([dxy_bound](const RVecF trk_dxy) {
				for (const Float_t dxy : trk_dxy) {
					if (TMath::Abs(dxy) > dxy_bound) return false;
					}
				return true;
				}, {"trk_dxy"}
				)
		.Foreach(
			[&hist] (const RVecLorCyl p4) {

			hist->Fill(p4[0].M()*1e3, p4[1].M()*1e3);
			hist->Fill(p4[2].M()*1e3, p4[3].M()*1e3);
		}, {"phi_four_momentum"}
		);


		Double_t x_min = hist->GetXaxis()->GetXmin();
		Double_t x_max = hist->GetXaxis()->GetXmax();
		Double_t y_min = hist->GetYaxis()->GetXmin();
		Double_t y_max = hist->GetYaxis()->GetXmax();
		RooRealVar x_var("x_var", "mass", x_min, x_max);
		RooRealVar y_var("x_var", "mass", y_min, y_max);
	

		RooArgList vars(x_var, y_var);
		RooDataHist data("data", "2D data", vars, Import(*hist));
		
	
		RooRealVar meanX("meanX", "Mean X", phi_mass*1e3, 1000, 1040);
		RooRealVar sigmaX("sigmaX", "Sigma X", 20, 0.1, 50);
		RooRealVar meanY("meanY", "Mean Y", phi_mass * 1e3, 1000, 1040);
		RooRealVar sigmaY("sigmaY", "Sigma Y", 20, 0.1, 50);

		RooGaussian gaussX("gaussX", "Gaussian X", x_var, meanX, sigmaX);
		RooGaussian gaussY("gaussY", "Gaussian Y", y_var, meanY, sigmaY);

		RooProdPdf gauss2D("gauss2D", "2D Gaussian", RooArgList(gaussX, gaussY));
		
		RooUniform uniformX("uniformX", "Uniform X", x_var);
		RooUniform uniformY("uniformY", "Uniform Y", y_var);
		RooProdPdf bkg2D("bkg2D", "2D background", RooArgList(uniformX, uniformY));
   
		RooRealVar nsig("nsig", "Signal events", hist->Integral() * 0.5, 0, hist->Integral());
		RooRealVar nbkg("nbkg", "Background events", hist->Integral() * 0.5, 0, hist->Integral());
		

		RooProdPdf* Signal = &gauss2D;
		RooProdPdf* Background = &bkg2D;
		RooAddPdf* model = new RooAddPdf("model", "Signal + Background", RooArgList(*Signal, *Background), RooArgList(nsig, nbkg));
		
		model->fitTo(data);
		calc_chi2(x_var, data, model, nbins);

		hist->SetMarkerStyle(20);

		double binWidth = hist->GetXaxis()->GetBinWidth(1);
		TH2D* hModel = model->createHistogram("hModel", x_var, YVar(y_var));
		hModel->SetLineColor(kBlack);

		TH2* hSignal = Signal->createHistogram("hSignal", x_var, YVar(y_var, Binning(nbins)));
		hSignal->SetFillColor(kCyan);
		hSignal->SetLineColor(kCyan);
		hSignal->Scale(nsig.getVal());
	
		TH1* hBackground = Background->createHistogram("hBackground", x_var, Binning(nbins), YVar(y_var, Binning(nbins)));
		hBackground->Scale(nbkg.getVal());
		hBackground->SetLineColor(kRed);
		hBackground->SetFillColor(kRed);

		TLegend* leg1 = new TLegend(0.7,0.45, 0.9, 0.7);
		TLegend* leg2 = new TLegend(0.7,0.45, 0.9, 0.7);
		leg1->AddEntry(hist, "Data", "lep");
		leg1->AddEntry(hModel, "Model Fit", "l");
		leg2->AddEntry(hSignal, "Gaussian Signal", "f");
		leg2->AddEntry(hBackground, "Background", "f");

		//TLine* l1 = DrawLine(phi_mass * 1e3, 0, hist, kVertical, kRed, 2, kDashed);
		//TLine* l2 = DrawLine(phi_mass * 1e3, 0, hist, kHorizontal, kRed, 2, kDashed);
		TCanvas* c1 = new TCanvas("c1", "c1", 600, 800);
		RVecDraw dat1 = {{hist, "EP"}, {hModel, "surf"}, {leg1, ""}};
		
		TCanvas* c2 = new TCanvas("c2", "c2", 600, 800);
		RVecDraw dat2 = {{hSignal, "surf"}, {hBackground, "surf"}, {leg2, ""}};
		SaveCanvas(c1, dat1, TString::Format("media/root_files/phi_phi_reconstruction/phi_phi_plot/sbr/" + topo + "/Invariant_KK_mass.root"), "RECREATE");
		SaveCanvas(c2, dat2, TString::Format("media/root_files/phi_phi_reconstruction/phi_phi_plot/sbr/" + topo + "/Invariant_KK_mass.root"), "UPDATE");
		

	delete c1;
	delete c2;
	//delete l1;
	//delete l2;
	delete leg1;
	delete leg2;
	delete hist;
	delete Background;
	delete Signal;
	delete hModel;
	delete hSignal;
	delete hBackground;
	}
	return 0;
}

