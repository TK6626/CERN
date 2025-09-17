#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFitResult.h"


#include "../../../lib/computations.h"
#include "../../../lib/admin_utils.h"
#include "../../../lib/plotting_params.h"
#include "../../../lib/roofithelper.h"

// !./bashing/run_file.sh phi_phi_reconstruction/phi_phi_plot/projections

using namespace RooFit;
int main() {

	SetPlotStyle();
	
	// set up transparant colours for use later
	
	Float_t kaon_mass = m_kaon_char;
	Float_t phi_mass = m_phi;
	Float_t mass_bound = 10; 
	Float_t lower_bound = phi_mass - mass_bound; 
	Float_t upper_bound = phi_mass + mass_bound;


	RVecStr projection_axis = {"x"};
	RVecStr topology = {"40"};
	ROOT::EnableImplicitMT();
	
		
	for (TString topo : topology) {
	
		TFile* file = new TFile(TString::Format("media/root_files/phi_phi_reconstruction/phi_phi_plot/dxy/" + topo + "/0.10_cut.root"), "READ");
		TH2F* phi_phi_h = (TH2F*)file->Get("h");
		for (TString ax : projection_axis) {

			TCanvas* c = new TCanvas("c", "c", 600, 800);
			TH1D *h = nullptr;
			
			if (ax == "y"){
				Int_t lower_x_bin = phi_phi_h->GetXaxis()->FindFixBin(lower_bound);
				Int_t upper_x_bin = phi_phi_h->GetXaxis()->FindFixBin(upper_bound);
				h = phi_phi_h->ProjectionY("h_proj", lower_x_bin, upper_x_bin);
			}
			if (ax == "x") {
				Int_t lower_y_bin = phi_phi_h->GetYaxis()->FindFixBin(lower_bound);
				Int_t upper_y_bin = phi_phi_h->GetYaxis()->FindFixBin(upper_bound);
				h = phi_phi_h->ProjectionX("h_pro", lower_y_bin, upper_y_bin);
			}
			
			h->SetTitle(TString::Format(";X Projection M_{K^{+}K^{-}} (MeV/c^{2});Events [%.2g MeV/c^{2}]", h->GetBinWidth(1)));
			c->SetTopMargin(0.02);
			c->SetRightMargin(0.05);
			h->GetXaxis()->SetTitleOffset(1.2);
			gStyle->SetOptStat(00000001); 


			double x_min = h->GetXaxis()->GetXmin();
			double x_max = h->GetXaxis()->GetXmax();
			RooRealVar x("x", "mass", x_min, x_max);
			Float_t sf = 1e-3;	
			Int_t nbins = 300;		

			RooDataHist data("data", "Dataset from histogram", RooArgList(x), h);

			RooRealVar mean1("mean1", "mean of gauss1", 1019.7, 1019, 1021 );
			RooRealVar sigma1("sigma1", "width of gauss1", 7, 4, 10);
			RooGaussian* gauss1 = new RooGaussian("gauss1", "gauss1", x, mean1, sigma1);

			RooRealVar mean2("mean2", "mean of gauss2", 1240, 1200, 1300 );
			RooRealVar sigma2("sigma2", "width of gauss2", 40, 20, 60 );
			RooGaussian* gauss2 = new RooGaussian("gauss2", "gauss2", x, mean2, sigma2);
			
			RooRealVar x_0("x_0", "threshold", 988, 980,  995);
			RooRealVar p0("p0", "power", 0.21, 0.1, 0.5);
			RooRealVar a0("a0", "a0", 0, -1, 1);
			RooRealVar b0("b0", "b0", 0, -1e-1, 1e-1);
			a0.setConstant(true);
			b0.setConstant(true);
			RooRealVar c0("c0", "c0", -7e-3, -2, 0.0);
			RooAbsPdf* thresh = ThresholdBackgroundPdf("threshbkg", "threshbkg", x, x_0, p0, a0, b0, c0);

			RooAbsPdf* Background = thresh;

			Float_t h_int = h->Integral();
			RooRealVar nsig1("nsig1", "yield gauss1", h_int * 0.15, 0.03*h_int, 0.3*h_int);
			RooRealVar nsig2("nsig2", "yield gauss2", h_int * 0.10, 0, 0.35*h_int);
			RooRealVar nbkg("nbkg", "background yield", h_int * 0.85, h_int * 0.5, 0.95*h_int);

			RooAddPdf* model = new RooAddPdf("model", "sig+bg",
                                 RooArgList(*gauss1, *gauss2, *Background),
                                 RooArgList(nsig1, nsig2, nbkg));

			model->fitTo(data);
			calc_chi2(x, data, model, nbins);
			
			Float_t binwidth = h->GetXaxis()->GetBinWidth(1);
			h->SetMarkerStyle(20);
			h->SetLineColor(kBlack);

			TH1* hModel = model->createHistogram("hModel", x, Binning(nbins));
			hModel->SetLineWidth(4);
			hModel->SetLineColor(kCyan);
			
			TH1* hphi = gauss1->createHistogram("hgauss1", x, Binning(nbins));
			hphi->Scale(nsig1.getVal());
			hphi->SetFillColor(color_scheme[1]);
			
			TH1* hSignal = gauss2->createHistogram("hSignal", x, Binning(nbins));
			hSignal->Scale(nsig2.getVal());
			hSignal->SetFillColor(color_scheme[2]);
			
			TH1* hBackground = Background->createHistogram("hBackground", x, Binning(nbins));
			hBackground->Scale(nbkg.getVal());
			hBackground->SetFillColor(color_scheme[3]);
			

			Float_t mu = mean1.getVal();
			Float_t sig = sigma1.getVal();
			std::cout << mu << "\n" << sig << std::endl;

			Float_t n = 1;
			Float_t int_data = hphi->Integral(hphi->FindBin(mu - n * sig), hphi->FindBin(mu + n*sig));
			Float_t int_model = hModel->Integral(hModel->FindBin(mu -n * sig), hModel->FindBin(mu + n*sig));
			std::cout << "ratio = "<< int_data / int_model << std::endl;


			TLegend* l = new TLegend(0.65, 0.55, 0.90, 0.70); // bottom left, top right
			// l->SetHeader("Invariant Mass Fits", "C");
			l->AddEntry(h, "Data", "lep");
			l->AddEntry(hModel, "Model", "l");
			l->AddEntry(hphi, "#phi Resonance", "f");
			l->AddEntry(hSignal, "Unknown Resonance", "f");
			l->AddEntry(hBackground, "Background", "f");
			l->SetTextSize(0.035);
			
			RVecDraw dat = {{h, "E1P"}, {l, "HIST"}, {hBackground, "HIST"}, {hphi, "HIST"}, {hSignal, "HIST"}, {h, "E1P"}, {hModel, "C"}};
			SaveCanvas(c, dat, TString::Format("media/root_files/phi_phi_reconstruction/"+ax+"_projection/"+ topo +"/mass_bound=%.4gMeV.root", mass_bound), "UPDATE");
		
			delete c;
			delete h;
			delete l;
			delete hModel;
			delete hphi;
			delete hSignal;
			delete hBackground;
	}
			delete file;
	}

	return 0;
}
