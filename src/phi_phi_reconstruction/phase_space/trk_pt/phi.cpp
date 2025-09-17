
#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFitResult.h"


#include "../../../../lib/computations.h"
#include "../../../../lib/admin_utils.h"
#include "../../../../lib/plotting_params.h"
#include "../../../../lib/roofithelper.h"

// !./bashing/run_file.sh phi_phi_reconstruction/phase_space/trk_pt/phi

using namespace RooFit;
int main() {

	SetPlotStyle();
	
	// set up transparant colours for use later
	
	Float_t kaon_mass = m_kaon_char;
	Float_t phi_mass = m_phi/1e3;
	RVecF mass = {5e5}; 
	mass /=1e3;

	RVecStr topology = {"40"};
	ROOT::EnableImplicitMT();
	
		
	for (TString topo : topology) {
		for (Float_t mass_bound : mass) {	
		
			TString file = TString::Format("data/phi_phi_reconstruction/uncut_SNR_" + topo + ".root");
			RDF df_df("tree", file);

			TCanvas* c = new TCanvas("c", "c", 600, 800);
			
			Int_t nbins = 500;		
			TH1* h = new TH1F("hist", "", nbins, 0, 1.5);
			h->SetTitle(TString::Format(";#phi p_t (GeV/c);Events [%.2g GeV/c^{2}]", h->GetBinWidth(1)));
			c->SetTopMargin(0.02);
			c->SetRightMargin(0.05);
			h->GetXaxis()->SetTitleOffset(1.2);
			gStyle->SetOptStat(1);

			df_df.Foreach([mass_bound, phi_mass, &h] (RVecLorCyl p) {
				if ((TMath::Abs(p[0].M() - phi_mass) < mass_bound)
					|| (TMath::Abs(p[1].M() - phi_mass) < mass_bound)
				) {
					h->Fill(p[0].Pt());
					h->Fill(p[1].Pt());
				} else if ( (TMath::Abs(p[2].M() - phi_mass) < mass_bound)
					|| (TMath::Abs(p[3].M() - phi_mass) < mass_bound)
				) {
					h->Fill(p[2].Pt());
					h->Fill(p[3].Pt());
					}}, {"phi_four_momentum"}
				);



			double x_min = h->GetXaxis()->GetXmin();
			double x_max = h->GetXaxis()->GetXmax();
			RooRealVar x("x", "mass", x_min, x_max);
			Float_t sf = 1e-3;	

			RooDataHist data("data", "Dataset from histogram", RooArgList(x), h);

			RooRealVar mean1("mean1", "mean of gauss1", 11e-3, 10e-3, 12e-3);
			RooRealVar sigma1("sigma1", "width of gauss1", 5e-3, 4e-3, 6e-3);
			RooGaussian* gauss1 = new RooGaussian("gauss1", "gauss1", x, mean1, sigma1);

			RooRealVar x_0("x_0", "threshold", -0.4, -0.5,  0.0);
			RooRealVar p0("p0", "power", 9, 4, 13);
			RooRealVar a0("a0", "a0", 1.0, -3.0, 3.0);
			RooRealVar b0("b0", "b0", -1.0, -3.0, 3.0);
		//	a0.setConstant(true);
		//	b0.setConstant(true);
			//x_0.setConstant(true);
			RooRealVar c0("c0", "c0", -14, -20, -5);
			RooAbsPdf* thresh = ThresholdBackgroundPdf("threshbkg", "threshbkg", x, x_0, p0, a0, b0, c0);

			RooAbsPdf* Background = thresh;

			Float_t h_int = h->Integral();
			RooRealVar nsig1("nsig1", "yield gauss1", h_int * 10e-3, 0, 100e-3*h_int);
			RooRealVar nbkg("nbkg", "background yield", h_int * 0.98, h_int * 0.85, h_int);

			RooAddPdf* model = new RooAddPdf("model", "sig+bg",
                                 RooArgList(*gauss1, *Background),
                                 RooArgList(nsig1, nbkg));

			model->fitTo(data);
			calc_chi2(x, data, model, nbins);
			
			Float_t binwidth = h->GetXaxis()->GetBinWidth(1);
			h->SetMarkerStyle(20);
			h->SetLineColor(kBlack);

			TH1* hModel = model->createHistogram("hModel", x, Binning(nbins));
			hModel->SetLineWidth(2);
			hModel->SetLineColor(kCyan);
			
			TH1* hGauss = gauss1->createHistogram("hgauss1", x, Binning(nbins));
			hGauss->Scale(nsig1.getVal());
			hGauss->SetFillColor(color_scheme[1]);
			
			TH1* hBackground = Background->createHistogram("hBackground", x, Binning(nbins));
			hBackground->Scale(nbkg.getVal());
			hBackground->SetFillColor(color_scheme[3]);
	
			std::cout << "mu = " << mean1.getVal() << std::endl;
			std::cout << "sig = " << sigma1.getVal() << std::endl;
			std::cout << "a0 = " << a0.getVal() << std::endl;
			std::cout << "b0 = " << b0.getVal() << std::endl;
			std::cout << "c0 = " << c0.getVal() << std::endl;
			std::cout << "p0 = " << p0.getVal() << std::endl;
			std::cout << "x_0 = " << x_0.getVal() << std::endl;


/*
			Float_t mu = mean1.getVal();
			Float_t sig = sigma1.getVal();
			std::cout << mu << "\n" << sig << std::endl;

			Float_t n = 1;
			Float_t int_data = hGauss->Integral(hGauss->FindBin(mu - n * sig), hGauss->FindBin(mu + n*sig));
			Float_t int_model = hModel->Integral(hModel->FindBin(mu -n * sig), hModel->FindBin(mu + n*sig));
			std::cout << "ratio = "<< int_data / int_model << std::endl;
*/

			TLegend* l = new TLegend(0.65, 0.55, 0.90, 0.70); // bottom left, top right
			l->AddEntry(h, "Data", "lep");
			l->AddEntry(hModel, "Model", "l");
			l->AddEntry(hGauss, "Signal", "f");
			l->AddEntry(hBackground, "Background", "f");
			l->SetTextSize(0.035);

			RVecDraw dat = {
				{h, "HIST"},
				{l, "HIST"},
				{hBackground, "HIST"},
				{hGauss, "HIST"},
				{h, "E1P"},
				{hModel, "C"}
			};
			SaveCanvas(c, dat, TString::Format("media/root_files/phi_phi_reconstruction/phase_space/mass_cuts/trk_pt/phi/"+ topo +"/mass_bound=%.4gMeV.root", mass_bound * 1e3), "UPDATE");
		
			delete c;
			delete h;
			delete l;
			delete hModel;
			delete hGauss;
			delete hBackground;
	}
	}

	return 0;
}
