#include <iostream>
#include <vector>

#include "TFile.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLegend.h"
#include "TH1.h"
#include "TString.h"

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

#include "../../lib/computations.h"
#include "../../lib/custom_definitions.h"
#include "../../lib/apply_cuts.cpp"
#include "../../lib/plotting_params.h"
#include "../../lib/admin_utils.h"

using namespace RooFit;

// !bashing/run_file.sh test/roofit



/** Looking at phase space distribution 
 * of individual and total pt
 */
enum class FitModelType {
	BackgroundOnly, 		// ARGUS
   	SignalOnly,				// Gaussian
   	SignalAndBackground		// ARGUS and gaussian
};

struct FitPlotCfg {
	Int_t nbins = 200;
	FitModelType modeltype = FitModelType::SignalAndBackground;
};
TH1F* FilterData(RDF& df, TH1F* hist, const char* branch, const Float_t val, const Float_t phi_mass);
void FitAndSave(TH1F* hist, const TString outname, const TString label, const struct FitPlotCfg&);


int main() {
	ROOT::EnableImplicitMT();

	const Float_t phi_mass = m_phi / 1e3; // GeV
	const char* branch = "trk_pt";

	const TString topologies[2] = {"20", "40"};

	// Mass cuts in GeV
	RVecF mass_cuts = {10, 15, 20, 25, 30};
	mass_cuts *= 1e-3;
	
	

	FitPlotCfg plot{200, FitModelType::SignalAndBackground};
	for (TString topo : topologies) {
		RDF df("tree", TString::Format("data/phi_phi_reconstruction/uncut_SNR_" + topo + ".root"));

		for (Float_t val : mass_cuts) {
			TH1F* hist = new TH1F("hist", "", plot.nbins, 0, 1000);
			TH1F *h= FilterData(df, hist, branch, val, phi_mass);
			
			const TString outname = TString::Format("media/root_files/phi_phi_reconstruction/phase_space/mass_cuts/%s/compare_topologies/fit/" + topo + "/mass_bound = %.3g MeV.root", branch, val*1e3);
			const TString label   = TString::Format("Topology "+ topo +", mass cut %.0f MeV", val*1e3);

			FitAndSave(h, outname, label, plot);
		}
	}

	return 0;
}


void FitAndSave(
		TH1F* hist, 
		const TString outname, 
		const TString label,
		const FitPlotCfg& cfg )
	{
	// Convert TH1F to RooDataHist

	double x_min = hist->GetXaxis()->GetXmin();
	double x_max = hist->GetXaxis()->GetXmax();
	RooRealVar m("m", "mass", x_min, x_max);
	

	RooDataHist data("data", "Dataset from histogram", RooArgList(m), hist);
	// Define the Gaussian signal
	RooRealVar mean("mean", "mean of gauss", 100, 50, 160);
	RooRealVar sigma("sigma", "width of gauss", 15, 0.01, 60);
	RooGaussian gauss("gauss", "signal gaussian", m, mean, sigma);

	// Define the Argus background
	RooRealVar m0("m0", "endpoint", x_max, x_min, x_max*2);
	RooRealVar p("p", "power", 0.7, 1e-5, 2);
	RooRealVar c("c", "curvature", 10, 0.1, 20);
	RooRealVar shift("shift", "pt offset", 10, 0, 50);
	ShiftedArgusPdf argus_class("argus", "shifted ARGUS", m, m0, c, p, shift);
	auto& argus = argus_class.getPdf();
	// Combine signal and background
	RooRealVar nsig("nsig", "signal yield", hist->Integral() * 0.5, 0, hist->Integral());
	RooRealVar nbkg("nbkg", "background yield", hist->Integral() * 0.5, 0, hist->Integral());
   
   RooAddPdf* model = nullptr; // pointer to assign later

	switch(cfg.modeltype) {
		case FitModelType::BackgroundOnly:
			model = new RooAddPdf("model", "bg", RooArgList(argus), RooArgList(nbkg));
			break;
		case FitModelType::SignalOnly:
			model = new RooAddPdf("model", "sig", RooArgList(gauss), RooArgList(nsig));
			break;
		case FitModelType::SignalAndBackground:
			model = new RooAddPdf("model", "sig+bg", RooArgList(gauss,argus), RooArgList(nsig,nbkg));
			break;
	}

	// Fit the model
	model->fitTo(data);
	double binWidth = hist->GetXaxis()->GetBinWidth(1);
	TH1* hModel = model->createHistogram("hModel", m, Binning(cfg.nbins));
	hist->SetMarkerStyle(20);

	TH1* hGauss = gauss.createHistogram("hGauss", m, Binning(cfg.nbins));
	hGauss->SetFillColor(kCyan);
	hGauss->Scale(nsig.getVal());
	
	TH1* hArgus = argus.createHistogram("hArgus", m, Binning(cfg.nbins));
	hArgus->Scale(nbkg.getVal());
	hArgus->SetFillColor(kRed);

	TCanvas* can = new TCanvas("can", "c", 600, 800);
	RVecDraw dat = {{hModel, "C"}, {hArgus, "HIST"}, {hGauss, "HIST"}, {hist, "EP"}};
	SaveCanvas(can, dat, outname, "RECREATE");

	delete can;
	delete hist;
	delete hModel;
	delete hGauss;
	delete hArgus;
}


TH1F* FilterData(
		RDF& df,
		TH1F* hist,
		const char* branch,
		const Float_t val,
		const Float_t phi_mass
		)
{
	// allocate dataset

	df.Filter([val, phi_mass](const RVecLorCyl& p) {
			return ((TMath::Abs(p[0].M() - phi_mass) < val) &&
					(TMath::Abs(p[1].M() - phi_mass) < val));
		}, {"phi_four_momentum"})
	  .Foreach([&](const RVecF& p4) {
			for (Float_t p : p4) {
				hist->Fill(p * 1e3);
			}
		}, {branch});

	return hist; // hand ownership to caller
}



