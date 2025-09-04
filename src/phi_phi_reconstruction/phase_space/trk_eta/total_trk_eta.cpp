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
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooAbsPdf.h"
#include "RooExtendPdf.h"

#include "../../../../lib/computations.h"
#include "../../../../lib/custom_definitions.h"
#include "../../../../lib/apply_cuts.cpp"
#include "../../../../lib/plotting_params.h"
#include "../../../../lib/admin_utils.h"
#include "../../../../lib/roofithelper.h"
using namespace RooFit;

// !bashing/run_file.sh phi_phi_reconstruction/phase_space/trk_eta/total_trk_eta



/** Looking at phase space distribution 
 * of individual and total eta
 */


TH1F* FilterData(RDF& df, TH1F* hist, const char* branch, const Float_t val, const Float_t phi_mass);
void FitAndSave(TH1F* hist, const TString outname, const TString label, const struct FitPlotCfg&);


int main() {
	ROOT::EnableImplicitMT();

	const Float_t phi_mass = m_phi / 1e3; // GeV
	const char* branch = "trk_eta";

	const RVecStr topologies = {"40"};

	// Mass cuts in GeV
	RVecF mass_cuts = {10, 15, 20, 25, 30};
	mass_cuts *= 1e-3;
	
	

	FitPlotCfg plot{40, -4, 4, FitModelType::Gauss};
	for (TString topo : topologies) {
		RDF df("tree", TString::Format("data/phi_phi_reconstruction/uncut_SNR_" + topo + ".root"));

		for (Float_t val : mass_cuts) {
			TH1F* hist = new TH1F("hist", "Total trk_eta; #eta (MeV/c);Events", plot.nbins, plot.xmin, plot.xmax);
			TH1F *h= FilterData(df, hist, branch, val, phi_mass);
			
			const TString outname = TString::Format("media/root_files/phi_phi_reconstruction/phase_space/mass_cuts/%s/" + topo + "/total_eta/mass_bound = %.3g MeV.root", branch, val*1e3);
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
	RooRealVar x_var("x_var", "eta", x_min, x_max);
	

	RooDataHist data("data", "Dataset from histogram", RooArgList(x_var), hist);
	// Define the Gaussian signal
	RooRealVar mean("mean", "mean of gauss", 0, -0.5, 0.5);
	RooRealVar sigma("sigma", "width of gauss", 0.1, 1e-3, 0.5);
	RooGaussian gauss("gauss", "signal gaussian", x_var, mean, sigma);
	
	RooRealVar mu1("mu1", "mean of first Gaussian", 0.5, 0.3, 0.7);
	RooRealVar sigma1("sigma1", "sigma of first Gaussian", 0.2, 0.02, 0.5);
	RooRealVar mu2("mu2", "mean of second Gaussian",1, 0.7, 1.2);
	RooRealVar sigma2("sigma2", "sigma of second Gaussian", 0.2, 0.02, 0.5);
	RooRealVar frac("frac", "fraction of first Gaussian", 0.5, 0., 1.);
	BiGaussianPdf bigauss_class("bigauss", "BiGaussian PDF", x_var, mu1, sigma1, mu2, sigma2, frac);
		
	// Define constant background
	RooPolynomial constbkg("constBkg", "Constant background", x_var);

	// Define the Argus background
	RooRealVar x_var0("x_var0", "endpoint", x_max, x_min, x_max*2);
	RooRealVar p("p", "power", 1, 1e-5, 6);
	RooRealVar c("c", "curvature", 5, 0.1, 20);
	RooRealVar shift("shift", "eta offset", 0.2, 0, 600);
	ShiftedArgusPdf argus_class("argus", "shifted ARGUS", x_var, x_var0, c, p, shift);
	
	// Combine signal and background
	RooRealVar nsig("nsig", "signal yield", hist->Integral() * 0.5, 0, hist->Integral());
	RooRealVar nbkg("nbkg", "background yield", hist->Integral() * 0.5, 0, hist->Integral());
   
	RooAddPdf* model = nullptr; // pointer to assign later
	RooAbsPdf* Signal = nullptr;
	RooAbsPdf* Background = nullptr;

	switch(cfg.modeltype) {
		case FitModelType::Argus:
			Background = &constbkg;
			Signal = &(argus_class.getPdf()); 
			model = new RooAddPdf("model", "sig+bg", RooArgList(*Signal, *Background), RooArgList(nsig, nbkg));
			break;
		case FitModelType::Gauss:
			Background = &constbkg;
			Signal = &gauss;
			model = new RooAddPdf("model", "sig_bg", RooArgList(*Signal, *Background), RooArgList(nsig, nbkg));
			break;
		case FitModelType::GaussAndConst:
			Background = &constbkg;
			Signal = &gauss;
			model = new RooAddPdf("model", "sig+bg", RooArgList(*Signal, *Background), RooArgList(nsig,nbkg));
			break;
		case FitModelType::GaussAndArgus:
			Background = &(argus_class.getPdf());
			Signal = &gauss;
			model = new RooAddPdf("model", "sig+bg", RooArgList(*Signal, *Background), RooArgList(nsig,nbkg));
			break;
		case FitModelType::BiGaussAndArgus:
			Signal = &(bigauss_class.getPdf());
			Background = &(argus_class.getPdf()); 
        model = new RooAddPdf("model", "sig+bg", RooArgList(*Signal, *Background), RooArgList(nsig, nbkg));
			break;
	}
	
	// Fit the model and get Chi squared value
	model->fitTo(data);
	calc_chi2(x_var, data, model, cfg.nbins);

	double binWidth = hist->GetXaxis()->GetBinWidth(1);
	TH1* hModel = model->createHistogram("hModel", x_var, Binning(cfg.nbins));
	hist->SetMarkerStyle(20);

	TH1* hSignal = Signal->createHistogram("hSignal", x_var, Binning(cfg.nbins));
	hSignal->SetFillColor(kCyan);
	hSignal->Scale(nsig.getVal());
	
	TH1* hBackground = Background->createHistogram("hBackground", x_var, Binning(cfg.nbins));
	hBackground->Scale(nbkg.getVal());
	hBackground->SetFillColor(kRed);

	TLegend* leg = new TLegend(0.7,0.45, 0.9, 0.7);
	leg->AddEntry(hist, "Data", "lep");
	leg->AddEntry(hModel, "Model Fit", "l");
	leg->AddEntry(hSignal, "Gaussian Signal", "f");
	leg->AddEntry(hBackground, "Background", "f");

	TCanvas* can = new TCanvas("can", "c", 600, 800);
	RVecDraw dat = {{hist, "EP"}, {hModel, "C"}, {hBackground, "HIST"}, {hSignal, "HIST"}, {hist, "EP"}, {leg, ""}};
	SaveCanvas(can, dat, outname, "RECREATE");

	delete can;
	delete hist;
	delete hModel;
	delete hSignal;
	delete hBackground;
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
		hist->Fill(Sum(p4));
	}, {branch});

	return hist; // hand ownership to caller
}



