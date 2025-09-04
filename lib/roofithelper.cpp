#include "roofithelper.h"


void calc_chi2(RooRealVar& x, RooDataHist& h, RooAbsPdf* model, Int_t nbins){
	RooPlot* frame = x.frame();
	h.plotOn(frame);
	model->plotOn(frame);
	Double_t chi2PerNDF = frame->chiSquare();  // returns chi² / ndf
	Int_t nParams = model->getParameters(h)->getSize();
	Int_t ndf = nbins - nParams;
	Double_t chi2 = chi2PerNDF * ndf;
	std::cout << "Chi2 = " << chi2 << "\nChi2/ndf = " << chi2PerNDF << std::endl;
};

