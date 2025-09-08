#pragma once

#include "custom_definitions.h"


enum class FitModelType {
    Argus,            // ARGUS
    Gauss,            // Gaussian
    BiGaussAndArgus,  // BiGaussian and Argus Background
    GaussAndArgus,     // ARGUS and Gaussian
	GaussAndConst,
	GaussAndMassThreshold
};

// Struct for fit plot configuration
struct FitPlotCfg {
    Int_t nbins = 200;
    Float_t xmin = 0;
    Float_t xmax = 1000;
    FitModelType modeltype = FitModelType::GaussAndArgus;
};

void calc_chi2(RooRealVar& x, RooDataHist& h, RooAbsPdf* model, Int_t nbins);
